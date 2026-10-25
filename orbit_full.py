"""Resumable whole-orbit 85 km run; native Master grid, 50-frame windows.

python orbit_full.py fit --workers 6 --iterations 32
python orbit_full.py assemble

Filled-hole pass (separate orbit85_full_filled output; sigma=1 matching):
python orbit_full.py fit --fill-small-gaps --workers 6 --iterations 32
python orbit_full.py assemble --fill-small-gaps
python verify_orbit_full.py --fill-small-gaps
python orbit_full.py score --fill-small-gaps --q-threshold 0.6

Shifted-window sensitivity pass (separate output):
python orbit_full.py fit --fill-small-gaps --start-offset 12
python orbit_full.py assemble --fill-small-gaps --start-offset 12
python combine_orbit_passes.py
"""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import orbit_pilot as p
from orbit_pilot import np, Path, json, csv, netCDF4, plt, arr, savejson
from scipy.ndimage import map_coordinates, gaussian_filter
from concurrent.futures import ProcessPoolExecutor, as_completed
import argparse, traceback, time

OUT=p.ROOT/(p.PREFIX+'_full')
PILOT=p.ROOT/(p.PREFIX+'_pilot')

def rows_read():
    rows=list(csv.DictReader((PILOT/'frame_inventory.csv').open()))
    for r in rows:
        for k in r:r[k]=None if r[k]=='' else (int(r[k]) if k in ('frame','x0','y0') else float(r[k]))
    return rows

def worker(start,iterations,output,fill_small_gaps=False,motion_smoothing=0):
    global OUT
    OUT=Path(output)
    p.OUT=OUT
    path=OUT/f'window_{start:04d}/result.json'
    if path.exists():
        cached=json.loads(path.read_text())
        if cached.get('fill_small_gaps',False)==fill_small_gaps and cached.get('motion_smoothing',0)==motion_smoothing and (cached['status']!='fit' or cached.get('model_components',2)==1 or cached.get('checkpoint_iteration',64)==iterations):return cached
    try:return p.run_window(start,rows_read(),iterations,fill_small_gaps=fill_small_gaps,motion_smoothing=motion_smoothing)
    except Exception as e:
        error=dict(start=start,end=start+49,status='error',reason=str(e),traceback=traceback.format_exc())
        savejson(OUT/f'error_{start:04d}.json',error);return error

def classify(r):
    if r['status']!='fit':return 'insufficient_evidence'
    modes=r['modes'];slow_ok=np.linalg.norm(modes[0]['velocity'])<.18 and modes[0]['pair_count']>=3
    if r.get('model_components',2)==1:return 'single_slow' if slow_ok else 'ambiguous_motion'
    if r['heldout_improvement_fraction']>.05 and all(m['pair_count']>=3 for m in modes):return 'two_motion'
    return 'single_slow' if slow_ok else 'ambiguous_motion'

def window_starts(n,offset=0):
    if n<50 or not 0<=offset<25:raise ValueError('Need >=50 frames and offset in [0,24]')
    return sorted(set([0,n-50]+list(range(offset,n-49,25))))


def fit(workers,iterations,fill_small_gaps=False,motion_smoothing=0,start_offset=0):
    OUT.mkdir(exist_ok=True);rows=rows_read();starts=window_starts(len(rows),start_offset)
    if (OUT/'run_config.json').exists():
        previous=json.loads((OUT/'run_config.json').read_text())
        if previous['starts']!=starts:raise ValueError('Output contains a different window layout; use a separate output folder')
    cover=np.zeros(len(rows),int)
    for s in starts:cover[s:s+50]+=1
    assert np.all(cover>0)
    savejson(OUT/'run_config.json',dict(source=str(p.L1),master=str(p.MASTER),l3c=str(p.L3),reference_altitude_km=85,window_frames=50,stride_frames=25,iterations=iterations,workers=workers,starts=starts,source_frames=len(rows),time_gap_policy='Keep actual timestamps; all source gaps are retained. No inserted frames.',training='Even local indices; odd local indices are held out.',model_selection='One supported slow mode is accepted. Two modes require >5% holdout improvement and >=3 supporting pairs each; otherwise use initial single-slow template when slow mode meets checks.',note='Model selection uses this holdout, so reported improvement is a selection diagnostic rather than untouched test performance.'))
    config=json.loads((OUT/'run_config.json').read_text())
    config.update(source_label=p.SOURCE_LABEL,orbit=p.ORBIT,input_preparation='See input projection provenance; L1F is spatially resampled with no temporal interpolation.' if p.SOURCE_LABEL!='L1R' else 'Existing L1R grid',start_offset=start_offset,endpoint_windows=[0,len(rows)-50],fill_small_gaps=fill_small_gaps,motion_smoothing=motion_smoothing,gap_policy='Enclosed 8-connected missing components <=16 pixels, bounding box <=5x5, valid geometry; harmonic interpolation. Original finite values unchanged. Masks saved per window.',support_policy='At least eight genuine training observations; errors use original observed pixels only.')
    savejson(OUT/'run_config.json',config)
    results=[];t=time.monotonic()
    with ProcessPoolExecutor(max_workers=workers) as ex:
        jobs={ex.submit(worker,s,iterations,str(OUT),fill_small_gaps,motion_smoothing):s for s in starts}
        for future in as_completed(jobs):
            r=future.result();r['selected_model']=classify(r);results.append(r)
            savejson(OUT/'progress.json',dict(completed=len(results),total=len(starts),elapsed_seconds=time.monotonic()-t,counts={k:sum(x['selected_model']==k for x in results) for k in sorted(set(x['selected_model'] for x in results))}))
            print('COMPLETED',len(results),'/',len(starts),'start',r['start'],r['selected_model'],flush=True)
    results.sort(key=lambda r:r['start']);savejson(OUT/'window_summary.json',results)
    if any(r['status']=='error' for r in results):raise RuntimeError('Window errors recorded; correct and rerun before assembly')

def sample(a,v,dt,x,y):
    coords=[y-v[1]*dt,x-v[0]*dt];ok=np.isfinite(a)
    w=map_coordinates(ok.astype('f4'),coords,order=1,mode='constant',cval=0)
    b=map_coordinates(np.nan_to_num(a),coords,order=1,mode='constant',cval=0)
    return np.where(w>.999,b/np.maximum(w,1e-9),np.nan)

def load_model(r):
    it=r.get('checkpoint_iteration',32)
    if r['selected_model']=='single_slow' and r.get('model_components',2)>1:it=0
    z=np.load(OUT/f"window_{r['start']:04d}/checkpoint_{it:03d}.npz")
    return {k:z[k] for k in z.files}

def offset_graph(results,models):
    ids=[i for i,r in enumerate(results) if r['selected_model']=='two_motion'];edges=[]
    for ia,i in enumerate(ids):
        a=results[i];ca=models[i]
        for j in ids[ia+1:]:
            b=results[j]
            if b['start']-a['start']>50:break
            cb=models[j];x0=max(a['bounds'][0],b['bounds'][0]);x1=min(a['bounds'][1],b['bounds'][1]);y0=max(a['bounds'][2],b['bounds'][2]);y1=min(a['bounds'][3],b['bounds'][3])
            if x1<=x0 or y1<=y0:continue
            yy,xx=np.mgrid[y0:y1,x0:x1];epoch=(a['reference_epoch']+b['reference_epoch'])/2
            aa=sample(ca['fast'],ca['fast_velocity'],epoch-a['reference_epoch'],xx-a['bounds'][0],yy-a['bounds'][2]);bb=sample(cb['fast'],cb['fast_velocity'],epoch-b['reference_epoch'],xx-b['bounds'][0],yy-b['bounds'][2])
            sa=sample(ca['fast_support'].astype(float),ca['fast_velocity'],epoch-a['reference_epoch'],xx-a['bounds'][0],yy-a['bounds'][2]);sb=sample(cb['fast_support'].astype(float),cb['fast_velocity'],epoch-b['reference_epoch'],xx-b['bounds'][0],yy-b['bounds'][2]);ok=np.isfinite(aa)&np.isfinite(bb)&(sa>=8)&(sb>=8)
            if ok.sum()<1000 or min(np.std(aa[ok]),np.std(bb[ok]))<.03:continue
            corr=float(np.corrcoef(aa[ok],bb[ok])[0,1]);delta=float(np.median(aa[ok]-bb[ok]));rmse=float(np.sqrt(np.mean((aa[ok]-bb[ok]-delta)**2)))
            edges.append(dict(i=i,j=j,starts=[a['start'],b['start']],pixels=int(ok.sum()),correlation=corr,median_difference_kR=delta,offset_removed_rmse_kR=rmse,accepted=corr>=.5))
    graph={i:set() for i in ids}
    for e in edges:
        if e['accepted']:graph[e['i']].add(e['j']);graph[e['j']].add(e['i'])
    unseen=set(ids);components=[];offsets={};component={}
    while unseen:
        seed=min(unseen);todo=[seed];group=set()
        while todo:
            i=todo.pop()
            if i in group:continue
            group.add(i);todo.extend(graph[i]-group)
        unseen-=group;g=sorted(group);pos={i:k for k,i in enumerate(g)};es=[e for e in edges if e['accepted'] and e['i'] in group]
        A=np.zeros((len(es)+1,len(g)));b=np.zeros(len(es)+1);w=np.ones(len(es)+1)
        for n,e in enumerate(es):A[n,pos[e['i']]]=-1;A[n,pos[e['j']]]=1;b[n]=e['median_difference_kR'];w[n]=np.sqrt(e['pixels'])*e['correlation']/max(e['offset_removed_rmse_kR'],.05)
        A[-1]=1/len(g);w[-1]=max(1,float(np.max(w)))
        solution=np.linalg.lstsq(A*w[:,None],b*w,rcond=None)[0]
        # Robust reweighting on redundant overlaps; zero-mean gauge per connected group.
        for _ in range(5):
            residual=A@solution-b;robust=np.minimum(1,.15/np.maximum(np.abs(residual),1e-8));robust[-1]=1;ww=w*np.sqrt(robust);solution=np.linalg.lstsq(A*ww[:,None],b*ww,rcond=None)[0]
        cid=len(components)
        for i,v in zip(g,solution):offsets[i]=float(v);component[i]=cid
        components.append(dict(id=cid,window_starts=[results[i]['start'] for i in g],offsets_kR=solution.tolist(),gauge='zero mean offset within this connected component'))
    for e in edges:
        if e['accepted']:e['remaining_median_difference_kR']=e['median_difference_kR']+offsets[e['i']]-offsets[e['j']]
    savejson(OUT/'offset_alignment.json',dict(components=components,edges=edges,note='Fast += offset; slow -= offset. Reconstruction is unchanged. Disconnected components have independent arbitrary absolute offsets.'))
    return offsets,component,components

def local_fast_score(observed,slow,fast,supported,sigma=3.,min_coverage=.95,min_energy=.04):
    """Walkthrough energy ratio; common original-observation mask, no edge reflection.

    Energy is the Gaussian-weighted sum with missing values set to zero, as in
    the walkthrough (not renormalized). q is only defined at valid centers.
    """
    ok=np.isfinite(observed)&np.isfinite(slow)&np.isfinite(fast)&supported
    coverage=gaussian_filter(ok.astype(float),sigma,mode='constant',cval=0)
    energy=gaussian_filter(np.where(ok,(observed-slow)**2,0.),sigma,mode='constant',cval=0)
    error=gaussian_filter(np.where(ok,(observed-slow-fast)**2,0.),sigma,mode='constant',cval=0)
    supported_local=ok&(coverage>=min_coverage)
    eligible=supported_local&(energy>=min_energy)&(energy>0)
    raw=np.full(observed.shape,np.nan)
    np.divide(error,energy,out=raw,where=eligible)
    raw[eligible]=1-raw[eligible]
    return raw,np.clip(raw,0,1),energy,error,coverage,supported_local


def score_orbit(threshold=.6,sigma=3.,min_coverage=.95,min_energy=.04,all_frames=False):
    """Score dominant two-motion models at exact withheld L3C source epochs.

    Do not interpolate scores in time or silently substitute a neighboring
    model. The score is for the dominant model, not the blended component map.
    """
    results=json.loads((OUT/'window_summary.json').read_text());rows=rows_read()
    with netCDF4.Dataset(OUT/'orbit85_separation.nc') as d:
        status=np.asarray(d['Status'][:]);choice=np.asarray(d['Dominant_Window_Start'][:]);epoch=arr(d['Pixel_Epoch'][:])
        cloud=arr(d['Existing_MLCloud'][:]);fastmap=arr(d['Fast_Contrast'][:])
    alignment=json.loads((OUT/'offset_alignment.json').read_text())
    offsets={s:v for c in alignment['components'] for s,v in zip(c['window_starts'],c['offsets_kR'])}
    shape=status.shape;maps={k:np.full(shape,np.nan,'f4') for k in ['Q_Raw','Q_Score','Q_Residual_Energy','Q_Error_Energy','Q_Local_Coverage']}
    # Reasons: 0 outside, 1 not two-motion, 2 training timestamp, 3 unsupported,
    # 4 weak residual, 5 eligible. The binary mask keeps unavailable separate.
    reason=np.where(status==0,0,np.where(status==3,2,1)).astype('i1')
    observation=np.full(shape,np.nan,'f4');scored_windows=0
    with netCDF4.Dataset(p.L1) as source:
        for r in results:
            if r['selected_model']!='two_motion':continue
            c=load_model(r);x0,x1,y0,y1=r['bounds'];sl=np.s_[y0:y1,x0:x1]
            target=(choice[sl]==r['start'])&(status[sl]==3);e=epoch[sl]
            yy,xx=np.indices(e.shape,dtype=float);offset=offsets[r['start']]
            for frame in (range(r['start'],r['end']+1) if all_frames else r['heldout_frames']):
                select=target&(e==rows[frame]['epoch'])
                if not select.any():continue
                observed=np.full(e.shape,np.nan,'f4');row=rows[frame]
                ax,bx=max(x0,row['x0']),min(x1,row['x0']+300);ay,by=max(y0,row['y0']),min(y1,row['y0']+300)
                if ax<bx and ay<by:
                    a=arr(source['Radiance'][frame])[::-1];a[(a<0)|(a>65535)]=np.nan
                    observed[ay-y0:by-y0,ax-x0:bx-x0]=a[ay-row['y0']:by-row['y0'],ax-row['x0']:bx-row['x0']]
                dt=rows[frame]['epoch']-r['reference_epoch']
                slow=sample(c['slow'],c['slow_velocity'],dt,xx,yy)-offset
                fast=sample(c['fast'],c['fast_velocity'],dt,xx,yy)+offset
                sc=sample(c['slow_support'].astype(float),c['slow_velocity'],dt,xx,yy)
                fc=sample(c['fast_support'].astype(float),c['fast_velocity'],dt,xx,yy)
                raw,q,energy,error,coverage,supported=local_fast_score(observed,slow,fast,(sc>=8)&(fc>=8),sigma,min_coverage,min_energy)
                observation[sl][select]=observed[select]
                for name,a in zip(maps,[raw,q,energy,error,coverage]):maps[name][sl][select]=a[select]
                why=reason[sl];why[select]=3;why[select&supported]=4;why[select&np.isfinite(q)]=5
            scored_windows+=1
            print('Q scored window',r['start'],flush=True)
    q=maps['Q_Score'];mask=np.full(shape,-1,'i1');eligible=np.isfinite(q);mask[eligible]=(q[eligible].astype('f8')>threshold).astype('i1')
    assert np.array_equal(eligible,reason==5)
    settings=dict(threshold=threshold,gaussian_sigma_native_cells=sigma,min_local_coverage=min_coverage,min_weighted_residual_energy_kR2=min_energy,min_training_support=8,strict_threshold=True,definition='clip(1 - G[(I-S-F)^2 valid] / G[(I-S)^2 valid],0,1)',evaluation='Prepared input observations only, exact L3C Pixel_Epoch, odd local indices of dominant accepted two-motion window. No temporal interpolation. This validates templates within windows, not globally independent predictions.',gauge='Uses offset-aligned dominant window templates. The denominator and fast sign depend on the brightness-split gauge; q is not gauge-invariant. The total reconstruction error is invariant.',interpretation='Candidate fast-feature mask, not cloud probability. Both signs of fast contrast can qualify. Single-slow and unevaluable pixels remain unavailable.',defaults='Exploratory, not calibrated. Sigma=3 native cells is nominally 6 km, not identical to three pixels on the old latitude/longitude grid.')
    if all_frames:settings['evaluation']='All prepared input observations at exact L3C Pixel_Epoch, using dominant accepted two-motion model. Includes training frames: descriptive fit diagnostic, not withheld validation. No temporal score interpolation.'
    prefix='Q_AllFrames_' if all_frames else 'Q_'
    with netCDF4.Dataset(OUT/'orbit85_separation.nc','a') as d:
        d.setncattr('q_all_frames_settings_json' if all_frames else 'q_settings_json',json.dumps(settings))
        if all_frames:maps['Q_Observation']=observation
        for name,a in maps.items():
            name=name.replace('Q_',prefix,1)
            v=d[name] if name in d.variables else d.createVariable(name,'f4',('cross_track','along_track'),zlib=True,fill_value=np.nan)
            v[:]=a;v.units='kRayleigh2' if name.endswith(('Residual_Energy','Error_Energy')) else '1'
            v.model_provenance='Dominant_Window_Start; offset-aligned dominant model, not blended map; exact Pixel_Epoch'
        for name,a in [('Q_Candidate_Mask',mask),('Q_Status',reason)]:
            name=name.replace('Q_',prefix,1)
            v=d[name] if name in d.variables else d.createVariable(name,'i1',('cross_track','along_track'),zlib=True)
            v[:]=a
        d[prefix+'Candidate_Mask'].flag_values=np.array([-1,0,1],dtype='i1');d[prefix+'Candidate_Mask'].flag_meanings='unavailable below_or_equal_threshold candidate_fast_feature'
        d[prefix+'Status'].flag_values=np.arange(6,dtype='i1');d[prefix+'Status'].flag_meanings='outside_l3c no_selected_two_motion training_timestamp insufficient_local_support weak_residual_energy eligible'
        if all_frames:d[prefix+'Observation'].units='kRayleigh'
    summary=dict(settings=settings,two_motion_windows=scored_windows,eligible_pixels=int(eligible.sum()),candidate_pixels=int((mask==1).sum()),status_counts={str(i):int((reason==i).sum()) for i in range(6)})
    savejson(OUT/('q_all_frames_summary.json' if all_frames else 'q_summary.json'),summary)
    if all_frames:
        assert not np.any(reason==2), 'All selected two-motion timestamps must have been evaluated'
        plot_q_swaths(q,observation,threshold)
        (OUT/'Q_DISPLAY.md').write_text('# Radiance and q display\n\nq_radiance_swaths.png includes L3C radiance, existing MLCloud, slow and fast components, prepared input radiance at scoring timestamps, all-frame q, prepared input radiance masked by q > '+str(threshold)+', and model availability. Four q_radiance_segment_*.png files provide enlarged views.\n\nNarrow stripes in the withheld-only q_score_whole_orbit.png represent alternate timestamps used for training. The new display evaluates every timestamp directly, including training frames; it is a descriptive diagnostic, not independent validation. No gaps are cosmetically filled or scores interpolated. Q_AllFrames_* variables are separate from the original withheld-only Q_* variables in orbit85_separation.nc. Remaining gaps reflect single-slow models, insufficient support or weak signal.\n\nThe mask shows actual prepared input radiance, not a binary color block or the coadded L3C values. It uses the same radiance scale as the unmasked radiance panels. q uses the dominant two-motion window; slow and fast orbit panels are blended across compatible windows.\n',encoding='utf-8')
        return summary
    from matplotlib.colors import ListedColormap,BoundaryNorm
    fig,axes=plt.subplots(4,1,figsize=(20,10),sharex=True,layout='constrained')
    fastlim=max(.1,float(np.nanpercentile(np.abs(fastmap[status==3]),98)))
    for ax,a,title,cm,lo,hi in zip(axes[:3],[cloud,fastmap,q],['Existing MLCloud input (%)','Fast contrast (blended orbit map; signed kR)','Local q | dominant model, exact withheld epochs only'],['viridis','RdBu_r','viridis'],[0,-fastlim,0],[100,fastlim,1]):
        im=ax.imshow(a,origin='lower',aspect='auto',cmap=cm,vmin=lo,vmax=hi);ax.set_title(title);ax.set_facecolor('#ccc');fig.colorbar(im,ax=ax,pad=.01)
    im=axes[3].imshow(mask,origin='lower',aspect='auto',cmap=ListedColormap(['#ccc','#f5f5f5','#d16a1d']),norm=BoundaryNorm([-1.5,-.5,.5,1.5],3))
    cb=fig.colorbar(im,ax=axes[3],ticks=[-1,0,1],pad=.01);cb.ax.set_yticklabels(['Unavailable',f'q <= {threshold}',f'q > {threshold}'])
    axes[3].set(title='Candidate fast-feature mask (not cloud probability)',xlabel='Master along-track cell');fig.suptitle(f'Orbit {p.ORBIT:05d} | 85 km | local explained-energy score')
    fig.savefig(OUT/'q_score_whole_orbit.png',dpi=150);plt.close(fig)
    (OUT/'Q_README.md').write_text('# Local fast-component score\n\n'+json.dumps(summary,indent=2)+'\n\nUnavailable scores are NaN and masks are -1, not zero. The source timestamp must be held out from the dominant model. Narrow missing stripes are therefore expected. Energies use identical Gaussian weights on genuine observations with valid predictions and >=8 training observations per component. Values outside the image count as missing, unlike reflected-edge filtering. The >0.6 rule is strict. q compares the final slow component with that same slow component plus fast; this differs from comparing two separately fitted models.\n\nRun `python orbit_full.py score --fill-small-gaps --q-threshold 0.6` to refresh without refitting. Settings: --q-sigma, --q-min-coverage, --q-min-energy. Assembly also refreshes q automatically.\n',encoding='utf-8')
    return summary


def plot_q_swaths(q,observation,threshold):
    """Walkthrough-style masked observed radiance, alongside existing swaths."""
    from matplotlib.colors import ListedColormap,BoundaryNorm
    q=np.asarray(q,dtype='f8') # compare saved f4 values against the exact Python threshold
    with netCDF4.Dataset(OUT/'orbit85_separation.nc') as d:
        rad=arr(d['L3C_Radiance'][:]);cloud=arr(d['Existing_MLCloud'][:]);slow=arr(d['Slow_Component'][:]);fast=arr(d['Fast_Contrast'][:]);status=np.asarray(d['Status'][:])
    lim=np.nanpercentile(rad,[2,98]);flim=max(.1,float(np.nanpercentile(np.abs(fast[status==3]),98)))
    panels=[(rad,'Existing L3C radiance (coadded input)','gray',lim,'kR'),(cloud,'Existing MLCloud input','viridis',[0,100],'%'),(slow,'Slow component at L3C pixel times','gray',lim,'kR'),(fast,'Fast contrast (signed; blended orbit map)','RdBu_r',[-flim,flim],'kR'),(observation,'Prepared input radiance at scoring times | selected two-motion regions','gray',lim,'kR'),(q,'Local q | all-frame diagnostic, includes training observations','viridis',[0,1],'q'),(np.where(q>threshold,observation,np.nan),f'Prepared input radiance where q > {threshold} | candidate fast-pattern mask','gray',lim,'kR')]
    for name,x0,x1 in [('q_radiance_swaths',0,rad.shape[1])]+[(f'q_radiance_segment_{i+1}',a,b) for i,(a,b) in enumerate(list(zip(np.linspace(0,rad.shape[1],5,dtype=int)[:-1],np.linspace(0,rad.shape[1],5,dtype=int)[1:])))]:
        fig,axes=plt.subplots(8,1,figsize=(20,18),sharex=True,layout='constrained')
        for ax,(a,title,cm,limits,unit) in zip(axes,panels):
            im=ax.imshow(a[:,x0:x1],origin='lower',aspect='auto',extent=[x0,x1,0,a.shape[0]],cmap=cm,vmin=limits[0],vmax=limits[1],interpolation='nearest');ax.set(title=title,ylabel='Cross-track cell');ax.set_facecolor('#ccc');fig.colorbar(im,ax=ax,pad=.01,label=unit)
        im=axes[-1].imshow(status[:,x0:x1],origin='lower',aspect='auto',extent=[x0,x1,0,status.shape[0]],cmap=ListedColormap(['#eee','#aaa','#3b9c75','#6c63bd']),norm=BoundaryNorm([-.5,.5,1.5,2.5,3.5],4),interpolation='nearest')
        cb=fig.colorbar(im,ax=axes[-1],ticks=[0,1,2,3],pad=.01);cb.ax.set_yticklabels(['Outside','Insufficient','Single slow','Two motions']);axes[-1].set(title='Model availability',xlabel='Master along-track cell',ylabel='Cross-track cell')
        fig.suptitle(f'Orbit {p.ORBIT:05d} | 85 km | radiance and local q\nAll-frame q is a fit diagnostic, not independent validation; gray means unavailable or masked',fontsize=15)
        fig.savefig(OUT/(name+'.png'),dpi=120);plt.close(fig)


def assemble():
    results=json.loads((OUT/'window_summary.json').read_text());rows=rows_read();config=json.loads((OUT/'run_config.json').read_text())
    assert len(results)==len(config['starts']) and all(r['status']!='error' for r in results)
    models={i:load_model(r) for i,r in enumerate(results) if r['selected_model'] in ('two_motion','single_slow')}
    offsets,components,groups=offset_graph(results,models)
    z=np.load(PILOT/'l3c_reference.npz');rad=z['radiance'];cloud=z['mlcloud'];epoch=z['epoch'];shape=rad.shape;valid=np.isfinite(rad)&np.isfinite(epoch)
    # Status 0 outside L3C; 1 insufficient evidence/support; 2 single slow; 3 two motions.
    status=np.where(valid,1,0).astype('i1');choice=np.full(shape,-1,'i2');component_map=np.full(shape,-1,'i2');best=np.zeros(shape,'f4')
    slow=np.full(shape,np.nan,'f4');fast=slow.copy();rawfast=slow.copy();sumslow=np.zeros(shape,'f8');sumfast=sumslow.copy();sumraw=sumslow.copy();weights=sumslow.copy();count=np.zeros(shape,'i2');support=np.zeros(shape,'f4')
    preds={}
    for i,c in models.items():
        r=results[i];x0,x1,y0,y1=r['bounds'];e=epoch[y0:y1,x0:x1];dt=e-r['reference_epoch'];safe=np.where(np.isfinite(dt),dt,0);yy,xx=np.indices(e.shape,dtype=float)
        a=sample(c['slow'],c['slow_velocity'],safe,xx,yy);sc=sample(c['slow_support'].astype(float),c['slow_velocity'],safe,xx,yy)
        single=r['selected_model']=='single_slow';f=np.zeros(a.shape) if single else sample(c['fast'],c['fast_velocity'],safe,xx,yy);fc=sc if single else sample(c['fast_support'].astype(float),c['fast_velocity'],safe,xx,yy)
        keep=valid[y0:y1,x0:x1]&(e>=rows[r['start']]['epoch'])&(e<=rows[r['end']]['epoch'])&(sc>=8)&(fc>=8)&np.isfinite(a)&np.isfinite(f)
        # Soft overlap weights avoid central-window gaps; never extrapolate beyond observed times.
        w=np.exp(-.5*(safe/16.)**2)*np.minimum(sc,fc)/25;w=np.where(keep,w,0)
        offset=offsets.get(i,0.);preds[i]=(a-offset,f+offset,f,w,np.minimum(sc,fc))
        better=w>best[y0:y1,x0:x1];best[y0:y1,x0:x1][better]=w[better];choice[y0:y1,x0:x1][better]=r['start'];status[y0:y1,x0:x1][better]=2 if single else 3;component_map[y0:y1,x0:x1][better]=components.get(i,-1)
    # Blend only compatible models in the same offset-connected component.
    for i,(a,f,raw,w,sup) in preds.items():
        r=results[i];x0,x1,y0,y1=r['bounds'];sl=np.s_[y0:y1,x0:x1];single=r['selected_model']=='single_slow'
        compatible=(status[sl]==2) if single else ((status[sl]==3)&(component_map[sl]==components[i]))
        ww=np.where(compatible,w,0);good=ww>0
        sumslow[sl]+=np.where(good,a,0)*ww;sumfast[sl]+=np.where(good,f,0)*ww;sumraw[sl]+=np.where(good,raw,0)*ww;weights[sl]+=ww;count[sl]+=good;support[sl]+=np.where(good,sup,0)*ww
    good=weights>0;slow[good]=sumslow[good]/weights[good];fast[good]=sumfast[good]/weights[good];rawfast[good]=sumraw[good]/weights[good];support[good]/=weights[good]
    assert np.array_equal(good,status>=2)
    # Fast=0 in selected single-layer model is a model convention, not measured cloud absence.
    full_epochs=np.array([r['epoch'] for r in rows]);coverage=np.zeros(len(rows),int)
    for r in results:coverage[r['start']:r['end']+1]+=1
    np.testing.assert_array_less(0,coverage)
    e=epoch[valid];j=np.searchsorted(full_epochs,e);j=np.clip(j,1,len(full_epochs)-1);diff=np.minimum(np.abs(e-full_epochs[j]),np.abs(e-full_epochs[j-1]));assert np.max(diff)<1e-6
    la,lo=p.master()
    with netCDF4.Dataset(OUT/'orbit85_separation.nc','w') as d:
        d.createDimension('cross_track',shape[0]);d.createDimension('along_track',shape[1]);d.reference_altitude_km=85.;d.master_file=str(p.MASTER);d.source_file=str(p.L1)
        d.source_label=p.SOURCE_LABEL;d.orbit=p.ORBIT
        d.description='Whole available orbit processed in overlapping windows. Motion-derived radiance evaluated at existing L3C Pixel_Epoch; not an official L3C reprocessing or reproduction of its upstream coadd. Single-layer fast=0 is a model convention. Missing estimates are NaN. Absolute brightness split is unidentified; each alignment component has its own gauge.'
        d.small_gap_filling=int(config.get('fill_small_gaps',False));d.motion_smoothing_sigma_pixels=config.get('motion_smoothing',0)
        d.interpolation_provenance='Per-window interpolation_provenance.npz masks. Support and error metrics use genuine observations.'
        for name,data,units in [('Latitude',la,'degrees_north'),('Longitude',lo,'degrees_east'),('L3C_Radiance',rad,'kRayleigh'),('Existing_MLCloud',cloud,'percent'),('Slow_Component',slow,'kRayleigh'),('Fast_Contrast',fast,'kRayleigh'),('Fast_Before_Offset_Alignment',rawfast,'kRayleigh'),('Pixel_Epoch',epoch,'source Epoch seconds'),('Mean_Training_Support',support,'frames')]:
            v=d.createVariable(name,'f8' if name=='Pixel_Epoch' else 'f4',('cross_track','along_track'),zlib=True,fill_value=np.nan);v[:]=data;v.units=units
        for name,data in [('Status',status),('Dominant_Window_Start',choice),('Offset_Component',component_map),('Contributing_Windows',count)]:
            v=d.createVariable(name,'i2',('cross_track','along_track'),zlib=True);v[:]=data
        d['Status'].flag_values=np.array([0,1,2,3],dtype='i2');d['Status'].flag_meanings='outside_l3c insufficient_evidence_or_support single_slow two_motion'
    counts={k:sum(r['selected_model']==k for r in results) for k in sorted(set(r['selected_model'] for r in results))};pix={str(k):int((status==k).sum()) for k in range(4)}
    summary=dict(frames_considered=len(rows),windows=len(results),window_classifications=counts,l3c_valid_pixels=int(valid.sum()),pixels_by_status=pix,supported_fraction_of_l3c=float(good.sum()/valid.sum()),offset_connected_components=len(groups),maximum_l3c_epoch_match_error_seconds=float(np.max(diff)),all_frames_in_at_least_one_window=bool((coverage>0).all()),two_motion_holdout_improvement_percent=[100*r['heldout_improvement_fraction'] for r in results if r['selected_model']=='two_motion'])
    summary.update(fill_small_gaps=config.get('fill_small_gaps',False),motion_smoothing=config.get('motion_smoothing',0),interpolated_window_frame_samples=sum(r.get('interpolated_frame_pixels',0) for r in results))
    savejson(OUT/'summary.json',summary)
    plot(rad,cloud,slow,fast,status,results,summary)
    report(results,summary,groups)
    print(json.dumps(summary,indent=2),flush=True)

def plot(rad,cloud,slow,fast,status,results,summary):
    from matplotlib.colors import ListedColormap,BoundaryNorm
    fig,axes=plt.subplots(5,1,figsize=(20,12),sharex=True,layout='constrained');lim=np.nanpercentile(rad,[2,98]);flim=max(.1,float(np.nanpercentile(np.abs(fast[status==3]),98)))
    panels=[(rad,'Existing L3C radiance','gray',lim,'kR'),(cloud,'Existing MLCloud map','viridis',[0,100],'%'),(slow,'Estimated slow component at L3C pixel times','gray',lim,'kR'),(fast,'Fast contrast | single-layer regions use modeled zero','RdBu_r',[-flim,flim],'kR')]
    for ax,(data,title,cm,limits,unit) in zip(axes,panels):
        im=ax.imshow(data,origin='lower',aspect='auto',cmap=cm,vmin=limits[0],vmax=limits[1]);ax.set(title=title,ylabel='Cross-track cell');ax.set_facecolor('#ddd');fig.colorbar(im,ax=ax,pad=.01,label=unit)
    colors=['#eeeeee','#b6b6b6','#3b9c75','#6c63bd'];im=axes[4].imshow(status,origin='lower',aspect='auto',cmap=ListedColormap(colors),norm=BoundaryNorm([-.5,.5,1.5,2.5,3.5],4));cb=fig.colorbar(im,ax=axes[4],ticks=[0,1,2,3],pad=.01);cb.ax.set_yticklabels(['Outside','Insufficient','Single slow','Two motions']);axes[4].set(title='Processing outcome | every source frame was included',ylabel='Cross-track cell',xlabel='Master along-track cell (nominal 2 km/cell; vertical axis stretched)')
    fig.suptitle(f"Orbit {p.ORBIT:05d} | 85 km | {summary['frames_considered']:,} frames / {summary['windows']} windows | supported map coverage {summary['supported_fraction_of_l3c']:.1%}",fontsize=16);fig.savefig(OUT/'l3c_cloud_whole_orbit.png',dpi=160);plt.close(fig)
    for n,(x0,x1) in enumerate(list(zip(np.linspace(0,rad.shape[1],5,dtype=int)[:-1],np.linspace(0,rad.shape[1],5,dtype=int)[1:]))):
        fig,axes=plt.subplots(3,1,figsize=(16,8),sharex=True,layout='constrained')
        for ax,data,title,cm,limits in zip(axes,[rad[:,x0:x1],fast[:,x0:x1],status[:,x0:x1]],['L3C radiance','Aligned fast contrast (modeled zero for single slow)','Outcome: gray insufficient, green single slow, purple two motions'],['gray','RdBu_r',ListedColormap(colors)],[lim,[-flim,flim],[0,3]]):
            im=ax.imshow(data,origin='lower',aspect='auto',extent=[x0,x1,0,rad.shape[0]],cmap=cm,vmin=limits[0],vmax=limits[1]);ax.set(title=title,ylabel='Cross-track cell');ax.set_facecolor('#ddd');fig.colorbar(im,ax=ax,pad=.01)
        axes[-1].set_xlabel('Master along-track cell');fig.suptitle(f'Orbit {p.ORBIT} • segment {n+1} of 4 • 85 km');fig.savefig(OUT/f'segment_{n+1}.png',dpi=150);plt.close(fig)

def report(results,s,groups):
    lines=['# Whole-orbit separation at 85 km','',f"Processed all {s['frames_considered']} source frames through {s['windows']} overlapping windows. Supported estimates cover {s['supported_fraction_of_l3c']:.1%} of valid L3C pixels. The remainder were processed but lack accepted motion fits or adequate spatial/temporal support.",'',f"Window classifications: {s['window_classifications']}",'','## Interpretation','','- One slow motion is a valid outcome. Fast=0 there means the selected model contains no fast layer; it does not prove cloud absence.','- Two-motion estimates require >5% held-out improvement and evidence from at least three frame pairs per mode. Motion labels alone do not establish physical cloud/airglow identity.','- Every 50-frame window uses even local frames for motion estimation and fitting and odd frames for model selection. Adjacent frames are correlated; these are exploratory selection diagnostics, not independent validation.','- Fits use 32 iterations for this first orbit pass; the pilot showed small additional error reduction from 32 to 64. Saved histories permit targeted refinement.','- Actual timestamps, missing radiance, and frame gaps are preserved. L1F runs use a floating-point projection onto the Master grid; all diagnostics refer to those resampled observations. No further spatial regridding is performed.','- Maps are evaluated at L3C Pixel_Epoch. Existing L3C is a coadded product, so these are aligned component estimates, not an official reprocessing of the L3C coadd.','- Overlapping windows are weighted toward their reference times and require at least eight training-frame contributions. Different model classes and disconnected offset groups are not blended together.',f'- Brightness offsets were aligned in {len(groups)} connected groups using same-time overlaps. Fast gains an offset and slow loses it, preserving their sum. Absolute offsets remain arbitrary between disconnected groups.','- Gray missing regions are insufficient evidence/support, not cloud-free areas.','','## Outputs','','- `l3c_cloud_whole_orbit.png`: orbit-wide L3C, cloud score, slow/fast estimates, and coverage status.','- `segment_1.png` through `segment_4.png`: enlarged views.','- `orbit85_separation.nc`: geolocated estimates, status, timestamps and window provenance.','- `window_summary.json`: all window outcomes and held-out histories.','- `offset_alignment.json`: overlap evidence, solved offsets and connected groups.','- `summary.json`: frame and pixel coverage counts.','','## Window results','','| Start | End | Selected model | Holdout RMSE, final candidate | Improvement vs slow baseline |','|---:|---:|---|---:|---:|']
    for r in results:
        h=r.get('history',[]);rmse=f"{h[-1]['heldout_rmse_kR']:.4f}" if h else '—';imp=f"{r['heldout_improvement_fraction']:.1%}" if h else '—';lines.append(f"| {r['start']} | {r['end']} | {r['selected_model']} | {rmse} | {imp} |")
    lines+=['','## Acceptance criteria and checks','','A supported motion mode requires at least 15 matching candidates within 0.085 native cells/s of its peak. A single mode is labeled slow only below 0.18 cells/s with at least three supporting frame pairs. These thresholds are heuristic. A well-supported lone faster mode is classified as ambiguous motion, rather than assuming it is cloud-free.','', 'Run `python verify_orbit_full.py` for numerical and product-integrity checks; the completed results are in `verification.json`.','','Reproduce: `python orbit_full.py fit --workers 6 --iterations 32`, then `python orbit_full.py assemble`. Completed windows with matching iteration settings are reused. Source files and pilot outputs are unchanged.']
    if s.get('fill_small_gaps'):
        lines=[line.replace('Actual timestamps, missing radiance, and frame gaps are preserved.','Actual timestamps and frame gaps are preserved. Small enclosed radiance holes are interpolated for fitting; originals are retained for evaluation.') for line in lines]
        lines+=['','## Small-hole processing','','Only enclosed missing components of at most 16 pixels and at most 5x5 cells, with valid geometry, are filled by harmonic interpolation. Original finite radiance is unchanged. Larger holes and footprint edges remain missing. Motion matching uses Gaussian smoothing with sigma '+str(s['motion_smoothing'])+' pixels; templates use the filled, unsmoothed radiance. Evaluation and minimum support use genuine observations. Each window saves an interpolation mask.','', 'Reproduce this run: `python orbit_full.py fit --fill-small-gaps --workers 6 --iterations 32`, then `python orbit_full.py assemble --fill-small-gaps`. Verify with `python verify_orbit_full.py --fill-small-gaps`.']
    config=json.loads((OUT/'run_config.json').read_text())
    if config.get('start_offset'):
        offset=config['start_offset']
        lines=[line.replace('--fill-small-gaps',f'--fill-small-gaps --start-offset {offset}') for line in lines]
        lines+=['',f'Regular window starts are shifted by {offset} frames. Endpoint starts 0 and {config["source_frames"]-50} are retained; they are shared with the original pass and are not independent comparisons.']
    if p.PREFIX!='orbit85':
        lines=[line for line in lines if not line.startswith('Reproduce') and not line.startswith('Run `python verify_orbit_full.py')]
    lines += ['',f'Input: {p.SOURCE_LABEL}; orbit {p.ORBIT}.',f"Reproduce: `python run_orbit.py --orbit {p.ORBIT} --source {'l1r' if p.SOURCE_LABEL=='L1R' else 'l1f'}`." if p.PREFIX!='orbit85' else '']
    (OUT/'README.md').write_text('\n'.join(lines),encoding='utf-8')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['fit','assemble','score']);ap.add_argument('--workers',type=int,default=6);ap.add_argument('--iterations',type=int,default=32);ap.add_argument('--fill-small-gaps',action='store_true');ap.add_argument('--motion-smoothing',type=float,default=None)
    ap.add_argument('--start-offset',type=int,default=0,help='Shift regular starts by 0..24 frames; retain endpoint windows; use separate output suffix')
    ap.add_argument('--q-threshold',type=float,default=.6);ap.add_argument('--q-sigma',type=float,default=3.);ap.add_argument('--q-min-coverage',type=float,default=.95);ap.add_argument('--q-min-energy',type=float,default=.04);args=ap.parse_args()
    if not 0<=args.start_offset<25:ap.error('--start-offset must be in [0,24]')
    if not 0<=args.q_threshold<=1:ap.error('--q-threshold must be in [0,1]')
    if not np.isfinite(args.q_sigma) or args.q_sigma<=0:ap.error('--q-sigma must be positive and finite')
    if not 0<args.q_min_coverage<=1:ap.error('--q-min-coverage must be in (0,1]')
    if not np.isfinite(args.q_min_energy) or args.q_min_energy<=0:ap.error('--q-min-energy must be positive and finite')
    if args.fill_small_gaps:OUT=p.ROOT/'orbit85_full_filled'
    if args.start_offset:OUT=OUT.with_name(OUT.name+f'_shift{args.start_offset:02d}')
    smoothing=(1 if args.fill_small_gaps else 0) if args.motion_smoothing is None else args.motion_smoothing
    if smoothing<0:ap.error('--motion-smoothing must be nonnegative')
    if not args.fill_small_gaps and smoothing:ap.error('Use --fill-small-gaps for the separate preprocessing output')
    if args.mode=='fit':fit(args.workers,args.iterations,args.fill_small_gaps,smoothing,args.start_offset)
    else:
        if args.mode=='assemble':assemble()
        print(json.dumps(score_orbit(args.q_threshold,args.q_sigma,args.q_min_coverage,args.q_min_energy),indent=2))
        print(json.dumps(score_orbit(args.q_threshold,args.q_sigma,args.q_min_coverage,args.q_min_energy,all_frames=True),indent=2))
