"""Compare shifted windows without selecting the maximum q or hiding conflicts."""
import orbit_full as f
from orbit_full import np,json,netCDF4,arr,plt,Path
from matplotlib.colors import ListedColormap,BoundaryNorm
import argparse

LABELS=['Outside','Neither eligible','One eligible, below','Both below','A only candidate','B only candidate','Both candidate','Disagreement']
COLORS=['#eee','#bbb','#ddd','#fff','#e6ab02','#377eb8','#2b9b62','#d64c72']

def agreement(qa,qb,valid,threshold=.6):
    a=np.isfinite(qa)&valid;b=np.isfinite(qb)&valid;ca=a&(qa>threshold);cb=b&(qb>threshold)
    s=np.where(valid,1,0).astype('i1');s[(a^b)&~(ca|cb)]=2;s[a&b&~ca&~cb]=3
    s[ca&~b]=4;s[cb&~a]=5;s[ca&cb]=6;s[a&b&(ca^cb)]=7
    return s,a.astype('i1')+2*b.astype('i1')

def read_pass(path):
    with netCDF4.Dataset(path/'orbit85_separation.nc') as d:
        names=['Latitude','Longitude','Pixel_Epoch','L3C_Radiance','Existing_MLCloud','Status','Dominant_Window_Start','Q_Score','Q_AllFrames_Score']
        out={k:arr(d[k][:]) for k in names};out['settings']=json.loads(d.q_settings_json)
        out['all_settings']=json.loads(d.q_all_frames_settings_json)
    out['results']=json.loads((path/'window_summary.json').read_text());out['config']=json.loads((path/'run_config.json').read_text())
    return out

def paired_errors(a,b,pa,pb):
    """Exact same genuine center observations, held out from BOTH models."""
    good=np.isfinite(a['Q_Score'])&np.isfinite(b['Q_Score'])
    yy,xx=np.where(good);times=a['Pixel_Epoch'][good];rows=f.rows_read()
    epoch=np.array([r['epoch'] for r in rows]);frame=np.searchsorted(epoch,times)
    np.testing.assert_array_equal(epoch[frame],times)
    observed=np.full(len(xx),np.nan)
    with netCDF4.Dataset(f.p.L1) as d:
        for i in np.unique(frame):
            take=frame==i;r=rows[i];image=arr(d['Radiance'][int(i)])[::-1]
            observed[take]=image[yy[take]-r['y0'],xx[take]-r['x0']]
    predictions=[]
    for data,path in [(a,pa),(b,pb)]:
        prediction=np.full(len(xx),np.nan);ids=data['Dominant_Window_Start'][good];old=f.OUT;f.OUT=path
        try:
            for r in data['results']:
                take=ids==r['start']
                if not take.any():continue
                assert np.isin(frame[take],r['heldout_frames']).all()
                c=f.load_model(r);dt=times[take]-r['reference_epoch'];x=xx[take]-r['bounds'][0];y=yy[take]-r['bounds'][2]
                prediction[take]=f.sample(c['slow'],c['slow_velocity'],dt,x,y)+f.sample(c['fast'],c['fast_velocity'],dt,x,y)
        finally:f.OUT=old
        predictions.append(prediction)
    assert np.isfinite(observed).all() and all(np.isfinite(p).all() for p in predictions)
    distinct=a['Dominant_Window_Start'][good]!=b['Dominant_Window_Start'][good]
    report={}
    for name,take in [('all_common',np.ones(len(xx),bool)),('distinct_windows_only',distinct)]:
        report[name]=dict(samples=int(take.sum()),rmse_A_kR=float(np.sqrt(np.mean((observed[take]-predictions[0][take])**2))) if take.any() else None,rmse_B_kR=float(np.sqrt(np.mean((observed[take]-predictions[1][take])**2))) if take.any() else None)
    return report

def main(pa,pb,out):
    out.mkdir(exist_ok=True);a=read_pass(pa);b=read_pass(pb)
    for key in ['Latitude','Longitude','Pixel_Epoch','L3C_Radiance']:np.testing.assert_array_equal(a[key],b[key])
    for key in ['source','master','l3c','iterations','fill_small_gaps','motion_smoothing']:assert a['config'].get(key)==b['config'].get(key),key
    keys=['threshold','gaussian_sigma_native_cells','min_local_coverage','min_weighted_residual_energy_kR2','min_training_support']
    for key in keys:
        assert a['settings'][key]==b['settings'][key]==a['all_settings'][key]==b['all_settings'][key],key
    threshold=a['settings']['threshold'];valid=np.isfinite(a['L3C_Radiance'])&np.isfinite(a['Pixel_Epoch'])
    common=set(a['config']['starts'])&set(b['config']['starts'])
    same=(a['Dominant_Window_Start']==b['Dominant_Window_Start'])&np.isin(a['Dominant_Window_Start'],list(common))&valid
    maps={};summary=dict(pass_A=str(pa),pass_B=str(pb),threshold=threshold,common_window_starts=sorted(common),same_window_pixels=int(same.sum()),paired_withheld_prediction_errors=paired_errors(a,b,pa,pb),modes={})
    sa=(a['Status']>=2)&valid;sb=(b['Status']>=2)&valid
    summary['supported_component_coverage']=dict(valid_l3c_pixels=int(valid.sum()),A_fraction=float(sa.sum()/valid.sum()),B_fraction=float(sb.sum()/valid.sum()),either_fraction=float((sa|sb).sum()/valid.sum()),newly_supported_from_B=int((sb&~sa).sum()))
    for mode,name in [('withheld','Q_Score'),('all_frames','Q_AllFrames_Score')]:
        qa,qb=a[name],b[name];s,available=agreement(qa,qb,valid,threshold);maps[mode]=(s,available)
        ea=np.isfinite(qa)&valid;eb=np.isfinite(qb)&valid
        summary['modes'][mode]=dict(eligible_A=int(ea.sum()),eligible_B=int(eb.sum()),eligible_either=int((ea|eb).sum()),newly_eligible_from_B=int((eb&~ea).sum()),categories={LABELS[i]:int((s==i).sum()) for i in range(8)},both_candidate_distinct_windows=int(((s==6)&~same).sum()),both_candidate_shared_window=int(((s==6)&same).sum()))
    with netCDF4.Dataset(out/'orbit85_shift_comparison.nc','w') as d:
        d.createDimension('cross_track',valid.shape[0]);d.createDimension('along_track',valid.shape[1]);dims=('cross_track','along_track')
        d.description='Window-placement sensitivity experiment. No maximum-score selection, no probability calibration. All-frame results include training observations; withheld contributors individually excluded their exact source timestamp. Adjacent overlapping windows are correlated. q depends on each pass brightness-split gauge.'
        for name,data in [(k,a[k]) for k in ['Latitude','Longitude','Pixel_Epoch','L3C_Radiance','Existing_MLCloud']]+[('Q_A_Heldout',a['Q_Score']),('Q_B_Heldout',b['Q_Score']),('Q_A_AllFrames',a['Q_AllFrames_Score']),('Q_B_AllFrames',b['Q_AllFrames_Score'])]:
            v=d.createVariable(name,'f8' if name=='Pixel_Epoch' else 'f4',dims,zlib=True,fill_value=np.nan);v[:]=data
        for mode,(s,available) in maps.items():
            v=d.createVariable(mode+'_Agreement','i1',dims,zlib=True);v[:]=s;v.flag_values=np.arange(8,dtype='i1');v.flag_meanings='outside neither_eligible one_eligible_below both_below A_only_candidate B_only_candidate both_candidate disagreement'
            v=d.createVariable(mode+'_Eligible_Passes','i1',dims,zlib=True);v[:]=available;v.flag_values=np.arange(4,dtype='i1');v.flag_meanings='neither A B both'
        v=d.createVariable('Same_Dominant_Window','i1',dims,zlib=True);v[:]=same.astype('i1')
        for name,data in [('Dominant_Window_A',a['Dominant_Window_Start']),('Dominant_Window_B',b['Dominant_Window_Start'])]:d.createVariable(name,'i2',dims,zlib=True)[:]=data
    f.savejson(out/'summary.json',summary)
    lim=np.nanpercentile(a['L3C_Radiance'],[2,98]);cm=ListedColormap(COLORS);norm=BoundaryNorm(np.arange(-.5,8.5),8)
    for mode,name in [('withheld','Q_Score'),('all_frames','Q_AllFrames_Score')]:
        s,_=maps[mode]
        for suffix,x0,x1 in [('whole_orbit',0,valid.shape[1])]+[(f'segment_{i+1}',l,r) for i,(l,r) in enumerate([(0,1950),(1950,3900),(3900,5850),(5850,7801)])]:
            fig,axes=plt.subplots(6,1,figsize=(20,14),sharex=True,layout='constrained')
            panels=[(a['L3C_Radiance'],'Existing L3C radiance','gray',lim),(a['Existing_MLCloud'],'Existing MLCloud input','viridis',[0,100]),(a[name],'Pass A: original starts | q','viridis',[0,1]),(b[name],'Pass B: starts shifted by 12 frames | q','viridis',[0,1])]
            for ax,(data,title,cmap,limits) in zip(axes,panels):
                im=ax.imshow(data[:,x0:x1],origin='lower',aspect='auto',extent=[x0,x1,0,valid.shape[0]],cmap=cmap,vmin=limits[0],vmax=limits[1],interpolation='nearest');ax.set(title=title);ax.set_facecolor('#bbb');fig.colorbar(im,ax=ax,pad=.01)
            im=axes[4].imshow(s[:,x0:x1],origin='lower',aspect='auto',extent=[x0,x1,0,valid.shape[0]],cmap=cm,norm=norm,interpolation='nearest');cb=fig.colorbar(im,ax=axes[4],ticks=range(8),pad=.01);cb.ax.set_yticklabels(LABELS,fontsize=8);axes[4].set_title('Agreement categories | single-pass candidates are separate from disagreements')
            stable=(s==6)&~same;data=np.where(stable,a['L3C_Radiance'],np.nan)
            im=axes[5].imshow(data[:,x0:x1],origin='lower',aspect='auto',extent=[x0,x1,0,valid.shape[0]],cmap='gray',vmin=lim[0],vmax=lim[1],interpolation='nearest');axes[5].set(title=f'L3C radiance where both distinct windows have q > {threshold} (not a cloud probability)',xlabel='Master along-track cell');axes[5].set_facecolor('#bbb');fig.colorbar(im,ax=axes[5],pad=.01,label='kR')
            fig.suptitle('Orbit 00004 | 85 km | '+('Withheld from each contributing window; correlated validation data' if mode=='withheld' else 'All-frame window-placement diagnostic; includes training observations'))
            fig.savefig(out/f'{mode}_{suffix}.png',dpi=120);plt.close(fig)
    (out/'README.md').write_text('''# Shifted-window comparison

Pass A uses starts 0,25,50,...; B uses 12,37,62,... plus the original endpoint windows. Endpoints shared by both passes do not provide distinct-window confirmation. Same_Dominant_Window flags them, and the final masked-radiance panel excludes them.

Categories distinguish both-above-threshold, only-one-eligible candidates, disagreement when both are eligible, below-threshold evidence, and no evidence. No maximum-q map or undifferentiated union is presented. Ineligible is not cloud-free. A candidate is not a calibrated cloud detection.

Withheld maps only use scores from models that excluded that source timestamp. All-frame maps include training observations and show descriptive fit sensitivity. Overlapping windows share data, so agreement is not independent replication. q is computed in each pass's aligned brightness gauge and is sensitive to the slow/fast brightness split. The paired RMSE comparison uses exactly the same genuine observations held out from both models; its total prediction is gauge-invariant. This comparison is not an untouched test because withheld errors also select models.

Radiance masked in comparison plots is the existing coadded L3C, used as geographic context; q itself uses single-frame L1R at Pixel_Epoch. The two q fields and eligibility provenance remain separate in orbit85_shift_comparison.nc.

Reproduce:
```
python orbit_full.py fit --fill-small-gaps --start-offset 12 --workers 6 --iterations 32
python orbit_full.py assemble --fill-small-gaps --start-offset 12
python combine_orbit_passes.py
```
''',encoding='utf-8')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--pass-a',type=Path,default=f.p.ROOT/'orbit85_full_filled');ap.add_argument('--pass-b',type=Path,default=f.p.ROOT/'orbit85_full_filled_shift12');ap.add_argument('--output',type=Path,default=f.p.ROOT/'orbit85_shift_comparison');args=ap.parse_args();main(args.pass_a,args.pass_b,args.output)
