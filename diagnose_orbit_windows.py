"""Instrument rejected-window matching; controlled preprocessing ablations.

Original orbit products are read-only. All experiments use even local frames.
"""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import orbit_full as full
import orbit_pilot as p
from orbit_pilot import np, json, csv, Path, arr, netCDF4, gaussian_filter, maximum_filter, plt, savejson
from concurrent.futures import ProcessPoolExecutor,as_completed
from collections import Counter
import argparse, warnings
OUT=p.ROOT/'orbit85_diagnostics'
VARIANTS={'baseline':dict(half=12,smooth=0,threshold=.08),
          'lower_texture':dict(half=12,smooth=0,threshold=.04),
          'smooth_1px':dict(half=12,smooth=1,threshold=.08),
          'larger_patch':dict(half=20,smooth=0,threshold=.08)}

def cube_for(start,rows):
    ref=rows[start+25];la,_=p.master();ny,nx=la.shape
    x0=max(0,ref['x0']-20);x1=min(nx,ref['x0']+320);y0=max(0,ref['y0']);y1=min(ny,ref['y0']+300)
    cube=np.full((50,y1-y0,x1-x0),np.nan,'f4');invalid_range=0;finite_source=0
    with netCDF4.Dataset(p.L1) as d:
        for j,i in enumerate(range(start,start+50)):
            r=rows[i];a=arr(d['Radiance'][i])[::-1];finite_source+=int(np.isfinite(a).sum());bad=np.isfinite(a)&((a<0)|(a>65535));invalid_range+=int(bad.sum());a[bad]=np.nan
            ax,bx=max(x0,r['x0']),min(x1,r['x0']+300);ay,by=max(y0,r['y0']),min(y1,r['y0']+300)
            if ax<bx and ay<by:cube[j,ay-y0:by-y0,ax-x0:bx-x0]=a[ay-r['y0']:by-r['y0'],ax-r['x0']:bx-r['x0']]
    times=np.array([rows[i]['epoch']-ref['epoch'] for i in range(start,start+50)])
    return cube,times,dict(finite_source_pixels=finite_source,pixels_rejected_by_radiance_range=invalid_range,bounds=[x0,x1,y0,y1])

def match(cube,times,variant):
    cfg=VARIANTS[variant];h={};noise=[]
    for i in range(0,50,2):
        a=cube[i];diff=np.diff(a,axis=1);v=diff[np.isfinite(diff)];noise.append(float(np.median(np.abs(v-np.median(v)))/(.67448975*np.sqrt(2))))
        if cfg['smooth']:
            good=np.isfinite(a);w=gaussian_filter(good.astype(float),cfg['smooth']);a=np.where(w>.98,gaussian_filter(np.nan_to_num(a),cfg['smooth'])/np.maximum(w,1e-9),np.nan)
        h[i]=p.highpass(a)
    counts=Counter();records=[];stds=[];maxcorr=[]
    for start in (0,12,24):
        for lag in (8,16,24):
            end=start+lag
            if end>=50:continue
            dt=times[end]-times[start];radius=int(np.ceil(dt*1.1));half=cfg['half']
            for y in range(half+radius,cube.shape[1]-half-radius,24):
                for x in range(half+radius,cube.shape[2]-half-radius,24):
                    counts['patches_considered']+=1;t=h[start][y-half:y+half+1,x-half:x+half+1]
                    if not np.isfinite(t).all():counts['reference_missing']+=1;continue
                    std=float(t.std());stds.append(std)
                    if std<cfg['threshold']:counts['texture_below_threshold']+=1;continue
                    counts['patches_searched']+=1
                    c=p.ncc(t,h[end][y-half-radius:y+half+radius+1,x-half-radius:x+half+radius+1]);counts['candidate_positions']+=c.size;counts['valid_candidate_positions']+=int(np.isfinite(c).sum())
                    if not np.isfinite(c).any():counts['no_valid_search_position']+=1;continue
                    best=float(np.nanmax(c));maxcorr.append(best)
                    if best<=.35:counts['correlation_below_threshold']+=1;continue
                    c0=np.nan_to_num(c,nan=-2);ys,xs=np.where((c0==maximum_filter(c0,5))&(c0>.35));accepted=0
                    for k in np.argsort(c0[ys,xs])[::-1][:4]:
                        py,px=int(ys[k]),int(xs[k])
                        if min(py,px)<=0 or max(py,px)>=2*radius:counts['peaks_on_search_boundary']+=1;continue
                        dy,dx=float(py-radius),float(px-radius)
                        for axis in (0,1):
                            z=c[py-1:py+2,px] if axis==0 else c[py,px-1:px+2];den=z[0]-2*z[1]+z[2]
                            delta=float(np.clip(.5*(z[0]-z[2])/den,-.5,.5)) if np.isfinite(z).all() and den<0 else 0.
                            if axis==0:dy+=delta
                            else:dx+=delta
                        records.append([dx/dt,dy/dt,float(c[py,px]),start,end,x,y]);accepted+=1
                    counts['patches_with_matches' if accepted else 'boundary_only_patch']+=1
    a=np.array(records).reshape(-1,7);counts['retained_candidates']=len(a)
    weights=(a[:,2]-.35)**2
    hist,edges,_=np.histogram2d(a[:,0],a[:,1],bins=np.linspace(-1.1,1.1,111),weights=weights);hist=gaussian_filter(hist,1)
    peaks=np.argwhere((hist==maximum_filter(hist,11))&(hist>.15*hist.max()));peaks=sorted(peaks,key=lambda q:hist[tuple(q)],reverse=True)
    modes=[];peakstats=[]
    for peak in peaks:
        center=np.array([edges[peak[0]]+.01,edges[peak[1]]+.01]);near=np.linalg.norm(a[:,:2]-center,axis=1)<.085;n=int(near.sum());peakstats.append(dict(center=center.tolist(),support=n))
        if n<15:continue
        v=np.average(a[near,:2],axis=0,weights=weights[near])
        if any(np.linalg.norm(v-m['velocity'])<.12 for m in modes):continue
        modes.append(dict(velocity=v.tolist(),support=n,median_ncc=float(np.median(a[near,2])),pair_count=len(set(map(tuple,a[near,3:5])))))
    modes=sorted(modes[:2],key=lambda m:np.linalg.norm(m['velocity']))
    if len(modes)>1 and np.linalg.norm(np.subtract(modes[0]['velocity'],modes[1]['velocity']))<.15:modes=modes[:1]
    if len(a)<20:reason='fewer_than_20_candidates';modes=[]
    elif not modes:reason='no_cluster_with_15_candidates'
    elif len(modes)==1 and np.linalg.norm(modes[0]['velocity'])>=.18:reason='only_non_slow_mode'
    elif len(modes)==1:reason='single_slow_mode'
    else:reason='two_modes'
    stats=dict(variant=variant,parameters=cfg,counts=dict(counts),modes=modes,reason=reason,peak_stats=peakstats[:8],median_noise_proxy_kR=float(np.median(noise)),texture_std_percentiles_kR=np.percentile(stds,[0,10,50,90,100]).tolist() if stds else [],best_ncc_percentiles=np.percentile(maxcorr,[0,10,50,90,100]).tolist() if maxcorr else [])
    return stats,a

def diagnose_one(start,rejected):
    folder=OUT/f'window_{start:04d}';folder.mkdir(parents=True,exist_ok=True)
    if (folder/'matching.json').exists():return json.loads((folder/'matching.json').read_text())
    cube,times,source=cube_for(start,full.rows_read());v=cube[::2][np.isfinite(cube[::2])]
    source['radiance_percentiles_kR']=np.percentile(v,[0,.1,1,50,99,99.9,100]).tolist();source['finite_fraction_in_analysis_cube']=float(np.isfinite(cube).mean())
    variants={}
    for variant in VARIANTS if rejected else ['baseline']:
        stats,a=match(cube,times,variant);variants[variant]=stats
        np.save(folder/f'{variant}_candidates.npy',a)
    old=Path(full.OUT/f'window_{start:04d}/motion_candidates.csv')
    if old.exists() and old.stat().st_size>50:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore');old_a=np.loadtxt(old,delimiter=',',skiprows=1).reshape(-1,7)
        new=np.load(folder/'baseline_candidates.npy');np.testing.assert_allclose(new,old_a,atol=1e-7,rtol=1e-6)
    result=dict(start=start,rejected=rejected,source=source,variants=variants,baseline_candidate_reproduction='PASS');savejson(folder/'matching.json',result);return result

def scan(workers):
    OUT.mkdir(exist_ok=True);original=json.loads((full.OUT/'window_summary.json').read_text());rejected=[r['start'] for r in original if r['selected_model'] not in ('single_slow','two_motion')]
    results=[]
    with ProcessPoolExecutor(max_workers=workers) as ex:
        jobs=[ex.submit(diagnose_one,r['start'],r['start'] in rejected) for r in original]
        for future in as_completed(jobs):
            results.append(future.result());savejson(OUT/'progress.json',dict(completed=len(results),total=len(original)));print('Diagnosed',len(results),'/',len(original),flush=True)
    results.sort(key=lambda r:r['start']);savejson(OUT/'matching_summary.json',results)
    print({v:dict(Counter(r['variants'][v]['reason'] for r in results if r['rejected'])) for v in VARIANTS},flush=True)

def refit_one(start,variant):
    folder=OUT/f'fits_{variant}';folder.mkdir(exist_ok=True)
    diag=json.loads((OUT/f'window_{start:04d}/matching.json').read_text())['variants'][variant]
    candidates=np.load(OUT/f'window_{start:04d}/{variant}_candidates.npy');weights=(candidates[:,2]-.35)**2;hist,edges,_=np.histogram2d(candidates[:,0],candidates[:,1],bins=np.linspace(-1.1,1.1,111),weights=weights)
    original_motion=p.motions
    def supplied(cube,times,diagnostic_folder=None):return diag['modes'],candidates,gaussian_filter(hist,1),edges
    try:
        p.motions=supplied;p.OUT=folder;r=p.run_window(start,full.rows_read(),32)
    finally:p.motions=original_motion
    r['selected_model']=full.classify(r);r['variant']=variant
    savejson(folder/f'window_{start:04d}/diagnostic_result.json',r);return r

def refit(workers):
    diags=json.loads((OUT/'matching_summary.json').read_text());jobspec=[]
    # Prespecified variants, no holdout-based choice of which settings to fit.
    for r in diags:
        if not r['rejected']:continue
        for v in ('smooth_1px','larger_patch','lower_texture'):
            s=r['variants'][v]
            if s['reason'] in ('single_slow_mode','two_modes') and s['modes']!=r['variants']['baseline']['modes']:jobspec.append((r['start'],v))
    savejson(OUT/'refit_plan.json',jobspec);results=[]
    with ProcessPoolExecutor(max_workers=workers) as ex:
        jobs=[ex.submit(refit_one,*j) for j in jobspec]
        for job in as_completed(jobs):
            r=job.result();results.append(r);savejson(OUT/'refit_progress.json',dict(completed=len(results),total=len(jobspec)));print('Refit',len(results),'/',len(jobspec),r['start'],r['variant'],r['selected_model'],flush=True)
    savejson(OUT/'refit_summary.json',results)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['scan','refit']);ap.add_argument('--workers',type=int,default=6);a=ap.parse_args()
    (scan if a.mode=='scan' else refit)(a.workers)
