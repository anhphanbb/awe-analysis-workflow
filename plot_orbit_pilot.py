"""Plot the pilot against L3C on verified native Master cells and pixel times."""
from orbit_pilot import *

def sample(a, v, dt, x=None, y=None, support=False):
    if x is None:y,x=np.indices(a.shape,dtype=float)
    coords=[y-v[1]*dt,x-v[0]*dt]
    valid=map_coordinates(np.isfinite(a).astype(float),coords,order=1,mode='constant',cval=0)
    z=map_coordinates(np.nan_to_num(a),coords,order=1,mode='constant',cval=0)
    return np.where(valid>.999,z/np.maximum(valid,1e-9),np.nan)

def main():
    results=json.loads((OUT/'pilot_summary.json').read_text());inv=json.loads((OUT/'inventory.json').read_text())
    z=np.load(OUT/'l3c_reference.npz');rad=z['radiance'];cloud=z['mlcloud'];epochs=z['epoch']
    fastmap=np.full(rad.shape,np.nan,dtype='f4');slowmap=fastmap.copy();choice=np.full(rad.shape,-1,dtype='i2');best=np.full(rad.shape,np.inf,dtype='f4')
    comparisons=[];fits=[r for r in results if r['status']=='fit'];limits=[]
    for r in results:
        if r['status']=='fit':continue
        x0,x1,y0,y1=r['bounds'];fig,ax=plt.subplots(1,3,figsize=(15,5),layout='constrained')
        for a,data,title,cm,lim,unit in zip(ax[:2],[rad[y0:y1,x0:x1],cloud[y0:y1,x0:x1]],['Existing L3C radiance','Existing L3C MLCloud'],['gray','viridis'],[np.nanpercentile(rad[y0:y1,x0:x1],[2,98]),[0,100]],['kR','%']):
            im=a.imshow(data,origin='lower',extent=[x0,x1,y0,y1],cmap=cm,vmin=lim[0],vmax=lim[1]);a.set(title=title,xlabel='Master along-track cell',ylabel='Master cross-track cell');fig.colorbar(im,ax=a,shrink=.7,label=unit)
        ax[2].axis('off');ax[2].text(.5,.5,'INSUFFICIENT MOTION EVIDENCE\n\nNo supported motion mode\npassed the current checks.\n\nCloud presence is undetermined.',ha='center',va='center',fontsize=13)
        fig.suptitle(f"85 km pilot | frames {r['start']}–{r['end']}");fig.savefig(OUT/f"window_{r['start']:04d}/comparison.png",dpi=150);plt.close(fig)
    for r in fits:
        folder=OUT/f"window_{r['start']:04d}";c=np.load(folder/f"checkpoint_{r.get('checkpoint_iteration',64):03d}.npz");maps=np.load(folder/'reference_maps.npz')
        x0,x1,y0,y1=r['bounds'];e=epochs[y0:y1,x0:x1];dt=e-r['reference_epoch'];safe=np.where(np.isfinite(dt),dt,0)
        fast=sample(c['fast'],c['fast_velocity'],safe);slow=sample(c['slow'],c['slow_velocity'],safe)
        fs=sample(c['fast_support'].astype(float),c['fast_velocity'],safe);ss=sample(c['slow_support'].astype(float),c['slow_velocity'],safe)
        # Evaluate at actual L3C source-frame times, only in each window's central 25 frames.
        ts=maps['relative_seconds'];keep=np.isfinite(e)&np.isfinite(rad[y0:y1,x0:x1])&(dt>=ts[12])&(dt<=ts[37])&(fs>=8)&(ss>=8)
        fast=np.where(keep,fast,np.nan);slow=np.where(keep,slow,np.nan)
        if r.get('model_components',2)==1:fast[:]=np.nan
        r['l3c_matched_pixels']=int(np.isfinite(fast).sum())
        eligible=keep&np.isfinite(fast)&np.isfinite(slow)&(np.abs(dt)<best[y0:y1,x0:x1])
        fastmap[y0:y1,x0:x1][eligible]=fast[eligible];slowmap[y0:y1,x0:x1][eligible]=slow[eligible];choice[y0:y1,x0:x1][eligible]=r['start'];best[y0:y1,x0:x1][eligible]=np.abs(dt[eligible])
        np.savez_compressed(folder/'l3c_time_matched.npz',fast=fast,slow=slow,epoch=e,supported=keep)
        lim=np.nanpercentile(maps['observed'],[2,98]);contrast=max(.05,float(np.nanpercentile(np.abs(maps['fast']),98)));limits.append(contrast)
        fig,ax=plt.subplots(2,3,figsize=(15,9),layout='constrained')
        panels=[(rad[y0:y1,x0:x1],'Existing L3C radiance','gray',lim,'kR'),(cloud[y0:y1,x0:x1],'Existing L3C MLCloud','viridis',[0,100],'%'),(fast,'Fast contrast at L3C pixel times','RdBu_r',[-contrast,contrast],'kR'),(maps['observed'],f"L1R reference frame {r['reference_frame']}",'gray',lim,'kR'),(maps['slow'],'Slow template at reference time','gray',lim,'kR'),(maps['residual'],'Reference-frame reconstruction residual','RdBu_r',[-.8,.8],'kR')]
        for a,(data,title,cm,lims,unit) in zip(ax.ravel(),panels):
            im=a.imshow(data,origin='lower',extent=[x0,x1,y0,y1],aspect='equal',cmap=cm,vmin=lims[0],vmax=lims[1]);a.set(title=title,xlabel='Master along-track cell',ylabel='Master cross-track cell');a.set_facecolor('#ddd');fig.colorbar(im,ax=a,shrink=.75,label=unit)
        if r.get('model_components',2)==1:ax[0,2].text(.5,.5,'Single slow-motion model\nNo fast component estimated',transform=ax[0,2].transAxes,ha='center',va='center')
        fig.suptitle(f"85 km pilot | frames {r['start']}–{r['end']} | {r['quality']}\nFast contrast has an arbitrary zero; gray = outside support/time coverage",fontsize=14)
        fig.savefig(folder/'comparison.png',dpi=150);plt.close(fig)
    # Compare adjacent overlapping models at the same physical time and grid cells.
    for a,b in zip(results[::2],results[1::2]):
        if a['status']!='fit' or b['status']!='fit':continue
        ca=np.load(OUT/f"window_{a['start']:04d}/checkpoint_{a.get('checkpoint_iteration',64):03d}.npz");cb=np.load(OUT/f"window_{b['start']:04d}/checkpoint_{b.get('checkpoint_iteration',64):03d}.npz")
        x0=max(a['bounds'][0],b['bounds'][0]);x1=min(a['bounds'][1],b['bounds'][1]);y0=max(a['bounds'][2],b['bounds'][2]);y1=min(a['bounds'][3],b['bounds'][3]);yy,xx=np.mgrid[y0:y1,x0:x1]
        epoch=(a['reference_epoch']+b['reference_epoch'])/2
        rec=dict(windows=[a['start'],b['start']],comparison_epoch=epoch)
        for layer in ('slow','fast'):
            if layer=='fast' and min(a.get('model_components',2),b.get('model_components',2))==1:
                rec[layer]={'pixels':0};continue
            aa=sample(ca[layer],ca[layer+'_velocity'],epoch-a['reference_epoch'],xx-a['bounds'][0],yy-a['bounds'][2]);bb=sample(cb[layer],cb[layer+'_velocity'],epoch-b['reference_epoch'],xx-b['bounds'][0],yy-b['bounds'][2])
            sa=sample(ca[layer+'_support'].astype(float),ca[layer+'_velocity'],epoch-a['reference_epoch'],xx-a['bounds'][0],yy-a['bounds'][2]);sb=sample(cb[layer+'_support'].astype(float),cb[layer+'_velocity'],epoch-b['reference_epoch'],xx-b['bounds'][0],yy-b['bounds'][2])
            ok=np.isfinite(aa)&np.isfinite(bb)&(sa>=8)&(sb>=8)
            if ok.sum()<100:rec[layer]={'pixels':int(ok.sum())};continue
            diff=aa[ok]-bb[ok];offset=float(np.median(diff));corr=float(np.corrcoef(aa[ok],bb[ok])[0,1]);rec[layer]=dict(pixels=int(ok.sum()),median_offset_kR=offset,offset_removed_rmse_kR=float(np.sqrt(np.mean((diff-offset)**2))),correlation=corr)
        comparisons.append(rec)
    savejson(OUT/'overlap_consistency.json',comparisons)
    savejson(OUT/'pilot_summary.json',results)
    la,lo=master()
    with netCDF4.Dataset(OUT/'pilot_on_mastergrid.nc','w') as d:
        d.createDimension('cross_track',rad.shape[0]);d.createDimension('along_track',rad.shape[1])
        d.description='Exploratory six-window pilot only. Fast component is signed contrast, not cloud probability or absolute cloud radiance. Values evaluated at existing L3C Pixel_Epoch; L3C coadding is not reproduced. Window gauges are independent; no global offset stitching.'
        d.reference_altitude_km=85.;d.master_file=str(MASTER)
        for name,data,unit in [('Latitude',la,'degrees_north'),('Longitude',lo,'degrees_east'),('L3C_Radiance',rad,'kRayleigh'),('Existing_MLCloud',cloud,'percent'),('Pilot_Fast_Contrast',fastmap,'kRayleigh'),('Pilot_Slow',slowmap,'kRayleigh'),('Pixel_Epoch',epochs,'source Epoch seconds')]:
            v=d.createVariable(name,'f8' if name=='Pixel_Epoch' else 'f4',('cross_track','along_track'),zlib=True,fill_value=np.nan);v[:]=data;v.units=unit
        v=d.createVariable('Pilot_Window_Start','i2',('cross_track','along_track'),zlib=True,fill_value=-1);v[:]=choice
    fig,ax=plt.subplots(3,1,figsize=(18,8),sharex=True,layout='constrained');clim=float(np.nanpercentile(np.abs(fastmap),98)) if np.isfinite(fastmap).any() else 1.
    for a,data,title,cm,lim,unit in zip(ax,[rad,cloud,fastmap],['Existing L3C radiance','Existing MLCloud map','Motion-derived fast contrast — pilot coverage only'],['gray','viridis','RdBu_r'],[np.nanpercentile(rad,[2,98]),[0,100],[-clim,clim]],['kR','%','kR']):
        im=a.imshow(data,origin='lower',aspect='auto',cmap=cm,vmin=lim[0],vmax=lim[1]);a.set(title=title,ylabel='Cross-track cell');a.set_facecolor('#ddd');fig.colorbar(im,ax=a,pad=.01,label=unit)
    for r in results:
        if r['status']=='fit':continue
        x0,x1,y0,y1=r['bounds'];ax[-1].plot([x0,x1,x1,x0,x0],[y0,y0,y1,y1,y0],color='#6f478c',ls='--',lw=1)
    ax[-1].text(6600,440,'Insufficient motion evidence',color='#6f478c',ha='center',fontsize=9)
    ax[-1].set_xlabel('Master along-track cell (nominal 2 km/cell; vertical axis stretched)');fig.suptitle('Orbit 00004 | 85 km | common Mastergrid | fast component evaluated at L3C pixel times',fontsize=14);fig.savefig(OUT/'l3c_cloud_and_pilot.png',dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(10,6),layout='constrained')
    for r in fits:
        h=r['history'];line=ax.plot([h0['iteration'] for h0 in h],[h0['heldout_rmse_kR'] for h0 in h],marker='o',label=f"{r['start']}–{r['end']}")[0]
        ax.axhline(r['slow_only_heldout_rmse_kR'],color=line.get_color(),ls=':',alpha=.5)
    ax.set(xscale='log',xlabel='Iterations',ylabel='Held-out RMSE (kR)',title='Held-out reconstruction | dotted lines: slow-only baselines');ax.legend();ax.grid(alpha=.2);fig.savefig(OUT/'heldout_convergence.png',dpi=160);plt.close(fig)
    lines=['# Orbit 4: 85 km pilot','',f"Inventory: {inv['frames']} frames over {inv['duration_seconds']/60:.2f} minutes; {len(inv['gaps'])} retained timestamp gaps. Master shape: 551 × 7801. All valid L1R coordinates and L3C coordinates matched the Master within {inv['master_max_l1r_coordinate_error_deg']:.2g} degrees (float32 rounding).",'', '## Six-window experiment','', '| Frames | Classification | Slow-only holdout RMSE | Two-layer RMSE (64) | Improvement | L3C-matched pixels |','|---|---|---:|---:|---:|---:|']
    for r in results:
        if r['status']=='fit':lines.append(f"| {r['start']}–{r['end']} | {r['quality']} | {r['slow_only_heldout_rmse_kR']:.3f} | {r['history'][-1]['heldout_rmse_kR']:.3f} | {r['heldout_improvement_fraction']:.1%} | {r['l3c_matched_pixels']} |")
        else:lines.append(f"| {r['start']}–{r['end']} | unresolved: {r['reason']} | — | — | — | — |")
    lines+=['','## Method and limits','','- Native L1R radiance is placed directly on verified Master cells; no new spatial interpolation. Source files are untouched.','- Each window has 50 observations; only even local indices estimate velocities and train templates. Odd indices evaluate predictions. Neighbouring frames remain temporally correlated; this is an exploratory holdout, not independent atmospheric validation.','- Two strongest distinct supported motion modes are selected without fixed slow/fast speed bins. The slower is labeled slow; neither label proves a physical identity. Full candidates and diagnostic histograms are saved.','- Constant translation is expressed in native along/cross-track grid cells, nominally 2 km. This is a local approximation over each window, not an Earth-fixed velocity retrieval.','- Minimum template support is 8 training frames. Unsupported updates retain previous values; displayed cloud estimates require final support of at least 8 for both layers.','- Held-out sample masks are fixed from iteration 1 onward. Baseline is the initial slow-only temporal-median template at the same slow velocity.','- The classification pilot_supported requires >5% held-out improvement and at least 3 frame pairs per mode; it is a heuristic flag, not calibrated confidence.','- L3C is the existing coadded product. Fast contrast is evaluated at its per-pixel representative Epoch, restricted to the central 25 frames of each window. It is not a reconstruction of the upstream 19-frame coadd.','- MLCloud is the existing percentage-valued product. It is not interchangeable with the signed motion-derived fast radiance.','- Each fast template has median zero; absolute radiance partition remains ambiguous. The orbit pilot map chooses the closest reference time without blending or solving between-window offsets.','- Only the pilot regions are processed; gray areas in its map are unprocessed or unsupported, not cloud-free.','', '## Overlap checks','']
    for c in comparisons:
        f=c['fast'];lines.append(f"- Windows {c['windows']}: fast-component correlation {f.get('correlation',float('nan')):.3f}; offset-removed RMSE {f.get('offset_removed_rmse_kR',float('nan')):.3f} kR; {f['pixels']} supported common pixels.")
    lines+=['','## Files','','- `l3c_cloud_and_pilot.png`: orbit comparison.','- `window_*/comparison.png`: detailed L3C, MLCloud, fast contrast and reconstruction panels.','- `heldout_convergence.png`, `overlap_consistency.json`: diagnostics.','- `pilot_on_mastergrid.nc`: geolocated L3C reference and sparse pilot fields, with window provenance.','- `inventory.json`, `frame_inventory.csv`: complete source inventory.','- `window_*/checkpoint_*.npz`: fitted templates at 16/32/64 iterations.','','Reproduce using Python 3.12 with NumPy, SciPy, matplotlib and netCDF4:','```','python orbit_pilot.py all','python plot_orbit_pilot.py','```','The scripts use the existing workspace `.packages` directory. Completed windows are reused; remove only an explicitly chosen window output to refit it.']
    lines += ['', '## Single-motion scenes and interpretation', '',
        'A second layer is optional. The estimator now accepts one supported motion mode and fits a single translated template in that case. A mode below 0.18 native cells/s with evidence from at least three frame pairs is provisionally labeled single slow motion; this speed threshold is heuristic. The displayed fast field is left missing when no fast component is estimated, rather than presenting zeros as measured cloud absence.',
        'A known single-motion synthetic scene passes the one-mode check. In these real pilot windows, the four successful fits have two supported modes. The two late windows do not meet the current minimum of 15 candidates within a 0.085-cell/s neighborhood for even one mode. They therefore remain insufficient-evidence cases, not evidence for cloud presence or absence.',
        'The next whole-orbit pass should retain three outcomes: one supported slow motion, two supported motions with predictive benefit, and insufficient evidence. Resolve brightness offsets between adjacent templates before constructing a continuous fast-contrast mosaic. The current pilot has not performed that stitching.',
        '', 'Verification: run `python verify_orbit_pilot.py`; see `verification.json` for the completed numerical and product checks.']
    (OUT/'README.md').write_text('\n'.join(lines),encoding='utf-8')
    print(json.dumps(dict(overlap=comparisons,pilot_pixels=int(np.isfinite(fastmap).sum()),results=[{k:r[k] for k in ('start','quality','heldout_improvement_fraction','l3c_matched_pixels') if k in r} for r in results]),indent=2),flush=True)

if __name__=='__main__':main()
