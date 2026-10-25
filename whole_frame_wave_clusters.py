"""Find four strongest juwavelet clusters in a whole fitted slow frame."""
import slow_wave_phase_speed as s
from slow_wave_phase_speed import np, plt, Path, json, f, arr, netCDF4
from scipy.ndimage import distance_transform_edt
import argparse
from wave_cluster_limited import identify_cluster2d

def main(frame=955):
    out=f.p.ROOT/f'whole_frame_clusters_{frame:04d}';out.mkdir(exist_ok=True)
    transform=s.wavelet_module(Path('C:/Users/domin/juwavelet'))
    rows=f.rows_read();row=rows[frame];results=json.loads((f.OUT/'window_summary.json').read_text())
    r=min([r for r in results if r['start']<=frame<=r['end'] and r['selected_model'] in ('single_slow','two_motion')],key=lambda r:abs(r['reference_frame']-frame))
    c=f.load_model(r);yy,xx=np.indices((300,300));x=xx+row['x0']-r['bounds'][0];y=yy+row['y0']-r['bounds'][2];dt=row['epoch']-r['reference_epoch']
    with netCDF4.Dataset(f.p.L1) as d: obs=arr(d['Radiance'][frame])[::-1]
    obs[(obs<0)|(obs>65535)]=np.nan
    slow=f.sample(c['slow'],c['slow_velocity'],dt,x,y)
    support=f.sample(c['slow_support'].astype(float),c['slow_velocity'],dt,x,y)
    valid=np.isfinite(obs)&np.isfinite(slow)&(support>=8)
    # Remove a plane; taper smoothly over 20 cells (40 nominal km) at all
    # missing-data/footprint edges. Do not extrapolate outside the footprint.
    A=np.column_stack([np.ones(valid.sum()),xx[valid],yy[valid]])
    b=np.linalg.lstsq(A,slow[valid],rcond=None)[0]
    residual=np.where(valid,slow-(b[0]+b[1]*xx+b[2]*yy),0.)
    distance=distance_transform_edt(np.pad(valid,1))[1:-1,1:-1]
    taper=.5-.5*np.cos(np.pi*np.clip(distance/20,0,1))
    scale=np.std(residual[distance>=20]);proc=residual*taper/scale
    minimum=12.;dj=.125;js=int(np.floor(np.log2(200/minimum)/dj))+1
    print('Transforming full frame',flush=True)
    bank=transform.decompose2d(proc.T,dx=2,dy=2,s0=minimum,dj=dj,js=js,jt=36,aspect=1,mode='scaled',dtype=np.complex64)
    print('Grouping top four clusters',flush=True)
    amps,idxs,labels=identify_cluster2d(bank,min_amp=.25,thr=.1,max_clusters=4)
    coeff=bank['decomposition'];maps=[];info=[]
    for k,amp in enumerate(amps):
        weights=np.zeros(coeff.shape[:2]);partial=np.zeros_like(coeff)
        for si in range(js):
            for ti in range(36):
                mask=labels[si,ti]==k
                partial[si,ti][mask]=coeff[si,ti][mask]
                weights[si,ti]=np.abs(coeff[si,ti][mask]).sum()
        total=weights.sum();wavelength=float(np.sum(weights*bank['period'][:,None])/total)
        angle=float(np.degrees(np.angle(np.sum(weights*np.exp(2j*bank['theta'][None,:])))/2)%180)
        rec=transform.reconstruct2d(dict(bank,decomposition=partial)).real.T*scale
        maps.append(np.where(valid,rec,np.nan));si,ti,px,py=idxs[k]
        info.append(dict(rank=k+1,peak_amplitude=float(amp),weighted_wavelength_km=wavelength,wave_axis_degrees=angle,peak_wavelength_km=float(bank['period'][si]),peak_axis_degrees=float(np.degrees(bank['theta'][ti])),peak_pixel_xy=[int(px),int(py)],peak_distance_from_edge_km=float(distance[py,px]*2),wavelength_at_search_edge=bool(si in (0,js-1))))
        del partial
    assert len(maps)==4,'Fewer than four clusters exceed threshold; report fewer rather than invent clusters.'
    np.savez_compressed(out/'clusters.npz',observed=obs,slow=np.where(valid,slow,np.nan),detrended=residual,taper=taper,processed=proc,clusters=np.array(maps),valid=valid,period_km=bank['period'],axis_degrees=np.degrees(bank['theta']),spectrum=np.mean(np.abs(coeff)**2,axis=(2,3)))
    summary=dict(frame=frame,window=[r['start'],r['end']],source='fitted slow spatial template evaluated at frame time',ranking='descending peak wavelet coefficient magnitude, matching juwavelet identify_cluster2d',clusters=info,settings=dict(min_wavelength_km=minimum,max_wavelength_km=float(bank['period'][-1]),angle_step_degrees=5,min_amp=.25,cluster_threshold=.1,taper_width_km=40,spacing_km=2),notes='Wavelengths use nominal Mastergrid spacing. Angles are wave-normal axes relative to native +x, modulo 180 degrees; a single frame does not determine travel direction. Boundary taper reduces but does not eliminate edge effects. Clusters are not necessarily independent physical waves.')
    (out/'summary.json').write_text(json.dumps(summary,indent=2))
    fig,axes=plt.subplots(2,3,figsize=(16,10),layout='constrained')
    images=[np.where(valid,slow,np.nan),np.where(valid,residual,np.nan)]+maps
    titles=['Whole fitted slow frame','Plane removed (before taper)']+[f'Cluster {v["rank"]}: {v["weighted_wavelength_km"]:.1f} km | axis {v["wave_axis_degrees"]:.1f}°' for v in info]
    for i,(ax,a,title) in enumerate(zip(axes.flat,images,titles)):
        lim=np.nanpercentile(abs(a),99) if i else None
        im=ax.imshow(a,origin='lower',extent=[0,600,0,600],cmap='gray' if i==0 else 'RdBu_r',vmin=-lim if lim else None,vmax=lim)
        ax.plot(301,301,'k+',ms=7);ax.set(title=title,xlabel='Native x (nominal km)',ylabel='Native y (nominal km)');fig.colorbar(im,ax=ax,label='kRayleigh')
        if i>=2:
            px,py=info[i-2]['peak_pixel_xy'];ax.plot(2*px+1,2*py+1,'o',mfc='none',mec='black',ms=8)
    fig.suptitle(f'Frame {frame} — four strongest juwavelet clusters across the whole slow frame\n+ frame-center reference; ○ cluster peak; each panel has its own color scale')
    fig.savefig(out/'top_four_clusters.png',dpi=150);plt.close(fig)
    print(json.dumps(summary,indent=2),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--frame',type=int,default=955);args=ap.parse_args()
    f.OUT=f.p.ROOT/'orbit85_full_filled';main(args.frame)
