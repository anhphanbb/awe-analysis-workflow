"""Build the standalone orbit 19 analysis HTML from verified numerical results."""

from build_orbit19_actual import OUT,ROOT,np,plt,json

from PIL import Image

import base64,io,html



def data_uri(payload):return 'data:image/png;base64,'+base64.b64encode(payload).decode()

def figure_uri(fig):

    b=io.BytesIO();fig.savefig(b,format='png',dpi=145,bbox_inches='tight',facecolor='white');plt.close(fig);return data_uri(b.getvalue())

def map_uri(a,cmap,low,high):

    rgba=plt.get_cmap(cmap)(np.clip((np.nan_to_num(a,nan=low)-low)/(high-low),0,1),bytes=True)

    rgba[~np.isfinite(a)]=[225,230,233,255]

    b=io.BytesIO();Image.fromarray(rgba[::-1]).save(b,format='PNG',optimize=True);return data_uri(b.getvalue())



def main():

    meta=json.loads((OUT/'fit_summary.json').read_text());wave=json.loads((OUT/'grid128/step05_thr020/wavelet_summary.json').read_text())

    d=np.load(OUT/'projected_frames.npz');seq=np.load(OUT/'fitted_sequence.npz');w=np.load(OUT/'grid128/step05_thr020/wavelet_results.npz');wd=np.load(OUT/'grid128/step05_thr020/projected_frames.npz')

    history=json.loads((OUT/'actual_results_128/metrics.json').read_text())

    rad=d['radiance'];slow=seq['slow'];fast=seq['fast'];res=seq['residual']

    np.testing.assert_allclose((rad-slow-fast)[np.isfinite(res)],res[np.isfinite(res)],atol=2e-6)

    low,high=np.nanpercentile(rad,[2,98]);flim=max(.05,float(np.nanpercentile(abs(fast),98)));rlim=max(.05,float(np.nanpercentile(abs(res),98)))

    frames=[]

    for j in range(len(rad)):

        frames.append(dict(frame=int(d['frame_index'][j]),utc=meta['utc'][j],images=[map_uri(rad[j],'gray',low,high),map_uri(slow[j],'gray',low,high),map_uri(fast[j],'RdBu_r',-flim,flim),map_uri(res[j],'RdBu_r',-rlim,rlim)]))

    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})

    fig,ax=plt.subplots(figsize=(9,2.7),layout='constrained')

    ax.plot([r['iteration'] for r in history],[r['reconstruction_rmse_kR'] for r in history],color='#006e73',lw=2)

    ax.set(xlabel='Completed iterations',ylabel='Reconstruction RMSE (kR)');ax.grid(alpha=.15)

    convergence=figure_uri(fig)

    fig,ax=plt.subplots(figsize=(10,4.5),layout='constrained')
    spectrum=w['spectrum'];period=w['period_km']
    im=ax.imshow(spectrum,origin='lower',aspect='auto',extent=[-2.5,177.5,-.5,len(period)-.5],cmap='magma')
    tick=np.arange(0,len(period),4);ax.set(yticks=tick,yticklabels=[f'{period[k]:.0f}' for k in tick],xlabel='Wave-normal axis from east (degrees)',ylabel='Wavelength (km)')
    fig.colorbar(im,ax=ax,label='Mean normalized coefficient power',shrink=.85)
    spectral=figure_uri(fig)

    extent=[d['longitude'][0]-.01,d['longitude'][-1]+.01,d['latitude'][0]-.01,d['latitude'][-1]+.01]

    wave_extent=[wd['longitude'][0]-.01,wd['longitude'][-1]+.01,wd['latitude'][0]-.01,wd['latitude'][-1]+.01]
    clusterlim=max(.01,float(np.nanpercentile(abs(w['clusters']),99)))

    clusters=[]

    for i,c in enumerate(wave['clusters']):

        fig,ax=plt.subplots(figsize=(5,4.6),layout='constrained')

        im=ax.imshow(w['clusters'][i],origin='lower',extent=wave_extent,aspect=1/np.cos(np.deg2rad(d['latitude'].mean())),cmap='RdBu_r',vmin=-clusterlim,vmax=clusterlim)

        ax.set(xlabel='Longitude (°E)',ylabel='Latitude (°N)');ax.set_facecolor('#e1e6e9')

        ax.plot(c['peak_longitude'],c['peak_latitude'],'o',mfc='none',mec='#202b35',ms=8)

        fig.colorbar(im,ax=ax,label='Reconstructed contrast (kR)',shrink=.8)

        c=dict(c,image=figure_uri(fig));clusters.append(c)

    fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')

    lim=float(np.nanpercentile(abs(w['detrended']),98))

    for ax,a,title,cmap,vmin,vmax in [(axes[0],w['detrended'],'Slow field after plane removal','RdBu_r',-lim,lim),(axes[1],w['taper'],'40 km boundary taper','viridis',0,1)]:

        im=ax.imshow(a,origin='lower',extent=wave_extent,cmap=cmap,vmin=vmin,vmax=vmax,aspect=1/np.cos(np.deg2rad(d['latitude'].mean())));ax.set(title=title,xlabel='Longitude (°E)',ylabel='Latitude (°N)');fig.colorbar(im,ax=ax,shrink=.8)

    preparation=figure_uri(fig)

    cluster_cards=''.join(f'''<article class="card cluster"><div class="cluster-head"><span class="rank">{c['rank']:02}</span><div><h3>{c['weighted_wavelength_km']:.1f} km <span>· {c['wave_normal_axis_degrees']:.1f}° axis</span></h3><p>Amplitude-weighted cluster estimates</p></div></div><img src="{c['image']}" alt="Wavelet cluster {c['rank']}: reconstructed signed contrast"><p class="caption">Peak coefficient at {c['peak_wavelength_km']:.1f} km / {c['peak_axis_degrees']:.0f}°. Peak is {c['peak_edge_distance_km']:.0f} km from the nearest missing-data edge. {'Inside the taper region: edge-sensitive.' if c['peak_edge_distance_km']<40 else ''} {'Peak reaches a search limit.' if c['search_edge'] else ''}</p></article>''' for c in clusters)

    latrange=f"{d['latitude'][0]:.2f} to {d['latitude'][-1]:.2f}°N";lonrange=f"{d['longitude'][0]:.2f} to {d['longitude'][-1]:.2f}°E"

    # Crop the concentric arcs near the top of the same numerical slow field.

    from matplotlib.patches import Rectangle

    crop_bounds=(-59.4,-56.4,-4.8,-1.8)

    xsel=(d['longitude']>=crop_bounds[0])&(d['longitude']<=crop_bounds[1])

    ysel=(d['latitude']>=crop_bounds[2])&(d['latitude']<=crop_bounds[3])

    crop=slow[25][np.ix_(ysel,xsel)]

    np.savez_compressed(OUT/'concentric_waves_crop.npz',slow=crop,latitude=d['latitude'][ysel],longitude=d['longitude'][xsel],frame=1467)

    fig,axes=plt.subplots(1,2,figsize=(12,6),layout='constrained')

    for ax,a,ext,title in [(axes[0],slow[25],extent,'Full fitted slow component'),(axes[1],crop,crop_bounds,'Concentric-wave region · same radiance scale')]:

        im=ax.imshow(a,origin='lower',extent=ext,cmap='gray',vmin=low,vmax=high,aspect=1/np.cos(np.deg2rad(d['latitude'].mean())))

        ax.set(title=title,xlabel='Longitude (°E)',ylabel='Latitude (°N)');ax.set_facecolor('#e1e6e9')

    axes[0].add_patch(Rectangle((crop_bounds[0],crop_bounds[2]),crop_bounds[1]-crop_bounds[0],crop_bounds[3]-crop_bounds[2],fill=False,edgecolor='#d68a26',lw=2))

    fig.colorbar(im,ax=axes,label='Slow component radiance (kR)',shrink=.8)

    fig.savefig(OUT/'concentric_waves_crop.png',dpi=170,bbox_inches='tight')

    crop_image=figure_uri(fig)

    payload=dict(frames=frames,reference=25)

    replacements={'__CROP__':crop_image,'__DATA__':json.dumps(payload,separators=(',',':')),'__CLUSTERS__':cluster_cards,'__CONVERGENCE__':convergence,'__SPECTRAL__':spectral,'__PREPARATION__':preparation,

        '__FIRST__':meta['first_utc'][11:23],'__LAST__':meta['last_utc'][11:23],'__REF__':meta['reference_utc'][11:23],

        '__DURATION__':f"{meta['duration_seconds']:.1f}",'__INITIAL__':f"{meta['initial_rmse_kR']:.3f}",'__FINAL__':f"{meta['final_rmse_kR']:.3f}",

        '__PIXELS__':f"{meta['evaluation_pixels_per_frame']:,}",'__LAT__':latrange,'__LON__':lonrange,

        '__LOW__':f'{low:.2f}','__HIGH__':f'{high:.2f}','__FAST__':f'{flim:.2f}','__RESIDUAL__':f'{rlim:.2f}',

        '__SLOWV__':', '.join(f'{v:+.4f}' for v in meta['slow_velocity_pixels_per_s']),

        '__FASTV__':', '.join(f'{v:+.4f}' for v in meta['fast_velocity_pixels_per_s']),

        '__SUPPORT__':f"{meta['slow_peak_support']} slow / {meta['fast_peak_support']} fast candidates",'__NCLUSTERS__':str(len(clusters)),

        '__METADATA__':html.escape(json.dumps(dict(fit=meta,wavelet=wave),indent=2)),

        '__ASPECT__':str(rad.shape[2]*np.cos(np.deg2rad(d['latitude'].mean()))/rad.shape[1])}

    page=(ROOT/'orbit19_actual_template.html').read_text(encoding='utf-8')

    for key,value in replacements.items():page=page.replace(key,value)

    assert '__DATA__' not in page and len(frames)==50

    target=OUT/'orbit19_052340.html';target.write_text(page,encoding='utf-8')

    (OUT/'html_verification.json').write_text(json.dumps(dict(status='PASS',frames=len(frames),reference_frame=1467,clusters=len(clusters),reconstruction_identity='PASS',analysis_scope='concentric-wave crop only',wavelet_input_shape=list(w['slow'].shape),standalone=True,bytes=target.stat().st_size),indent=2))

    print('HTML COMPLETE',target,target.stat().st_size,flush=True)



if __name__=='__main__':main()

