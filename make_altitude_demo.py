"""Compare supplied 85 km / 15 km geolocation at the original demo's epochs."""
import orbit_pilot as p
from orbit_pilot import np, netCDF4, arr, plt, json
from scipy.interpolate import LinearNDInterpolator
from scipy.ndimage import gaussian_filter
from small_gap_fill import fill_small
from PIL import Image
import base64, zipfile

OUT=p.ROOT/'parallax_demo_85km_15km'
SOURCE15=p.L1.parent/'remap_85km_to_15km'/(p.L1.stem+'_remap15km.nc')

def main():
    OUT.mkdir(exist_ok=True)
    old=np.load(p.ROOT/'parallax_demo/projected_frames.npz')
    with netCDF4.Dataset(p.L1) as a, netCDF4.Dataset(SOURCE15) as b:
        epochs=arr(a['Epoch'][:]);np.testing.assert_array_equal(epochs,arr(b['Epoch'][:]))
        assert str(a.OH_layer_altitude).startswith('85.') and str(b.OH_layer_altitude).startswith('15.')
        matches=[];missing=[]
        for frame,e in zip(old['frame_index'],old['epoch']):
            found=np.flatnonzero(epochs==e)
            if len(found)==1:matches.append((int(frame),int(found[0])))
            elif not len(found):missing.append(int(frame))
            else:raise ValueError('Duplicate source epoch')
        labels,indices=np.asarray(matches).T;times=epochs[indices]
        coords={h:(arr(d['Latitude'][indices]),arr(d['Longitude'][indices])) for h,d in [(85,a),(15,b)]}
        original=arr(a['Radiance'][indices]);np.testing.assert_array_equal(original,arr(b['Radiance'][indices]))
        method=b.remap_method
    original[(original<0)|(original>65535)]=np.nan
    filled=original.copy();masks=np.zeros(original.shape,bool)
    for j in range(len(indices)):
        geom=np.ones(original.shape[1:],bool)
        for lat,lon in coords.values():geom&=np.isfinite(lat[j])&np.isfinite(lon[j])
        filled[j],masks[j],_=fill_small(original[j],geom)
    np.testing.assert_array_equal(filled[np.isfinite(original)],original[np.isfinite(original)])
    x,y=old['longitude'],old['latitude'];xx,yy=np.meshgrid(x,y);step=float(x[1]-x[0])
    cubes={};textures={}
    for h,(lat,lon) in coords.items():
        cube=[]
        for j,i in enumerate(indices):
            good=np.isfinite(lat[j])&np.isfinite(lon[j])
            # Retain missing-radiance vertices so large gaps are not bridged.
            interp=LinearNDInterpolator(np.column_stack((lon[j][good],lat[j][good])),filled[j][good])
            cube.append(interp(xx,yy).astype('f4'))
            print(f'{h} km: projected frame {i}',flush=True)
        cubes[h]=np.asarray(cube);res=[]
        for image in cube:
            ok=np.isfinite(image);weight=gaussian_filter(ok.astype(float),10)
            smooth=gaussian_filter(np.nan_to_num(image),10)/np.maximum(weight,1e-10)
            res.append(np.where(ok&(weight>.95),image-smooth,np.nan))
        textures[h]=np.asarray(res)
    common={h:np.isfinite(cubes[h]).all(axis=0) for h in cubes}
    limits=np.nanpercentile(np.concatenate([a.ravel() for a in cubes.values()]),[2,98])
    limit=float(np.nanpercentile(np.concatenate([np.abs(a).ravel() for a in textures.values()]),98))
    np.savez_compressed(OUT/'projected_frames.npz',radiance_85km=cubes[85],radiance_15km=cubes[15],spatial_residual_85km=textures[85],spatial_residual_15km=textures[15],latitude=y,longitude=x,frame_index=indices,matching_original_demo_index=labels,epoch=times,common_coverage_85km=common[85],common_coverage_15km=common[15])
    np.savez_compressed(OUT/'source_interpolation_provenance.npz',original_radiance=original,filled_radiance=filled,interpolated=masks,frame_index=indices,epoch=times)
    plt.rcParams.update({'font.size':10,'figure.facecolor':'#f5f7fa','axes.facecolor':'#dedfe2'})
    fig,axes=plt.subplots(2,2,figsize=(12,10),layout='constrained');ims=[]
    extent=[x[0]-step/2,x[-1]+step/2,y[0]-step/2,y[-1]+step/2]
    for row in range(2):
        for col,h in enumerate([85,15]):
            ax=axes[row,col];data=cubes[h] if row==0 else textures[h]
            lo,hi=limits if row==0 else [-limit,limit]
            im=ax.imshow(data[0],origin='lower',extent=extent,cmap='gray' if row==0 else 'RdBu_r',vmin=lo,vmax=hi,interpolation='nearest',aspect=1/np.cos(np.deg2rad(y.mean())))
            ims.append((im,data));ax.contour(xx,yy,common[h].astype(float),levels=[.5],colors=['#d59a19'],linewidths=.7)
            ax.set(title=f'{h} km mapping | '+('Radiance' if row==0 else 'Radiance minus spatial smooth'),xlabel='Longitude (degrees E)',ylabel='Latitude (degrees N)');ax.grid(alpha=.15)
            fig.colorbar(im,ax=ax,shrink=.8,label='kRayleigh')
    title=fig.suptitle('');paths=[];ref=int(np.flatnonzero(labels==955)[0])
    for j,i in enumerate(indices):
        for im,data in ims:im.set_data(data[j])
        title.set_text(f'AWE altitude comparison | v02 frame {i} | time from reference {times[j]-times[ref]:+.1f} s')
        path=OUT/f'frame_{i:04d}.png';fig.savefig(path,dpi=100);paths.append(path)
    plt.close(fig)
    images=[Image.open(path).convert('RGB') for path in paths]
    durations=[int(round(dt*120/1.1/10)*10) for dt in np.diff(times)]+[1000]
    images[0].save(OUT/'altitude_comparison.gif',save_all=True,append_images=images[1:],duration=durations,loop=0)
    images[ref].save(OUT/'preview.png')
    embedded=['data:image/png;base64,'+base64.b64encode(path.read_bytes()).decode() for path in paths]
    html='''<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>AWE: 85 km versus 15 km mapping</title><style>body{font:16px system-ui;background:#eef1f6;color:#182536;max-width:1250px;margin:24px auto;padding:0 20px}h1{font-size:26px}img{width:100%}.controls{position:sticky;top:0;background:#e3e9f2;padding:10px;display:flex;gap:10px;flex-wrap:wrap;align-items:center}button,select{font:inherit;padding:7px}input{flex:1}p{line-height:1.5}</style></head><body><h1>AWE: 85 km versus 15 km mapping</h1><p>The same observations, shown at two assumed mapping altitudes. Left: 85 km. Right: 15 km. Top: radiance. Bottom: small-scale spatial contrast. Both columns use identical map coordinates and color scales.</p><div class="controls"><button id="play">Play</button><button id="back">←</button><input id="slider" type="range" min="0" max="__MAX__" value="__REF__"><button id="next">→</button><select id="speed"><option value="120">Normal</option><option value="300">Slow</option><option value="1100">Approx. real time</option></select><button id="blink">Blink endpoints</button><span id="label"></span></div><img id="frame"><p>__COUNT__ exact timestamp matches to the original demo's frames 930–979. Original indices missing from v02: __MISSING__. Actual time gaps are retained in playback. Gold outlines coverage shared by every selected frame at that altitude. Gray indicates missing data. The view uses the original demo's fixed 0.02-degree grid and bounds.</p><p>Only enclosed missing-radiance holes ≤16 pixels and ≤5×5 pixels are filled before mapping; both altitudes use identical filled values. Larger gaps remain missing. Spatial contrast subtracts a normalized Gaussian smooth (sigma 0.2 degrees); it is not a separated layer or a motion correction.</p><p>The 15 km coordinates come from the supplied line-of-sight remapping. Comparing apparent motion can help investigate altitude-dependent parallax, but does not by itself establish a feature's altitude. The loop reset is not physical motion.</p><p><a href="altitude_comparison.gif">Animated GIF</a> · <a href="metadata.json">Metadata</a></p><script>const images=__IMAGES__,ids=__IDS__,times=__TIMES__;let n=__REF__,playing=false,blinking=false,timer=null;const slider=document.getElementById('slider'),frame=document.getElementById('frame'),label=document.getElementById('label');function show(v){n=(v+images.length)%images.length;slider.value=n;frame.src=images[n];label.textContent=`v02 ${ids[n]} · ${(times[n]-times[__REF__]).toFixed(1)} s`}function stop(){clearTimeout(timer);playing=false;blinking=false;document.getElementById('play').textContent='Play'}function tick(){let dt=n<images.length-1?(times[n+1]-times[n])/1.1:5;timer=setTimeout(()=>{show(blinking?(n===0?images.length-1:0):n+1);if(playing)tick()},blinking?650:dt*Number(document.getElementById('speed').value))}document.getElementById('play').onclick=()=>{if(playing)stop();else{playing=true;document.getElementById('play').textContent='Pause';tick()}};slider.oninput=()=>{stop();show(+slider.value)};document.getElementById('next').onclick=()=>{stop();show(n+1)};document.getElementById('back').onclick=()=>{stop();show(n-1)};document.getElementById('speed').onchange=stop;document.getElementById('blink').onclick=()=>{stop();show(0);playing=true;blinking=true;document.getElementById('play').textContent='Pause';tick()};show(n);</script></body></html>'''
    for key,value in {'MAX':len(indices)-1,'REF':ref,'COUNT':len(indices),'MISSING':str(missing),'IMAGES':json.dumps(embedded),'IDS':json.dumps(indices.tolist()),'TIMES':json.dumps(times.tolist())}.items():html=html.replace('__'+key+'__',str(value))
    (OUT/'index.html').write_text(html,encoding='utf-8')
    assert '__IMAGES__' not in html and Image.open(OUT/'altitude_comparison.gif').n_frames==len(indices)
    meta=dict(source_85km=str(p.L1),source_15km=str(SOURCE15),mapping_method_15km=method,frames=len(indices),v02_indices=indices.tolist(),missing_original_indices=missing,reference_v02_frame=int(indices[ref]),duration_seconds=float(times[-1]-times[0]),grid_spacing_degrees=step,radiance_limits_kR=limits.tolist(),spatial_contrast_limit_kR=limit,filled_source_samples=int(masks.sum()),checks='PASS: exact epochs; identical input radiance at both altitudes; original finite radiance unchanged; animation frame count; fixed shared map grid and scales',interpolation='Linear interpolation retains NaN vertices; conservative small-hole filling before projection; no motion registration')
    p.savejson(OUT/'metadata.json',meta)
    (OUT/'README.md').write_text('''# AWE 85 km versus 15 km demo

Open index.html in a browser. Images are embedded, so the interactive animation works from this single HTML file. Other download links use adjacent files in this ZIP.

The original demo epochs are matched exactly to v02. Both mappings use identical radiance and timestamps. See metadata.json for the skipped original index and source details. Small enclosed holes are filled before projection, with masks retained in the workspace file source_interpolation_provenance.npz. Projected data and spatial contrasts are retained in workspace file projected_frames.npz. Numerical arrays are omitted from this portable demo ZIP. Larger missing regions remain missing.

Columns compare mapping altitude, not separated physical layers. All panels use the same fixed geographic grid; scales are shared within each row. The common-coverage outline is computed independently for each altitude. Numerical source arrays retain the original disk orientation. Reproduce using make_altitude_demo.py in the parent workspace.
''',encoding='utf-8')
    with zipfile.ZipFile(p.ROOT/'parallax_demo_85km_15km.zip','w',zipfile.ZIP_DEFLATED) as z:
        for name in ['index.html','altitude_comparison.gif','preview.png','metadata.json','README.md']:
            path=OUT/name;z.write(path,arcname=OUT.name+'/'+path.name)

    print(json.dumps(meta,indent=2),flush=True)

if __name__=='__main__':main()
