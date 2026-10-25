import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent / '.packages'))
import numpy as np
import netCDF4
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import LinearNDInterpolator
from scipy.ndimage import gaussian_filter
from PIL import Image
import json

ROOT = Path(__file__).parent
OUT = ROOT / 'parallax_demo'
OUT.mkdir(exist_ok=True)
SOURCE = r'Y:\soc\l1r\2023\11\awe_l1r_q20_2023326T0547_00004_v01.nc'
indices = np.arange(930,980)
with netCDF4.Dataset(SOURCE) as ds:
    lat = np.asarray(ds['Latitude'][930:980])
    lon = np.asarray(ds['Longitude'][930:980])
    rad = np.asarray(ds['Radiance'][930:980])
    epoch = np.asarray(ds['Epoch'][930:980])
    altitude = ds.OH_layer_altitude
step = .02
x = np.arange(np.ceil(np.nanmin(lon[25])/step)*step, np.nanmax(lon[25]), step)
y = np.arange(np.ceil(np.nanmin(lat[25])/step)*step, np.nanmax(lat[25]), step)
xx, yy = np.meshgrid(x,y)
cube=[]
for j,i in enumerate(indices):
    good=np.isfinite(lat[j]) & np.isfinite(lon[j]) & np.isfinite(rad[j]) & (rad[j]>=0) & (rad[j]<=65535)
    interp=LinearNDInterpolator(np.column_stack((lon[j][good],lat[j][good])),rad[j][good])
    cube.append(interp(xx,yy).astype('float32'))
    print(f'Projected {i}',flush=True)
cube=np.array(cube)
common=np.all(np.isfinite(cube),axis=0)
residual=[]
for a in cube:
    ok=np.isfinite(a)
    weight=gaussian_filter(ok.astype(float),10)
    smooth=gaussian_filter(np.nan_to_num(a),10)/np.maximum(weight,1e-10)
    residual.append(np.where(ok & (weight>.95),a-smooth,np.nan))
residual=np.array(residual)
vmin,vmax=np.nanpercentile(cube,[2,98])
limit=np.nanpercentile(np.abs(residual),98)
np.savez_compressed(OUT/'projected_frames.npz',radiance=cube,spatial_residual=residual,latitude=y,longitude=x,frame_index=indices,epoch=epoch,common_coverage=common)
plt.rcParams.update({'font.size':10,'figure.facecolor':'#f5f7fa','axes.facecolor':'#dedfe2'})
fig,axes=plt.subplots(1,2,figsize=(12,6.3),layout='constrained')
extent=[x[0]-step/2,x[-1]+step/2,y[0]-step/2,y[-1]+step/2]
ims=[]
for ax,a,cmap,lo,hi,title in zip(axes,[cube[0],residual[0]],['gray','RdBu_r'],[vmin,-limit],[vmax,limit],['Calibrated radiance','Small-scale structure (radiance − spatial smooth)']):
    im=ax.imshow(a,origin='lower',extent=extent,cmap=cmap,vmin=lo,vmax=hi,interpolation='nearest',aspect=1/np.cos(np.deg2rad(y.mean())))
    ims.append(im)
    ax.contour(xx,yy,common.astype(float),levels=[.5],colors=['#e4a72b'],linewidths=.8)
    ax.set(xlabel='Longitude (°E)',ylabel='Latitude (°N)',title=title)
    ax.grid(alpha=.2)
    fig.colorbar(im,ax=ax,shrink=.72,label='kRayleigh')
title=fig.suptitle('')
fig.text(.5,.005,'Fixed 0.02° grid • file geolocation at 87 km • gold outline: coverage shared by all 50 frames',ha='center',fontsize=10)
paths=[]
for j,i in enumerate(indices):
    for im,a in zip(ims,[cube[j],residual[j]]): im.set_data(a)
    title.set_text(f'AWE Q20 | frame {i} (zero-based) | Δt = {epoch[j]-epoch[25]:+.1f} s from frame 955')
    p=OUT/f'frame_{i:04d}.png'
    fig.savefig(p,dpi=110)
    paths.append(p)
plt.close(fig)
images=[Image.open(p).convert('RGB') for p in paths]
images[0].save(OUT/'fixed_latlon.gif',save_all=True,append_images=images[1:],duration=120,loop=0)
images[25].save(OUT/'preview.png')
meta={'source':SOURCE,'frames_zero_based':[930,979],'reference_frame':955,'duration_seconds':float(epoch[-1]-epoch[0]),'grid_spacing_degrees':step,'reference_altitude':altitude,'radiance_limits_kR':[float(vmin),float(vmax)],'residual_limit_kR':float(limit),'common_pixels':int(common.sum()),'grid_shape':list(common.shape),'spatial_smoothing_sigma_grid_pixels':10}
(OUT/'metadata.json').write_text(json.dumps(meta,indent=2))
html='''<!doctype html><html><head><meta charset="utf-8"><title>AWE fixed latitude/longitude demo</title><style>body{font:16px system-ui;background:#eef1f6;color:#182536;max-width:1340px;margin:24px auto;padding:0 20px}h1{font-size:25px}img{width:100%;border-radius:10px}button,select{font:inherit;padding:7px 14px;margin:5px}input{width:55%}.note{line-height:1.5;color:#43536a}</style></head><body><h1>AWE · 50 frames on a fixed latitude/longitude grid</h1><p>Frames 930–979 (zero-based), centered near 955. About 54 seconds of observations.</p><img id="frame" src="frame_0955.png"><div><button id="play">Play</button><button id="back">←</button><button id="next">→</button><input id="slider" type="range" min="0" max="49" value="25"><span id="label">955</span><select id="speed"><option value="240">Slow</option><option value="120" selected>Normal</option><option value="60">Fast</option></select><button id="blink">Blink 930 / 979</button></div><p class="note">Both panels use the same fixed 0.02° map grid and fixed color scales throughout. Gold outlines the area observed in all 50 frames. Gray is outside the valid footprint. The right panel subtracts a Gaussian spatial smooth (σ = 0.2°, roughly 22 km north–south); it does not register or shift the images.</p><p class="note">Coordinates come directly from the file’s 87 km OH-layer geolocation. Motion on this map is apparent motion relative to that reference surface; it alone does not establish parallax or a feature’s altitude. Genuine scene evolution and geolocation errors can also produce motion. Playback wraps from the last frame to the first.</p><p><a href="fixed_latlon.gif">Animated GIF</a> · <a href="projected_frames.npz">Projected numerical data</a> · <a href="metadata.json">Processing metadata</a></p><script>let n=25,timer=null;const slider=document.querySelector('#slider'),frame=document.querySelector('#frame'),label=document.querySelector('#label');function show(v){n=(v+50)%50;slider.value=n;label.textContent=930+n;frame.src='frame_'+String(930+n).padStart(4,'0')+'.png'}function stop(){clearInterval(timer);timer=null;document.querySelector('#play').textContent='Play'}document.querySelector('#play').onclick=()=>{if(timer){stop()}else{timer=setInterval(()=>show(n+1),+document.querySelector('#speed').value);document.querySelector('#play').textContent='Pause'}};slider.oninput=()=>{stop();show(+slider.value)};document.querySelector('#next').onclick=()=>{stop();show(n+1)};document.querySelector('#back').onclick=()=>{stop();show(n-1)};document.querySelector('#speed').onchange=stop;document.querySelector('#blink').onclick=()=>{stop();show(0);timer=setInterval(()=>show(n===0?49:0),600);document.querySelector('#play').textContent='Pause'};for(let i=930;i<980;i++){let im=new Image();im.src='frame_'+String(i).padStart(4,'0')+'.png'}</script></body></html>'''
(OUT/'index.html').write_text(html,encoding='utf-8')
print(json.dumps(meta),flush=True)
