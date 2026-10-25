"""Self-contained browser animation and GIF for missing-pixel investigation."""
from diagnose_orbit_windows import *
from PIL import Image
import base64,io
OUT=p.ROOT/'orbit85_diagnostics'/'animation_1650'

def png(data):
    b=io.BytesIO();Image.fromarray(data).save(b,format='PNG');return 'data:image/png;base64,'+base64.b64encode(b.getvalue()).decode()

def main():
    OUT.mkdir(exist_ok=True);start=1650;rows=full.rows_read();cube,times,meta=cube_for(start,rows);x0,x1,y0,y1=meta['bounds'];geom=np.zeros(cube.shape,bool)
    with netCDF4.Dataset(p.L1) as d:
        for j,i in enumerate(range(start,start+50)):
            r=rows[i];a=arr(d['Latitude'][i])[::-1];b=arr(d['Longitude'][i])[::-1];g=np.isfinite(a)&np.isfinite(b)
            ax,bx=max(x0,r['x0']),min(x1,r['x0']+300);ay,by=max(y0,r['y0']),min(y1,r['y0']+300)
            if ax<bx and ay<by:geom[j,ay-y0:by-y0,ax-x0:bx-x0]=g[ay-r['y0']:by-r['y0'],ax-r['x0']:bx-r['x0']]
    raw=[];smooth=[]
    for a in cube:
        raw.append(p.highpass(a));good=np.isfinite(a);w=gaussian_filter(good.astype(float),1);b=np.where(w>.98,gaussian_filter(np.nan_to_num(a),1)/np.maximum(w,1e-9),np.nan);smooth.append(p.highpass(b))
    raw=np.array(raw);smooth=np.array(smooth);lo,hi=np.nanpercentile(cube,[2,98]);lim=float(np.nanpercentile(np.abs(raw),98));colors=np.array([[75,85,99],[204,68,170],[231,163,62],[59,156,117]],dtype='uint8')
    mask=np.zeros(cube.shape,'uint8');mask[geom&~np.isfinite(cube)]=1;mask[np.isfinite(cube)&~np.isfinite(raw)]=2;mask[np.isfinite(raw)]=3
    data=[];gifframes=[]
    def rgba(a,cmap,vmin,vmax):
        z=(plt.get_cmap(cmap)(np.clip((np.nan_to_num(a,nan=vmin)-vmin)/(vmax-vmin),0,1))[:,:,:3]*255).astype('uint8');z[~np.isfinite(a)]=[75,85,99];return z
    for j in range(50):
        ims=[rgba(cube[j],'gray',lo,hi),colors[mask[j]],rgba(raw[j],'RdBu_r',-lim,lim),rgba(smooth[j],'RdBu_r',-lim,lim)]
        bits=[base64.b64encode(np.packbits(np.isfinite(a[j])[::-1].ravel(),bitorder='big').tobytes()).decode() for a in (cube,raw,smooth)]
        data.append(dict(frame=start+j,seconds=float(times[j]-times[0]),images=[png(im[::-1]) for im in ims],bits=bits,counts=np.bincount(mask[j].ravel(),minlength=4).tolist()))
        fig,axes=plt.subplots(2,2,figsize=(10,8),layout='constrained')
        titles=['Observed radiance','Why pixels are excluded','Original matching texture','After 1-pixel smoothing (matching only)']
        for ax,im,title in zip(axes.ravel(),ims,titles):ax.imshow(im,origin='lower');ax.set_title(title,fontsize=10);ax.set_xticks([]);ax.set_yticks([])
        fig.suptitle(f'85 km | L1R frame {start+j} | +{times[j]-times[0]:.1f} s | fixed Mastergrid',fontsize=13)
        fig.supxlabel('Mask: gray = outside footprint · magenta = missing radiance · orange = filter edge · green = usable texture\nFixed scales throughout. Playback is accelerated; loop reset is not physical motion.',fontsize=8)
        buf=io.BytesIO();fig.savefig(buf,format='png',dpi=100);plt.close(fig);buf.seek(0);gifframes.append(Image.open(buf).convert('RGB').copy())
    gifframes[0].save(OUT/'missing_pixels_animation.gif',save_all=True,append_images=gifframes[1:],duration=[160]*49+[1000],loop=0,optimize=False)
    html='''<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>AWE — Investigate missing pixels and weak texture</title><style>
body{font:16px system-ui;margin:24px auto;max-width:1050px;padding:0 16px;background:#101821;color:#e8eef6}h1{font-size:26px}p{line-height:1.5;color:#bdcbdc}.controls{position:sticky;top:0;background:#182332;padding:14px;z-index:2;border-radius:8px;display:flex;gap:12px;align-items:center;flex-wrap:wrap}button,select{padding:8px;border:0;border-radius:5px}input{flex:1;min-width:180px}.grid{display:grid;grid-template-columns:1fr 1fr;gap:18px}canvas{width:100%;image-rendering:auto;cursor:crosshair}h2{font-size:16px;margin-bottom:8px}.legend{display:flex;flex-wrap:wrap;gap:14px;margin:15px 0}.legend span:before{content:'■';color:var(--c);margin-right:6px}#patch{background:#182332;padding:12px;border-radius:8px;line-height:1.6}@media(max-width:650px){.grid{grid-template-columns:1fr}}
</style></head><body><h1>Why did this patch fail?</h1><p>Frames 1650–1699, orbit 4, v02 at 85 km. All panels stay on the same Mastergrid cells. Click any panel to place the same <b>25 × 25 reference patch</b> in all views. Watch missing coverage move separately from the scene texture.</p>
<div class="controls"><button id="play">Pause</button><button id="prev">◀</button><input id="slider" type="range" min="0" max="49" value="0"><button id="next">▶</button><select id="speed"><option value="160">Fast</option><option value="400">Slow</option><option value="1100">Approx. real time</option></select><b id="label"></b></div>
<div class="legend"><span style="--c:#4b5563">Outside current valid coordinate footprint</span><span style="--c:#cc44aa">Missing radiance inside footprint</span><span style="--c:#e7a33e">Excluded by baseline filter support</span><span style="--c:#3b9c75">Usable baseline texture pixel</span></div>
<div class="grid"><div><h2>Original radiance <small id="radianceScale"></small></h2><canvas id="c0"></canvas></div><div><h2>Coverage and exclusion reason</h2><canvas id="c1"></canvas></div><div><h2>Original matching texture <small id="textureScale"></small></h2><canvas id="c2"></canvas></div><div><h2>Matching texture after 1-pixel smoothing</h2><canvas id="c3"></canvas></div></div><p id="patch"></p>
<p>Texture panels use the <b>same fixed blue–white–red scale</b>, with zero at white. Radiance uses one fixed scale for the full sequence. Smoothing is used only for matching; the fitted radiance data are not smoothed. The mask panel describes baseline preprocessing; patch counts also show support after smoothing.</p>
<p>A reference patch with even one missing texture pixel is rejected by the current matcher. All-valid does not guarantee a match: the later-frame search must also contain a valid patch with sufficient correlation. The yellow box is fixed in grid coordinates, not feature tracking. The jump at the loop boundary is a playback reset.</p>
<script>const D=__DATA__,W=340,H=300,X0=__X0__,Y0=__Y0__;let frame=0,playing=true,px=170,py=150,timer;const canv=Array.from({length:4},(_,i)=>document.getElementById('c'+i));canv.forEach(c=>{c.width=W;c.height=H});const images=D.map(d=>d.images.map(src=>{let i=new Image();i.src=src;return i}));const bits=D.map(d=>d.bits.map(s=>Uint8Array.from(atob(s),c=>c.charCodeAt(0))));function valid(b,x,y){let k=y*W+x;return (b[k>>3]>>(7-(k&7)))&1}function render(){let d=D[frame];document.getElementById('label').textContent=`Frame ${d.frame} · +${d.seconds.toFixed(1)} s`;document.getElementById('slider').value=frame;canv.forEach((c,k)=>{let ctx=c.getContext('2d');ctx.drawImage(images[frame][k],0,0);ctx.strokeStyle='#ffe36b';ctx.lineWidth=1;ctx.strokeRect(px-12.5,py-12.5,25,25)});let n=bits[frame].map(b=>{let count=0;for(let y=py-12;y<=py+12;y++)for(let x=px-12;x<=px+12;x++)count+=valid(b,x,y);return count});document.getElementById('patch').innerHTML=`<b>Patch center: Master cell (${X0+px}, ${Y0+H-1-py})</b><br>Valid radiance: ${n[0]}/625 · Original texture: ${n[1]}/625 · Smoothed texture: ${n[2]}/625<br>Original reference-patch support: <b>${n[1]===625?'PASS — all pixels valid':'REJECT — contains missing texture pixels'}</b>. Smoothed support: <b>${n[2]===625?'PASS':'REJECT'}</b>.`;}
function tick(){clearTimeout(timer);timer=setTimeout(()=>{if(playing){frame=(frame+1)%D.length;render()}tick()},Number(document.getElementById('speed').value))}function pause(){playing=false;document.getElementById('play').textContent='Play'}document.getElementById('play').onclick=()=>{playing=!playing;document.getElementById('play').textContent=playing?'Pause':'Play'};document.getElementById('slider').oninput=e=>{pause();frame=Number(e.target.value);render()};document.getElementById('prev').onclick=()=>{pause();frame=(frame+49)%50;render()};document.getElementById('next').onclick=()=>{pause();frame=(frame+1)%50;render()};canv.forEach(c=>c.onclick=e=>{let r=c.getBoundingClientRect();px=Math.max(12,Math.min(W-13,Math.floor((e.clientX-r.left)/r.width*W)));py=Math.max(12,Math.min(H-13,Math.floor((e.clientY-r.top)/r.height*H)));render()});document.getElementById('radianceScale').textContent='(__LO__ to __HI__ kR)';document.getElementById('textureScale').textContent='(±__LIM__ kR)';Promise.all(images.flat().map(i=>i.decode())).then(()=>{render();tick()});</script></body></html>'''
    html=html.replace('__DATA__',json.dumps(data)).replace('__X0__',str(x0)).replace('__Y0__',str(y0)).replace('__LO__',f'{lo:.2f}').replace('__HI__',f'{hi:.2f}').replace('__LIM__',f'{lim:.2f}')
    (OUT/'investigate_missing_pixels.html').write_text(html,encoding='utf-8')
    assert cube.shape==(50,300,340)
    assert Image.open(OUT/'missing_pixels_animation.gif').n_frames==50
    for j,d in enumerate(data):
        for k,a in enumerate((cube,raw,smooth)):
            decoded=np.unpackbits(np.frombuffer(base64.b64decode(d['bits'][k]),dtype='uint8'),bitorder='big')[:300*340].reshape(300,340)
            np.testing.assert_array_equal(decoded,np.isfinite(a[j])[::-1])
    savejson(OUT/'metadata.json',dict(source=str(p.L1),frames=[start,start+49],bounds=meta['bounds'],radiance_limits_kR=[float(lo),float(hi)],texture_limits_kR=[-lim,lim],gif_frames=50,verification='50 GIF frames; all packed interactive masks exactly match arrays; fixed scales across time',note='Mask separates coordinate footprint, missing radiance, and baseline filter-support exclusion. No physical cloud classification.'))
    print(str(OUT),flush=True)

if __name__=='__main__':main()
