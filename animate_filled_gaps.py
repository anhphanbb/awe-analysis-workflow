"""Fill small enclosed radiance gaps for a separate diagnostic animation."""
from diagnose_orbit_windows import *
from small_gap_fill import fill_small
from PIL import Image
import base64,io
OUT=p.ROOT/'orbit85_diagnostics'/'animation_1650_filled'

def texture(a):
    good=np.isfinite(a);w=gaussian_filter(good.astype(float),1);b=np.where(w>.98,gaussian_filter(np.nan_to_num(a),1)/np.maximum(w,1e-9),np.nan)
    return p.highpass(b)

def main():
    OUT.mkdir(exist_ok=True);rows=full.rows_read();cube,times,meta=cube_for(1650,rows);x0,x1,y0,y1=meta['bounds'];geom=np.zeros(cube.shape,bool)
    with netCDF4.Dataset(p.L1) as d:
        for j,i in enumerate(range(1650,1700)):
            r=rows[i];g=np.isfinite(arr(d['Latitude'][i])[::-1])&np.isfinite(arr(d['Longitude'][i])[::-1]);ax,bx=max(x0,r['x0']),min(x1,r['x0']+300);ay,by=max(y0,r['y0']),min(y1,r['y0']+300)
            if ax<bx and ay<by:geom[j,ay-y0:by-y0,ax-x0:bx-x0]=g[ay-r['y0']:by-r['y0'],ax-r['x0']:bx-r['x0']]
    new=[];masks=[];counts=[];sizes=[]
    for a,g in zip(cube,geom):
        b,m,sz=fill_small(a,g);new.append(b);masks.append(m);counts.append(int(m.sum()));sizes.extend(sz)
    new=np.array(new);masks=np.array(masks);before=np.array([texture(a) for a in cube]);after=np.array([texture(a) for a in new])
    np.testing.assert_array_equal(new[np.isfinite(cube)],cube[np.isfinite(cube)]);assert not np.any(masks&np.isfinite(cube));assert geom[masks].all();assert np.array_equal(np.isfinite(new),np.isfinite(cube)|masks)
    np.savez_compressed(OUT/'filled_analysis_frames.npz',radiance=new,interpolated=masks,frame_index=np.arange(1650,1700),relative_seconds=times-times[0],bounds=meta['bounds'])
    lo,hi=np.nanpercentile(cube,[2,98]);lim=float(np.nanpercentile(np.abs(before),98));gifs=[];images=[]
    for j in range(50):
        fig,ax=plt.subplots(2,2,figsize=(10,8),layout='constrained')
        for a,im,title,cm,limits in zip(ax.ravel(),[cube[j],new[j],before[j],after[j]],['Original radiance',f'Tiny gaps filled ({counts[j]} estimated pixels)','Matching texture: original missing pixels','Matching texture: after small-gap filling'],['gray','gray','RdBu_r','RdBu_r'],[[lo,hi],[lo,hi],[-lim,lim],[-lim,lim]]):
            cmap=plt.get_cmap(cm).copy();cmap.set_bad('#4b5563');a.imshow(im,origin='lower',cmap=cmap,vmin=limits[0],vmax=limits[1]);a.set_title(title,fontsize=10);a.set_xticks([]);a.set_yticks([])
        fig.suptitle(f'85 km | frame {1650+j} | +{times[j]-times[0]:.1f} s | fixed Mastergrid',fontsize=13)
        fig.supxlabel(f'Fill only enclosed gaps ≤16 cells and ≤5 × 5 cells, inside valid geometry. Larger gaps stay missing.\nBoth textures use identical 1-pixel smoothing. Fixed scales: radiance {lo:.2f}–{hi:.2f} kR; texture ±{lim:.2f} kR.',fontsize=8)
        b=io.BytesIO();fig.savefig(b,format='png',dpi=100);plt.close(fig);images.append('data:image/png;base64,'+base64.b64encode(b.getvalue()).decode());b.seek(0);gifs.append(Image.open(b).convert('RGB').copy())
    gifs[0].save(OUT/'small_gaps_filled.gif',save_all=True,append_images=gifs[1:],duration=[160]*49+[1000],loop=0,optimize=False);gifs[25].save(OUT/'preview.png')
    html='''<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>AWE — Tiny gaps filled: before and after</title><style>body{max-width:1050px;margin:24px auto;padding:0 16px;background:#101821;color:#eef3fa;font:16px system-ui}p{line-height:1.5}img{width:100%;height:auto}.controls{position:sticky;top:0;background:#263548;padding:12px;display:flex;gap:10px;align-items:center;flex-wrap:wrap}input{flex:1}button,select{padding:8px}h1{font-size:25px}</style></head><body><h1>Tiny missing pixels: before and after filling</h1><p>Frames 1650–1699. We fill only enclosed missing-radiance components with at most 16 cells and a bounding box no larger than 5 × 5 cells, entirely inside valid geolocation. Each estimated pixel is solved from its immediate neighbors (discrete harmonic interpolation). Valid measurements, larger masks and footprint boundaries remain unchanged.</p><div class="controls"><button id="play">Pause</button><button id="back">◀</button><input id="slider" type="range" min="0" max="49" value="0"><button id="next">▶</button><select id="speed"><option value="160">Fast</option><option value="400">Slow</option><option value="1100">Approx. real time</option></select><b id="label"></b></div><img id="image"><p>The top row compares radiance; the bottom row compares matching texture using the same one-pixel smoothing and high-pass filter. Gray areas remain missing. Filling restores support around small holes, so the change in the texture footprint is larger than the number of estimated pixels. All color scales remain fixed through time. Loop reset is not physical motion.</p><p>Filled values are estimates, not new observations. They are saved with a separate interpolation mask. This is a separate animation experiment; the original source files and orbit products have not been replaced.</p><script>const images=__IMAGES__,times=__TIMES__;let f=0,playing=true;const img=document.getElementById('image');function draw(){img.src=images[f];document.getElementById('slider').value=f;document.getElementById('label').textContent=`Frame ${1650+f} · +${times[f].toFixed(1)} s`}function pause(){playing=false;document.getElementById('play').textContent='Play'}document.getElementById('play').onclick=()=>{playing=!playing;document.getElementById('play').textContent=playing?'Pause':'Play'};document.getElementById('slider').oninput=e=>{pause();f=Number(e.target.value);draw()};document.getElementById('back').onclick=()=>{pause();f=(f+49)%50;draw()};document.getElementById('next').onclick=()=>{pause();f=(f+1)%50;draw()};function tick(){setTimeout(()=>{if(playing){f=(f+1)%50;draw()}tick()},Number(document.getElementById('speed').value))}draw();tick();</script></body></html>'''
    html=html.replace('__IMAGES__',json.dumps(images)).replace('__TIMES__',json.dumps((times-times[0]).tolist()));(OUT/'small_gaps_filled.html').write_text(html,encoding='utf-8')
    assert Image.open(OUT/'small_gaps_filled.gif').n_frames==50
    report=dict(frames=[1650,1699],max_gap_cells=16,max_bounding_box_cells=5,connectivity=8,method='4-neighbor discrete harmonic interpolation using unchanged finite boundary samples',total_interpolated_frame_pixels=int(masks.sum()),per_frame_filled_range=[min(counts),max(counts)],total_components=len(sizes),gap_size_range=[min(sizes),max(sizes)],original_valid_texture_pixels=int(np.isfinite(before).sum()),filled_valid_texture_pixels=int(np.isfinite(after).sum()),checks='PASS: original finite measurements unchanged; filled pixels only inside valid geometry; finite footprint expands only by marked cells; interpolation within boundary range; 50 animation frames',note='Counts sum across 50 frames, not unique map cells. Filled values remain flagged and original source data are untouched.')
    savejson(OUT/'metadata.json',report);print(json.dumps(report,indent=2))

if __name__=='__main__':main()
