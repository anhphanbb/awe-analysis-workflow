"""Self-contained illustrated walkthrough of the existing separation run."""
import separate_layers as s
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from pathlib import Path
import json, base64

out=s.OUT/'walkthrough_assets';out.mkdir(exist_ok=True)
metrics=json.loads((s.OUT/'metrics.json').read_text())
va=np.array(metrics['slow_velocity_pixels_per_s']);vc=np.array(metrics['fast_velocity_pixels_per_s'])
rows=np.loadtxt(s.OUT/'motion_candidates.csv',delimiter=',',skiprows=1)
H=np.array([s.highpass(a) for a in s.I])
plt.rcParams.update({'font.size':10,'figure.facecolor':'#f6f8fc','axes.facecolor':'#dce1e9'})
examples=[]
for name,v in [('Fast',vc),('Slow',va)]:
    candidates=rows[(rows[:,2]==25)&(rows[:,3]==41)&(rows[:,7]==0)]
    near=candidates[np.linalg.norm(candidates[:,4:6]-v,axis=1)<.07]
    row=near[np.argmax(near[:,6])]
    cx,cy,start,end=map(int,row[:4]);r=18;h=15
    template=H[start,cy-h:cy+h+1,cx-h:cx+h+1]
    search=H[end,cy-h-r:cy+h+r+1,cx-h-r:cx+h+r+1]
    corr=s.ncc(template,search)
    py,px=np.unravel_index(np.nanargmax(corr),corr.shape)
    dx,dy=px-r,py-r
    match=search[py:py+31,px:px+31]
    direct=np.corrcoef(template.ravel(),match.ravel())[0,1]
    assert abs(direct-corr[py,px])<1e-9
    extent=[s.lon[0]-.01,s.lon[-1]+.01,s.lat[0]-.01,s.lat[-1]+.01]
    fig,axes=plt.subplots(1,2,figsize=(11,5.5),layout='constrained')
    for ax,j in zip(axes,[start,end]):
        ax.imshow(s.I[j],origin='lower',extent=extent,cmap='gray',vmin=6.6,vmax=12.4,aspect=1/np.cos(np.deg2rad(s.lat.mean())))
        ax.add_patch(Rectangle((s.lon[cx]-31*.02/2,s.lat[cy]-31*.02/2),31*.02,31*.02,fill=False,edgecolor='#ffd34d',lw=2,label='Original patch position'))
        if j==end:
            ax.add_patch(Rectangle((s.lon[cx]-(31+2*r)*.02/2,s.lat[cy]-(31+2*r)*.02/2),(31+2*r)*.02,(31+2*r)*.02,fill=False,edgecolor='#39cce0',lw=1.5,label='67 × 67 search area'))
            ax.add_patch(Rectangle((s.lon[cx+dx]-31*.02/2,s.lat[cy+dy]-31*.02/2),31*.02,31*.02,fill=False,edgecolor='#ff6464',lw=2,label='Best integer match'))
        ax.set(title=f'{name} example · frame {s.frames[j]}',xlabel='Longitude (°E)',ylabel='Latitude (°N)');ax.legend(fontsize=8,loc='upper left')
    fig.savefig(out/f'{name.lower()}_location.png',dpi=120,bbox_inches='tight');plt.close(fig)
    def serial(a):return np.where(np.isfinite(a),a,0).tolist()
    def nullable(a):return [[float(v) if np.isfinite(v) else None for v in line] for line in a]
    examples.append({'name':name,'cx':cx,'cy':cy,'lon':float(s.lon[cx]),'lat':float(s.lat[cy]),'dt':float(s.times[end]-s.times[start]),'template':nullable(template),'search':nullable(search),'corr':nullable(corr),'best':[int(dx),int(dy)],'subpixel':[float(row[4]*(s.times[end]-s.times[start])),float(row[5]*(s.times[end]-s.times[start]))],'velocity':row[4:6].tolist(),'score':float(corr[py,px]),'limit':float(np.nanpercentile(np.abs(np.concatenate([template.ravel(),search.ravel()])),98))})

# Capture the actual first alternating update, not a simulated illustration.
initial,count=s.stack_layer(s.I,va)
initialseq=np.array([s.warp(initial,va,t) for t in s.times])
res=s.I-initialseq
firstfast,_=s.stack_layer(res,vc)
offset=np.nanmedian(firstfast);firstfast-=offset
firstseq=np.array([s.warp(firstfast,vc,t) for t in s.times])
updated,_=s.stack_layer(s.I-firstseq,va)
final=np.load(s.OUT/'separated_layers.npz')
extent=[s.lon[0]-.01,s.lon[-1]+.01,s.lat[0]-.01,s.lat[-1]+.01]
def panel(ax,a,title,contrast=False):
    im=ax.imshow(a,origin='lower',extent=extent,cmap='RdBu_r' if contrast else 'gray',vmin=-2 if contrast else 6.6,vmax=2 if contrast else 12.4,aspect=1/np.cos(np.deg2rad(s.lat.mean())))
    ax.set(title=title,xlabel='Longitude (°E)',ylabel='Latitude (°N)')
    return im
fig,axes=plt.subplots(2,3,figsize=(14,9),layout='constrained')
for ax,a,title,contrast in zip(axes.flat,[s.I[25],initial,res[25],firstfast,updated,final['slow'][25]],['1. Observed frame 955','2. Initial slow-aligned median','3. Observation − initial slow','4. First fast-aligned residual median','5. Slow template after one update','6. Slow template after eight updates'],[False,False,True,True,False,False]):
    im=panel(ax,a,title,contrast);fig.colorbar(im,ax=ax,shrink=.7,label='kRayleigh')
fig.savefig(out/'iterations.png',dpi=115,bbox_inches='tight');plt.close(fig)

# See the same fast feature before/after alignment, using residuals after initial slow removal.
ex=examples[0];cx=ex['cx'];cy=ex['cy'];jlist=[9,25,41];half=28
fig,axes=plt.subplots(2,3,figsize=(12,7),layout='constrained')
for col,j in enumerate(jlist):
    aligned=s.warp(res[j],vc,-s.times[j])
    for rr,a in enumerate([res[j],aligned]):
        crop=a[cy-half:cy+half+1,cx-half:cx+half+1]
        axes[rr,col].imshow(crop,origin='lower',cmap='RdBu_r',vmin=-2,vmax=2,extent=[-half-.5,half+.5,-half-.5,half+.5])
        axes[rr,col].axhline(0,color='#313131',alpha=.4,lw=.7);axes[rr,col].axvline(0,color='#313131',alpha=.4,lw=.7)
        axes[rr,col].set(title=f'Frame {s.frames[j]} · '+('fixed map' if rr==0 else 'shifted back to frame 955'),xlabel='Eastward pixel offset',ylabel='Northward pixel offset')
fig.suptitle('Align the fast-moving residual before taking its temporal median\nTop: fixed geographic crop. Bottom: same crop after undoing fast motion. All panels share ±2 kRayleigh.')
fig.savefig(out/'alignment.png',dpi=120,bbox_inches='tight');plt.close(fig)

def uri(p):return 'data:image/png;base64,'+base64.b64encode(Path(p).read_bytes()).decode()
images={name:uri(path) for name,path in {'preprocess':s.OUT/'matching_background_0955.png','fast':out/'fast_location.png','slow':out/'slow_location.png','motion':s.OUT/'motion_diagnostics.png','alignment':out/'alignment.png','iterations':out/'iterations.png','final':s.OUT/'preview.png','mask':s.OUT/'confidence.png'}.items()}
html=r'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>AWE — how the layers were separated</title><style>
body{font:17px/1.6 system-ui,sans-serif;background:#eef2f7;color:#233148;margin:0}main{max-width:1180px;margin:auto;padding:32px 24px 80px}h1{font-size:36px;line-height:1.2}h2{font-size:25px}h3{font-size:19px}section{background:white;border:1px solid #dce3ec;border-radius:16px;padding:25px;margin:24px 0}p{max-width:1000px}img{width:100%;height:auto;border-radius:8px}.eyebrow{color:#426c9e;letter-spacing:2px;font-size:13px;font-weight:bold}.note{background:#eef5ff;border-left:4px solid #426c9e;padding:12px 18px}button,select{font:inherit;padding:7px 13px;border-radius:7px;border:1px solid #9bacc1;background:white;cursor:pointer}button.active{background:#214d7e;color:white}.grid{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:20px}.tile{text-align:center;font-size:14px}canvas{width:100%;max-width:340px;aspect-ratio:1;image-rendering:pixelated;border:1px solid #acb8ca;background:#dce1e9}input[type=range]{width:230px}label{display:inline-block;margin:10px}table{border-collapse:collapse;width:100%}td,th{text-align:left;border-bottom:1px solid #dce3ec;padding:10px}.formula{font:16px/1.7 ui-monospace,monospace;background:#f1f4f8;padding:16px;overflow:auto}a{color:#245d98}nav a{margin-right:15px}#readout{font-variant-numeric:tabular-nums;background:#eef5ff;padding:14px}small{color:#53647b}@media(max-width:700px){.grid{grid-template-columns:1fr}h1{font-size:28px}section{padding:15px}}</style></head><body><main>
<div class="eyebrow">AWE Q20 · WORKED EXAMPLE · 50 FRAMES</div><h1>From a moving patch to two image layers</h1><p>A walkthrough of the actual calculation on frames 930–979. Explore one fast-moving patch and one slow-moving patch, move the comparison window yourself, then follow the reconstruction.</p><p class="note">This entire page is self-contained. Copy <b>walkthrough.html</b> to another computer to keep the figures and interactive patch comparison.</p><nav><a href="#pre">1. Preprocess</a><a href="#patch">2. Compare patches</a><a href="#motion">3. Measure velocities</a><a href="#align">4. Align</a><a href="#fit">5. Separate</a><a href="#check">6. Check</a></nav>
<section id="pre"><h2>1 · Make local texture easier to match</h2><p>First project every frame onto the same 0.02° grid using the file’s 87 km geolocation. Then subtract a Gaussian smooth with σ = 5 grid pixels. This preprocessing is used only to estimate motion; the layer reconstruction uses the original projected radiance.</p><img src="@@preprocess@@"><div class="formula">background = Gaussian(valid × radiance) / Gaussian(valid)<br>matching texture = radiance − background</div><p>Missing pixels contribute zero weight. Texture pixels are excluded if their Gaussian support is ≤98%. The smooth background contains both layers; it is not an airglow estimate.</p></section>
<section id="patch"><h2>2 · Pick a patch and search for it in a later frame</h2><p>Both examples use <b>frame 955 → frame 971</b>, a 16-frame separation. The yellow box is the original 31 × 31 pixel patch. The cyan box is the later frame’s 67 × 67 pixel search area; the red box marks the best integer-pixel match.</p><p><button id="fastButton" class="active">Fast-motion example</button> <button id="slowButton">Slow-motion example</button></p><p id="exampleInfo"></p><img id="location" src="@@fast@@"><p>For this lag the search radius is 18 pixels in each direction. A 31-pixel window slides inside the 67-pixel search area, producing <b>37 × 37 = 1,369 candidate positions</b>. Positions containing missing data are excluded.</p>
<div class="grid"><div class="tile"><h3>Reference texture · frame 955</h3><canvas id="template" width="310" height="310"></canvas></div><div class="tile"><h3>Search area · frame 971</h3><canvas id="search" width="335" height="335"></canvas><small>White dashed: zero motion. Yellow: current candidate.</small></div><div class="tile"><h3>Current candidate · 31 × 31</h3><canvas id="candidate" width="310" height="310"></canvas></div></div>
<p><label>Eastward shift Δx <input id="dx" type="range" min="-18" max="18" value="0"></label><label>Northward shift Δy <input id="dy" type="range" min="-18" max="18" value="0"></label><button id="zero">Zero displacement</button> <button id="best">Best match</button> <button id="scan">Scan one row</button></p><div id="readout"></div>
<div class="grid"><div class="tile"><h3>Normalized correlation surface</h3><canvas id="correlation" width="370" height="370"></canvas><small>Click to select a candidate. Right = east (+Δx); up = north (+Δy). Dark blue: −1, white: 0, red: +1. Gray: invalid. Black circle: best match; yellow square: current.</small></div><div style="grid-column:span 2"><h3>What “compare” means</h3><p>For each candidate, subtract the mean from each patch, multiply corresponding pixel values, sum the products, then divide by the two patches’ combined contrast:</p><div class="formula">NCC = Σ[(T − mean(T)) × (P − mean(P))]<br>      / √{Σ(T − mean(T))² × Σ(P − mean(P))²}</div><p>A score near +1 means similar texture even if overall brightness or contrast differs. Zero means little linear agreement; negative scores indicate reversed contrast. The original and candidate displays share a color scale, although NCC itself normalizes their means and contrast.</p></div></div>
<p id="subpixel"></p><p class="note">These examples were selected for high correlation near the measured fast and slow velocity modes. They demonstrate successful matches, not the quality of every patch. A single patch is not assumed to contain two separately measurable layers.</p></section>
<section id="motion"><h2>3 · Convert matches into velocities, then look for groups</h2><div class="formula">vx = Δx / elapsed seconds &nbsp; vy = Δy / elapsed seconds</div><p>The full run repeats matching at 8-, 16-, and 24-frame lags, with starting sequence indices 0, 12, and 25. Patch centers are spaced 18 pixels apart. Texture standard deviation must exceed 0.08 kRayleigh. Up to four distinct local correlation maxima above 0.35 are kept; search-boundary maxima are rejected.</p><img src="@@motion@@"><p>The plotted points are the strongest matches. The velocity estimator uses all retained peaks, with weight <code>(NCC − 0.35)² / (1 + peak rank)</code>. It finds a smoothed histogram mode separately below 0.18 pixels/s and above 0.22 pixels/s, then averages candidates within 0.09 pixels/s of that mode. These speed ranges are explicit modeling choices.</p><table><tr><th>Fitted layer</th><th>East (pixels/s)</th><th>North (pixels/s)</th><th>Apparent speed</th></tr><tr><td>Slow</td><td>+0.1269</td><td>+0.0754</td><td>0.32 km/s</td></tr><tr><td>Fast</td><td>+0.4677</td><td>−0.3970</td><td>1.34 km/s</td></tr></table><p>These are apparent map velocities, not physical winds. The reconstruction uses one constant velocity per layer, not a separate velocity for every patch.</p></section>
<section id="align"><h2>4 · Undo a layer’s motion so its features line up</h2><p>To place a frame at the reference time of frame 955, shift it by <b>−velocity × time difference</b>. If the velocity matches the pattern, that pattern becomes stationary. Taking a median through the aligned stack then preserves persistent features while reducing structures moving differently.</p><img src="@@alignment@@"><p>This example uses observations after subtracting the initial slow estimate. All crops cover the same reference pixel coordinates. Missing support remains missing; interpolation does not wrap at an edge. Imperfect alignment shows where the constant-velocity model is only approximate.</p></section>
<section id="fit"><h2>5 · Alternate between the two layers</h2><div class="formula">Observed(t) ≈ shift(Slow, vSlow × t) + shift(Fast, vFast × t) + Residual(t)</div><ol><li>Align original radiance with the slow velocity and take a temporal median: the initial slow template.</li><li>Subtract that slow prediction from every observed frame.</li><li>Align the residuals with the fast velocity and take their temporal median: the first fast template.</li><li>Subtract the fast prediction from the observations, then re-estimate the slow template.</li><li>Repeat the updates eight times. A template pixel requires at least 10 contributing frames.</li></ol><img src="@@iterations@@"><p>The six panels above show actual intermediate arrays from the first update and the final saved fit. For display, the first fast template uses the same zero-median convention as the algorithm.</p><p class="note">Absolute brightness cannot be uniquely split: adding a constant to one template and subtracting it from the other gives the same sum. The fast template is assigned median zero. Broad or persistent overlapping structure can remain mixed between layers.</p><img src="@@final@@"><p>The residual panel is the observation minus both predicted components. It shows unmodeled texture, brightness changes, velocity variation, and noise—not an additional identified physical layer.</p></section>
<section id="check"><h2>6 · Make an exploratory mask and check prediction</h2><img src="@@mask@@"><p>The mask score is <code>clip(1 − local error energy / local energy after slow subtraction, 0, 1)</code>. Both energies are Gaussian-averaged with σ = 3 pixels. Scores require >95% support and sufficient residual energy (>0.04 in the implemented weighted calculation). Scores above 0.6 produce the candidate mask. This measures local explanatory power, not cloud probability.</p><table><tr><th>Prediction on withheld odd sequence indices</th><th>RMSE (kRayleigh)</th></tr><tr><td>Slow-only template</td><td>0.669</td></tr><tr><td>Two-layer templates</td><td>0.423</td></tr></table><p>Templates are fitted to even sequence indices and evaluated at odd indices on the same valid pixels (1,074,100 samples). The predictive two-layer fit uses five updates. Velocities were estimated using all frames, so this is a template-prediction check, not fully independent validation. The ~37% improvement supports the two-motion model but does not establish a pure airglow retrieval.</p><p>Numerical checks also verify known-shift recovery, warp direction, edge handling, and that slow + fast + residual reconstructs the observation wherever all are valid.</p></section>
</main><script>
const examples=@@DATA@@,locations={Fast:'@@fast@@',Slow:'@@slow@@'};let e=examples[0],timer=null;const $=x=>document.getElementById(x);
function color(v,lim){if(v===null||!Number.isFinite(v))return '#dce1e9';let z=Math.max(-1,Math.min(1,v/lim));let q=Math.round(245-185*Math.abs(z));return z>=0?`rgb(235,${q},${q})`:`rgb(${q},${q},235)`}
function draw(id,a,lim){let c=$(id),ctx=c.getContext('2d'),n=a.length,m=a[0].length;ctx.clearRect(0,0,c.width,c.height);for(let y=0;y<n;y++)for(let x=0;x<m;x++){ctx.fillStyle=color(a[y][x],lim);ctx.fillRect(x*c.width/m,(n-1-y)*c.height/n,Math.ceil(c.width/m),Math.ceil(c.height/n))}}
function box(ctx,x,y,n,size,stroke,dash=[]){ctx.strokeStyle=stroke;ctx.lineWidth=2;ctx.setLineDash(dash);ctx.strokeRect(x*size/67,(67-y-n)*size/67,n*size/67,n*size/67);ctx.setLineDash([])}
function update(){let dx=+$('dx').value,dy=+$('dy').value,px=dx+18,py=dy+18;let candidate=e.search.slice(py,py+31).map(row=>row.slice(px,px+31));draw('template',e.template,e.limit);draw('search',e.search,e.limit);draw('candidate',candidate,e.limit);let ctx=$('search').getContext('2d');box(ctx,18,18,31,335,'white',[5,3]);box(ctx,px,py,31,335,'#ffe000');draw('correlation',e.corr,1);ctx=$('correlation').getContext('2d');ctx.strokeStyle='#ffe000';ctx.lineWidth=3;ctx.strokeRect(px*10,(36-py)*10,10,10);ctx.strokeStyle='black';ctx.beginPath();ctx.arc((e.best[0]+18.5)*10,(18.5-e.best[1])*10,7,0,Math.PI*2);ctx.stroke();let score=e.corr[py][px];$('readout').textContent=`Δx = ${dx} pixels east, Δy = ${dy} pixels north · NCC = ${score===null?'invalid (missing data)':score.toFixed(4)} · integer-candidate velocity = ${(dx/e.dt).toFixed(3)}, ${(dy/e.dt).toFixed(3)} pixels/s`;}
function stop(){clearInterval(timer);timer=null;$('scan').textContent='Scan one row'}function select(i){stop();e=examples[i];$('fastButton').classList.toggle('active',i===0);$('slowButton').classList.toggle('active',i===1);$('location').src=locations[e.name];$('exampleInfo').textContent=`${e.name} example: center (${e.lon.toFixed(2)}°E, ${e.lat.toFixed(2)}°N), grid pixel (${e.cx}, ${e.cy}); elapsed time ${e.dt.toFixed(4)} s.`;$('subpixel').textContent=`After finding the integer peak (${e.best[0]}, ${e.best[1]}), the code fits a three-point parabola along each axis, with refinement capped at ±0.5 pixel. This example gives (${e.subpixel[0].toFixed(3)}, ${e.subpixel[1].toFixed(3)}) pixels, or (${e.velocity[0].toFixed(4)}, ${e.velocity[1].toFixed(4)}) pixels/s. The retained score ${e.score.toFixed(4)} is the integer-peak NCC; the refinement is not a new correlation measurement.`;$('dx').value=0;$('dy').value=0;update()}
$('fastButton').onclick=()=>select(0);$('slowButton').onclick=()=>select(1);for(let id of ['dx','dy'])$(id).oninput=()=>{stop();update()};$('zero').onclick=()=>{stop();$('dx').value=0;$('dy').value=0;update()};$('best').onclick=()=>{stop();$('dx').value=e.best[0];$('dy').value=e.best[1];update()};$('scan').onclick=()=>{if(timer){stop();return}$('dy').value=e.best[1];$('dx').value=-18;update();$('scan').textContent='Stop scan';timer=setInterval(()=>{let v=+$('dx').value;if(v>=18){stop();return}$('dx').value=v+1;update()},130)};$('correlation').onclick=ev=>{stop();let rect=$('correlation').getBoundingClientRect();$('dx').value=Math.min(36,Math.max(0,Math.floor((ev.clientX-rect.left)/rect.width*37)))-18;$('dy').value=18-Math.min(36,Math.max(0,Math.floor((ev.clientY-rect.top)/rect.height*37)));update()};select(0);
</script></body></html>'''
html=html.replace('@@DATA@@',json.dumps(examples,allow_nan=False))
for key,value in images.items():html=html.replace('@@'+key+'@@',value)
assert '@@' not in html
(s.OUT/'walkthrough.html').write_text(html,encoding='utf-8')
(out/'examples.json').write_text(json.dumps(examples,indent=2,allow_nan=False))
p=s.OUT/'index.html';index=p.read_text(encoding='utf-8')
if 'walkthrough.html' not in index:index=index.replace('<body>','<body><p><a href="walkthrough.html">Open the detailed interactive walkthrough: follow a patch through every step →</a></p>',1);p.write_text(index,encoding='utf-8')
print(json.dumps([{k:e[k] for k in ['name','cx','cy','dt','best','subpixel','score']} for e in examples]),flush=True)
print('Saved self-contained walkthrough:',s.OUT/'walkthrough.html',flush=True)

