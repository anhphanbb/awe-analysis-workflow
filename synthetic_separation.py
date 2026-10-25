"""Known-truth demonstration of the existing alternating-median algorithm."""
import separate_layers as s
import numpy as np
import matplotlib.pyplot as plt
import json, base64
from PIL import Image

out=s.ROOT/'synthetic_demo';out.mkdir(exist_ok=True)
# Reuse the exact warp and stack_layer functions with a smaller synthetic grid.
s.ny=s.nx=160;s.yy,s.xx=np.indices((160,160),dtype=float)
s.times=(np.arange(50)-25)*1.1
va=np.array([.12,.07]);vc=np.array([.48,-.40])
x,y=s.xx,s.yy
def air(x,y):
    return 8+.65*np.sin(2*np.pi*(x+.40*y)/23)+.30*np.cos(2*np.pi*(y-.25*x)/39)
def cloud(x,y):
    return (2.5*np.exp(-((x-73)**2/100+(y-88)**2/150))
            +1.8*np.exp(-((x-94)**2/65+(y-68)**2/100))
            +1.3*np.exp(-((x-60)**2/55+(y-65)**2/55))
            +.9*np.exp(-((x-82)**2/600+(y-80)**2/420)))
rng=np.random.default_rng(42)
truthA=air(x,y);truthC=cloud(x,y)
truthseqA=np.array([air(x-va[0]*t,y-va[1]*t) for t in s.times])
truthseqC=np.array([cloud(x-vc[0]*t,y-vc[1]*t) for t in s.times])
observed=truthseqA+truthseqC+rng.normal(0,.06,truthseqA.shape)
# Generate observations analytically, independently of reconstruction interpolation.
slow,_=s.stack_layer(observed,va)
snapshots={0:(slow.copy(),np.zeros_like(slow))};stages={};metrics=[]
for iteration in range(1,9):
    slowseq=np.array([s.warp(slow,va,t) for t in s.times])
    before=observed-slowseq
    fast,_=s.stack_layer(before,vc)
    offset=np.nanmedian(fast);fast-=offset;slow+=offset
    fastseq=np.array([s.warp(fast,vc,t) for t in s.times])
    slow,_=s.stack_layer(observed-fastseq,va)
    if iteration in (1,2,4,8):
        snapshots[iteration]=(slow.copy(),fast.copy())
        stages[iteration]=(before[25].copy(),fast.copy(),(observed-fastseq)[25].copy(),slow.copy())

roi=np.s_[40:120,40:120]
plt.rcParams.update({'font.size':10,'figure.facecolor':'#f4f7fb','axes.facecolor':'#dce1e8'})
def plot(ax,a,title,kind):
    settings={'air':('gray',6.9,9.1),'cloud':('magma',-.4,4),'obs':('gray',6.9,12),'error':('RdBu_r',-1.2,1.2)}
    cmap,lo,hi=settings[kind]
    im=ax.imshow(a[roi],origin='lower',extent=[40,120,40,120],cmap=cmap,vmin=lo,vmax=hi,interpolation='nearest')
    ax.set(title=title,xlabel='x (pixels)',ylabel='y (pixels)')
    return im
fig,axes=plt.subplots(1,3,figsize=(12,4.5),layout='constrained')
for ax,a,title,kind in zip(axes,[truthA,truthC,observed[25]],['True slow layer: wave bands','True fast layer: cloud-like blobs','Observed: slow + fast + noise'],['air','cloud','obs']):
    fig.colorbar(plot(ax,a,title,kind),ax=ax,shrink=.75,label='Synthetic radiance units')
fig.savefig(out/'truth.png',dpi=130,bbox_inches='tight');plt.close(fig)

fig,axes=plt.subplots(5,4,figsize=(13,15),layout='constrained')
for row,(iteration,(a,c)) in enumerate(snapshots.items()):
    # Shared constant offset is unidentifiable. Score contrast after one equal-and-opposite correction.
    delta=np.mean((a-truthA)[roi]);aa=a-delta;cc=c+delta
    assert np.isfinite(a[roi]).all() and np.isfinite(c[roi]).all()
    model=np.array([s.warp(a,va,t)+s.warp(c,vc,t) for t in s.times])
    residual=observed-model
    record={'iteration':iteration,'slow_contrast_rmse':float(np.sqrt(np.mean((aa-truthA)[roi]**2))),'fast_contrast_rmse':float(np.sqrt(np.mean((cc-truthC)[roi]**2))),'observation_rmse':float(np.sqrt(np.mean(residual[:,40:120,40:120]**2))),'shared_offset_for_truth_comparison':float(delta)}
    metrics.append(record)
    label='Initialization' if iteration==0 else f'Iteration {iteration}'
    for ax,arr,title,kind in zip(axes[row],[a,c,aa-truthA,cc-truthC],[f'{label} · slow',f'{label} · fast',f'Slow error · RMSE {record["slow_contrast_rmse"]:.3f}',f'Fast error · RMSE {record["fast_contrast_rmse"]:.3f}'],['air','cloud','error','error']):
        plot(ax,arr,title,kind)
    f,axs=plt.subplots(2,3,figsize=(12,8),layout='constrained')
    for ax,arr,title,kind in zip(axs.flat,[truthA,a,aa-truthA,truthC,c,cc-truthC],['True slow layer',f'{label}: estimated slow','Slow contrast error','True fast layer',f'{label}: estimated fast','Fast contrast error'],['air','air','error','cloud','cloud','error']):
        f.colorbar(plot(ax,arr,title,kind),ax=ax,shrink=.75)
    f.suptitle(f'{label} · same 50 observations reused · known velocities fixed')
    f.savefig(out/f'iteration_{iteration}.png',dpi=115,bbox_inches='tight');plt.close(f)
fig.suptitle('Alternating separation: initialization and iterations 1, 2, 4, 8\nFixed display scales: slow 6.9–9.1; fast −0.4–4; errors ±1.2. Errors compare contrast after shared-offset correction.',fontsize=13)
fig.savefig(out/'iteration_comparison.png',dpi=125,bbox_inches='tight');plt.close(fig)
f,axs=plt.subplots(1,4,figsize=(15,4.5),layout='constrained')
for ax,arr,title,kind in zip(axs,stages[1],['Subtract initial slow\n(fast + remaining slow errors)','Align fast and take median\n(first fast template)','Subtract first fast\n(slow + remaining fast errors)','Align slow and take median\n(end of iteration 1)'],['cloud','cloud','air','air']):
    f.colorbar(plot(ax,arr,title,kind),ax=ax,shrink=.7)
f.savefig(out/'first_update.png',dpi=120,bbox_inches='tight');plt.close(f)
f,ax=plt.subplots(figsize=(8,4),layout='constrained')
for key,label in [('slow_contrast_rmse','Slow contrast error'),('fast_contrast_rmse','Fast contrast error'),('observation_rmse','Observation reconstruction error')]:
    ax.plot([m['iteration'] for m in metrics],[m[key] for m in metrics],'o-',label=label)
ax.axhline(.06,color='gray',ls='--',label='Input noise σ (reference only)');ax.legend();ax.set(xlabel='Completed iterations',ylabel='RMSE (synthetic units)',xticks=[0,1,2,4,8],title='Same fixed central region for every comparison');ax.grid(alpha=.2)
f.savefig(out/'errors.png',dpi=130,bbox_inches='tight');plt.close(f)
np.savez_compressed(out/'synthetic_data.npz',observed=observed,true_slow=truthA,true_fast=truthC,times=s.times,slow_velocity=va,fast_velocity=vc,iterations=np.array(list(snapshots)),estimated_slow=np.array([v[0] for v in snapshots.values()]),estimated_fast=np.array([v[1] for v in snapshots.values()]))
(out/'metrics.json').write_text(json.dumps(metrics,indent=2))
def embed(path):return 'data:image/png;base64,'+base64.b64encode(path.read_bytes()).decode()
html='''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>Synthetic two-layer separation</title><style>body{font:17px/1.6 system-ui;background:#edf2f8;color:#213149;max-width:1180px;margin:30px auto;padding:0 22px}section{background:white;border-radius:14px;padding:22px;margin:24px 0}img{width:100%;height:auto}button{font:inherit;padding:9px 16px;margin:5px;border:1px solid #99a9bd;border-radius:7px;cursor:pointer}.active{background:#1b527e;color:white}.note{border-left:4px solid #4480ad;padding:12px;background:#edf5ff}h1{line-height:1.2}code{font-size:15px}table{border-collapse:collapse;width:100%}td,th{padding:9px;border-bottom:1px solid #ddd;text-align:left}</style></head><body><h1>Watch two synthetic layers separate</h1><p>We know the true images and their velocities. This experiment isolates Step 5: alternating template estimation. No patch matching or velocity estimation is involved.</p><section><h2>1. Build the observation</h2><img src="@@truth@@"><p>The slow layer is a pair of wave patterns. The fast layer contains three compact blobs and one broader glow. We generate 50 frames at 1.1-second spacing, shifting each pattern analytically and adding independent Gaussian noise with σ = 0.06 synthetic radiance units.</p><p>Slow velocity: (0.12, 0.07) pixels/s. Fast velocity: (0.48, −0.40) pixels/s. The same known velocities are supplied to the reconstruction. These are illustrative images, not a physical cloud/airglow simulation.</p></section><section><h2>2. Start with a contaminated slow estimate</h2><p>Undo the slow motion in all 50 observations, then take a temporal median. The wave bands line up. The blobs move across the aligned stack, but the broad glow remains at some locations for enough frames to contaminate the median. Initialize the fast estimate to zero.</p></section><section><h2>3. One iteration means two updates</h2><ol><li>Subtract the current slow prediction from each observed frame.</li><li>Undo fast motion in those residuals and take their median to update the fast template.</li><li>Subtract this updated fast prediction from the original observations.</li><li>Undo slow motion and take their median to update the slow template.</li></ol><img src="@@first@@"><p>Every iteration reuses the same 50 observations; it does not add new frames. Templates are referenced to the middle frame. The panels show that reference time, while each median uses the full sequence.</p></section><section><h2>4. Compare iterations with the known truth</h2><p>Choose a completed iteration. Each selection shows true layer → estimated layer → error. Display scales remain fixed. The error panels account for the unidentifiable constant split between layers; estimates themselves remain in the algorithm’s zero-median fast convention.</p><div id="buttons"></div><img id="iteration" src="@@initial@@"><p id="numbers"></p><p class="note">The central 80 × 80 pixels are shown and scored at every iteration, avoiding edge-support changes. The algorithm uses the same missing-data-aware warping and temporal medians as the real-data code, with at least 10 contributing frames per template pixel.</p></section><section><h2>5. All requested iterations together</h2><img src="@@comparison@@"><img src="@@errors@@"><p>Interpret both kinds of error: a low observation reconstruction error does not guarantee that the individual layers are correct. Broad structures can remain shared between templates even when their sum closely matches the observations.</p><p>The absolute constant brightness split is not identifiable. For truth-error calculations only, we subtract the mean slow-layer error from the estimated slow layer and add that same constant to the fast layer. This preserves their sum and compares spatial contrast fairly. No spatially varying correction is applied.</p><p>This is a favorable test with known, constant velocities and patterns that do not evolve. Real-data errors can be larger. Eight iterations is a demonstration choice, not a convergence guarantee.</p></section><p>This HTML is self-contained and can be copied to another computer.</p><script>const pictures=@@PICS@@,metrics=@@METRICS@@;const buttons=document.getElementById('buttons');function select(n){document.getElementById('iteration').src=pictures[n];for(const b of buttons.children)b.classList.toggle('active',+b.dataset.n===n);let m=metrics.find(x=>x.iteration===n);document.getElementById('numbers').textContent=`${n===0?'Initialization':'Iteration '+n}: slow contrast RMSE ${m.slow_contrast_rmse.toFixed(3)}; fast contrast RMSE ${m.fast_contrast_rmse.toFixed(3)}; observation RMSE ${m.observation_rmse.toFixed(3)}.`}for(const n of [0,1,2,4,8]){let b=document.createElement('button');b.textContent=n===0?'Initialization':`Iteration ${n}`;b.dataset.n=n;b.onclick=()=>select(n);buttons.appendChild(b)}select(0);</script></body></html>'''
for key,file in [('truth','truth.png'),('first','first_update.png'),('initial','iteration_0.png'),('comparison','iteration_comparison.png'),('errors','errors.png')]:html=html.replace('@@'+key+'@@',embed(out/file))
html=html.replace('@@PICS@@',json.dumps({i:embed(out/f'iteration_{i}.png') for i in snapshots})).replace('@@METRICS@@',json.dumps(metrics))
assert '@@' not in html
(out/'synthetic_walkthrough.html').write_text(html,encoding='utf-8')
print(json.dumps(metrics,indent=2))
