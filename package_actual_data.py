from pathlib import Path
import sys,json,shutil,zipfile,hashlib
ROOT=Path(__file__).parent
sys.path.insert(0,str(ROOT/'.packages'))
import numpy as np
import netCDF4

out=ROOT/'AWE_actual_data';out.mkdir(exist_ok=True)
shutil.copy2(ROOT/'parallax_demo/projected_frames.npz',out/'projected_frames.npz')
shutil.copy2(ROOT/'parallax_demo/metadata.json',out/'projection_metadata.json')
shutil.copy2(ROOT/'layer_separation/metrics.json',out/'previous_fit_metrics.json')
shutil.copy2(ROOT/'layer_separation/motion_candidates.csv',out/'motion_candidates.csv')
source=Path(r'Y:\soc\l1r\2023\11\awe_l1r_q20_2023326T0547_00004_v01.nc')
with netCDF4.Dataset(source) as src,netCDF4.Dataset(out/'original_frames_0930_0979.nc','w') as dst:
    dst.setncatts(src.__dict__)
    dst.setncattr('subset_note','Original zero-based frame indices 930 through 979; original global Date_Start/End refer to the full source file, not this subset.')
    dst.setncattr('source_path',str(source))
    for name,dim in src.dimensions.items():dst.createDimension(name,50 if name=='time' else len(dim))
    for name,v in src.variables.items():
        attrs=v.__dict__.copy();fill=attrs.pop('_FillValue',None)
        kw={'fill_value':fill} if fill is not None else {}
        w=dst.createVariable(name,v.dtype,v.dimensions,zlib=True,**kw);w.setncatts(attrs)
        slices=tuple(slice(930,980) if d=='time' else slice(None) for d in v.dimensions)
        w[:]=v[slices]
    w=dst.createVariable('source_frame_index','i4',('time',));w[:]=np.arange(930,980)

old=(ROOT/'separate_layers.py').read_text(encoding='utf-8')
functions=old[old.index('def highpass'):old.index('def render')]
functions=functions.replace("OUT/'motion_candidates.csv'","OUTPUT/'recomputed_motion_candidates.csv'")
cells=[]
def add(kind,text):
    c={'cell_type':kind,'metadata':{},'id':f'cell-{len(cells)}','source':text.splitlines(keepends=True)}
    if kind=='code':c.update(execution_count=None,outputs=[])
    cells.append(c)
add('markdown','''# AWE actual-data layer separation

Extract the **entire ZIP**, then open this notebook from the extracted folder. No network drive or original 602 MB file is required.

Included are 50 actual frames, original zero-based indices **930–979**, from `awe_l1r_q20_2023326T0547_00004_v01.nc`. `projected_frames.npz` contains their calibrated radiance on the same fixed 0.02° latitude/longitude grid used in our demonstration. `original_frames_0930_0979.nc` also preserves the original 300 × 300 radiance and geolocation arrays plus the other source variables for these frames.

This experiment uses the **87 km OH-layer projection**. Its fast/slow labels identify apparent motion, not proven physical classes. The notebook defaults to previously fitted velocities; an optional cell re-estimates them. Contrast error against ground truth is unavailable for real data. Reconstruction RMSE alone cannot establish clean separation.
''')
add('code','# Run once if needed. netCDF4 is only needed to inspect the original-frame subset.\n%pip install numpy scipy matplotlib netCDF4\n')
add('code','''from pathlib import Path
import json, warnings
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter, map_coordinates, maximum_filter
from scipy.signal import fftconvolve

DATA_FOLDER = Path('.')  # Change if your notebook working directory is elsewhere.
ITERATIONS = 64         # Set to 256, 512, etc. Longer runs may take several minutes.
REESTIMATE_VELOCITIES = False
OUTPUT = DATA_FOLDER / f'actual_results_{ITERATIONS}'
OUTPUT.mkdir(exist_ok=True)
data = np.load(DATA_FOLDER / 'projected_frames.npz')
I = data['radiance'].astype(float)
lat, lon = data['latitude'], data['longitude']
frames = data['frame_index']
times = data['epoch'] - data['epoch'][25]
ny, nx = I.shape[1:]
yy, xx = np.indices((ny, nx), dtype=float)
previous = json.loads((DATA_FOLDER / 'previous_fit_metrics.json').read_text())
va = np.array(previous['slow_velocity_pixels_per_s'])
vc = np.array(previous['fast_velocity_pixels_per_s'])
assert I.shape == (50, 309, 317)
assert frames.tolist() == list(range(930, 980))
print('Frames:', frames[0], 'to', frames[-1], '| duration:', times[-1]-times[0], 's')
print('Slow / fast velocities [east, north], grid pixels/s:', va, vc)
''')
add('markdown','''## Inspect the actual observations
Values are in kRayleigh. Gray/blank areas are outside the projected footprint. The original Epoch values are retained; only relative seconds are used, avoiding ambiguity in the source time-system labels.
''')
add('code','''extent = [lon[0]-.01, lon[-1]+.01, lat[0]-.01, lat[-1]+.01]
aspect = 1/np.cos(np.deg2rad(lat.mean()))
lo, hi = np.nanpercentile(I, [2, 98])
fig, axes = plt.subplots(1, 3, figsize=(14, 5), layout='constrained')
for ax, j in zip(axes, [0, 25, 49]):
    im = ax.imshow(I[j], origin='lower', extent=extent, aspect=aspect,
                   cmap='gray', vmin=lo, vmax=hi)
    ax.set(title=f'Frame {frames[j]}', xlabel='Longitude (degrees E)', ylabel='Latitude (degrees N)')
fig.colorbar(im, ax=axes, label='kRayleigh', shrink=.7)
plt.show()
''')
add('markdown','''## Matching and alignment functions
The following cell contains the original patch-matching implementation and interpolation functions. Optional motion re-estimation uses the same heuristic thresholds: patch width 31, lags 8/16/24, NCC > 0.35, slow speed < 0.18 pixels/s, fast speed > 0.22 pixels/s, averaging radius 0.09 pixels/s. These are not calibrated physical limits.
''')
add('code',functions)
add('code','''if REESTIMATE_VELOCITIES:
    candidates = motions()
    va, slow_count, slow_ncc = mode_velocity(candidates, slow=True)
    vc, fast_count, fast_ncc = mode_velocity(candidates, slow=False)
    print('Re-estimated velocities:', va, vc)
else:
    print('Using velocities from the previous actual-data fit.')
''')
add('markdown','''## Alternate fast and slow template updates
Each iteration uses all 50 observations. First subtract slow and estimate fast from its aligned median; then subtract fast and estimate slow. A template update requires 10 usable frames. **Where support is insufficient, retain the previous template value** to prevent the boundary erosion found in the original implementation. Such retained values are not new measurements.

The fast template is assigned median zero. Absolute brightness division between the layers remains ambiguous. This boundary correction means results need not exactly reproduce the original uncorrected eight-iteration real-data fit.

Checkpoints at powers of two and the final iteration contain templates and support counts. RMSE uses a fixed region observed in all frames and supported by the first fitted model. The notebook fails rather than silently changing that evaluation region if support is later lost.
''')
run='''assert ITERATIONS >= 1
def predict(template, velocity):
    return np.array([warp(template, velocity, t) for t in times])

slow, slow_support = stack_layer(I, va)
fast = np.zeros_like(slow)
slowseq = predict(slow, va)
history = []
fixed_mask = None
save_at = {ITERATIONS} | {2**i for i in range(ITERATIONS.bit_length()) if 2**i <= ITERATIONS}
for iteration in range(1, ITERATIONS+1):
    candidate, fast_support = stack_layer(I-slowseq, vc)
    fast = np.where(np.isfinite(candidate), candidate, fast)
    offset = np.nanmedian(fast)
    fast -= offset
    slow += offset
    fastseq = predict(fast, vc)
    candidate, slow_support = stack_layer(I-fastseq, va)
    slow = np.where(np.isfinite(candidate), candidate, slow)
    slowseq = predict(slow, va)
    residual = I-slowseq-fastseq
    if fixed_mask is None:
        fixed_mask = np.all(np.isfinite(residual), axis=0) & data['common_coverage']
        if fixed_mask.sum() == 0: raise RuntimeError('No shared evaluation region')
    if not np.isfinite(residual[:, fixed_mask]).all():
        raise RuntimeError('Fixed evaluation region lost support')
    rmse = float(np.sqrt(np.mean(residual[:, fixed_mask]**2)))
    history.append({'iteration': iteration, 'reconstruction_rmse_kR': rmse,
                    'evaluation_pixels_per_frame': int(fixed_mask.sum())})
    if iteration in save_at:
        np.savez_compressed(OUTPUT / f'checkpoint_{iteration:04d}.npz',
            slow=slow, fast=fast, slow_support=slow_support, fast_support=fast_support,
            slow_velocity=va, fast_velocity=vc, latitude=lat, longitude=lon,
            frame_index=frames, relative_seconds=times, evaluation_mask=fixed_mask)
    if iteration % 8 == 0 or iteration in save_at:
        print(f'Iteration {iteration}: reconstruction RMSE = {rmse:.5f} kRayleigh', flush=True)
        (OUTPUT / 'metrics.json').write_text(json.dumps(history, indent=2))
print('Results saved in', OUTPUT.resolve())
'''
add('code',run)
add('code','''fig, ax = plt.subplots(figsize=(8,4), layout='constrained')
ax.plot([m['iteration'] for m in history], [m['reconstruction_rmse_kR'] for m in history])
ax.set(xlabel='Completed iterations', ylabel='Reconstruction RMSE (kRayleigh)',
       title='Fixed common coverage; training fit, not independent validation')
ax.grid(alpha=.25)
fig.savefig(OUTPUT / 'error_history.png', dpi=150)
plt.show()
''')
add('code','''ITERATION_TO_VIEW = ITERATIONS  # Any saved checkpoint: 1,2,4,8,..., or final.
FRAME_TO_VIEW = 955            # Original zero-based frame index, 930 through 979.
checkpoint = np.load(OUTPUT / f'checkpoint_{ITERATION_TO_VIEW:04d}.npz')
j = int(np.flatnonzero(frames == FRAME_TO_VIEW)[0])
a = warp(checkpoint['slow'], checkpoint['slow_velocity'], times[j])
c = warp(checkpoint['fast'], checkpoint['fast_velocity'], times[j])
r = I[j]-a-c
# Fixed scales allow comparison across checkpoints.
fig, axes = plt.subplots(2,2, figsize=(12,10), layout='constrained')
for ax, arr, title, cmap, low, high in zip(axes.flat, [I[j],a,c,r],
    ['Observed', 'Slow component', 'Fast signed contrast', 'Unexplained residual'],
    ['gray','gray','RdBu_r','RdBu_r'], [lo,lo,-2,-1], [hi,hi,2,1]):
    im=ax.imshow(arr, origin='lower', extent=extent, aspect=aspect, cmap=cmap,vmin=low,vmax=high)
    ax.set(title=title, xlabel='Longitude (degrees E)', ylabel='Latitude (degrees N)')
    fig.colorbar(im,ax=ax,shrink=.75,label='kRayleigh')
fig.suptitle(f'Frame {FRAME_TO_VIEW} | iteration {ITERATION_TO_VIEW}')
fig.savefig(OUTPUT / f'view_{ITERATION_TO_VIEW:04d}_{FRAME_TO_VIEW}.png',dpi=130)
plt.show()
''')
add('markdown','''## Optional: inspect the original, unprojected subset
`time[0]` here is original frame 930; `time[25]` is original frame 955. The subset preserves original global date metadata, which refers to the full source orbit. Use `source_frame_index` and `Epoch` for the selected interval. To change the projection, start with these original arrays. The runnable separation above uses the supplied fixed-grid frames.
''')
add('code','''import netCDF4
with netCDF4.Dataset(DATA_FOLDER / 'original_frames_0930_0979.nc') as original:
    print(original.variables.keys())
    print('Original radiance shape:', original['Radiance'].shape)
    raw_frame_955 = original['Radiance'][25]
    raw_lat_955 = original['Latitude'][25]
    raw_lon_955 = original['Longitude'][25]
''')
nb={'cells':cells,'metadata':{'kernelspec':{'display_name':'Python 3','language':'python','name':'python3'},'language_info':{'name':'python','version':'3.12'}},'nbformat':4,'nbformat_minor':5}
(out/'AWE_actual_data.ipynb').write_text(json.dumps(nb,indent=2),encoding='utf-8')
(out/'requirements.txt').write_text('numpy\nscipy\nmatplotlib\nnetCDF4\n')
(out/'README.txt').write_text('''Extract this entire ZIP into a folder. Open AWE_actual_data.ipynb in Jupyter and run cells in order.
Set ITERATIONS to 64, 256, or another positive integer. Use a new output folder to preserve earlier runs.
The notebook is self-contained; it needs only the accompanying data files and the packages in requirements.txt.
Original subset: 50 frames, source zero-based indices 930-979, original radiance/coordinates and ancillary variables.
Projected data: 50 x 309 x 317 radiance values in kRayleigh, fixed 0.02-degree grid using 87 km source geolocation.
Source: awe_l1r_q20_2023326T0547_00004_v01.nc.
No synthetic data are used. Ground-truth layer errors cannot be computed for this real sequence.
The optional re-estimation uses heuristic matching thresholds. Default velocities come from the previous actual-data fit.
Boundary handling retains previous template values when an update has insufficient support.
The PNG panels and residual-fit metrics are exploratory; do not assume recovered slow radiance is pure airglow.
''',encoding='utf-8')
print('Created notebook and original-frame subset.',flush=True)
