from pathlib import Path
import json,zipfile,hashlib
import separate_layers as s
import numpy as np
import netCDF4
out=s.ROOT/'AWE_v02_actual_data'
d=np.load(out/'projected_frames.npz')
s.I=d['radiance'].astype(float);s.lat=d['latitude'];s.lon=d['longitude'];s.frames=d['frame_index']
s.times=d['epoch']-d['epoch'][np.flatnonzero(d['matching_v01_frame_index']==955)[0]]
s.ny,s.nx=s.I.shape[1:];s.yy,s.xx=np.indices((s.ny,s.nx),dtype=float);s.OUT=out
rows=s.motions();va,na,ca=s.mode_velocity(rows,True);vc,nc,cc=s.mode_velocity(rows,False)
metrics={'slow_velocity_pixels_per_s':va.tolist(),'fast_velocity_pixels_per_s':vc.tolist(),'slow_peak_support':na,'fast_peak_support':nc,'slow_median_ncc':ca,'fast_median_ncc':cc,'source':'Re-estimated using only matched v02 frames and their actual Epoch differences'}
(out/'previous_fit_metrics.json').write_text(json.dumps(metrics,indent=2))
(out/'README.txt').write_text('''AWE Q20 orbit 4 v02, exact Epoch match to previous demonstration.
Source: E:\\soc\\l1r\\2023\\11\\awe_l1r_q20_2023326T0542_00004_v02.nc
49 exact matches to earlier v01 indices 930-979. v01 frame 958 has no exact v02 Epoch match and is omitted.
No replacement, duplicate, interpolation in time, or synthetic frame is inserted.
All 49 matches have exactly equal float64 Epoch values. See epoch_mapping.csv.
v02 source indices are 1154-1202. Reference v01 frame 955 corresponds to v02 index 1179.
v02 geolocation assumes 85 km, not the earlier 87 km. The same geographic output grid is used.
Missing radiance vertices are retained during spatial interpolation to preserve v02 masks.
Extract the entire ZIP. Open AWE_v02_actual_data.ipynb in Jupyter and run cells in order.
Set ITERATIONS to 256 or another count. Default is 64.
Default velocities have been freshly estimated from these v02 frames, not reused from v01.
The optional original-frame section reads the supplied subset; no E: drive is required locally.
Original global date metadata refer to the full source file; Epoch and mapping specify subset times.
Reconstruction residual is not ground-truth separation error. More iterations need not yield pure layers.
''',encoding='utf-8')
old=np.load(s.ROOT/'parallax_demo/projected_frames.npz')
np.testing.assert_equal(d['epoch'],old['epoch'][np.isin(old['frame_index'],d['matching_v01_frame_index'])])
with netCDF4.Dataset(r'E:\soc\l1r\2023\11\awe_l1r_q20_2023326T0542_00004_v02.nc') as src,netCDF4.Dataset(out/'original_v02_matched_frames.nc') as sub:
    for name,v in src.variables.items():
        sl=tuple(d['frame_index'] if dim=='time' else slice(None) for dim in v.dimensions)
        np.testing.assert_equal(np.asarray(v[sl]),np.asarray(sub[name][:]))
print('Verified exact Epoch matching and every source subset variable.',metrics,flush=True)
# Execute all notebook cells except dependency installation with a short fit.
nb=json.loads((out/'AWE_v02_actual_data.ipynb').read_text());ns={}
import matplotlib
matplotlib.use('Agg')
for i,c in enumerate(nb['cells']):
    if c['cell_type']!='code':continue
    code=''.join(c['source'])
    if '%pip' in code:continue
    code=code.replace("DATA_FOLDER = Path('.')",f'DATA_FOLDER = Path({str(out)!r})').replace('ITERATIONS = 64','ITERATIONS = 2')
    exec(compile(code,f'v02_cell_{i}','exec'),ns)
files=[p for p in out.iterdir() if p.is_file() and p.name!='manifest.json']
(out/'manifest.json').write_text(json.dumps({p.name:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in files},indent=2))
zip_path=s.ROOT/'AWE_v02_notebook_and_epoch_matched_frames.zip'
with zipfile.ZipFile(zip_path,'w',compression=zipfile.ZIP_DEFLATED) as z:
    for p in files+[out/'manifest.json']:z.write(p,arcname=out.name+'/'+p.name)
with zipfile.ZipFile(zip_path) as z:assert z.testzip() is None
print('PASS notebook execution and archive integrity. ZIP MB:',zip_path.stat().st_size/1e6)
