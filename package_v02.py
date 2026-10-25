from pathlib import Path
import sys,json,csv,shutil
ROOT=Path(__file__).parent
sys.path.insert(0,str(ROOT/'.packages'))
import numpy as np
import netCDF4
from scipy.interpolate import LinearNDInterpolator
out=ROOT/'AWE_v02_actual_data';out.mkdir(exist_ok=True)
source=Path(r'E:\soc\l1r\2023\11\awe_l1r_q20_2023326T0542_00004_v02.nc')
old=np.load(ROOT/'parallax_demo/projected_frames.npz')
with netCDF4.Dataset(source) as src:
    epochs=np.asarray(src['Epoch'][:]);matches=[];missing=[]
    for f,e in zip(old['frame_index'],old['epoch']):
        found=np.flatnonzero(epochs==e)
        if len(found)==1:matches.append((int(f),int(found[0]),float(e)))
        elif len(found)==0:missing.append({'v01_frame_index':int(f),'epoch':float(e)})
        else:raise RuntimeError('Duplicate exact Epoch')
    idx=np.array([m[1] for m in matches]);labels=np.array([m[0] for m in matches])
    lat=np.asarray(src['Latitude'][idx]);lon=np.asarray(src['Longitude'][idx]);rad=np.asarray(src['Radiance'][idx])
    altitude=src.OH_layer_altitude
    with netCDF4.Dataset(out/'original_v02_matched_frames.nc','w') as dst:
        dst.setncatts(src.__dict__);dst.setncattr('subset_note','Exact Epoch matches to v01 frames 930-979; see epoch_mapping.csv and metadata.json. Original global time attributes refer to full file.')
        for name,d in src.dimensions.items():dst.createDimension(name,len(idx) if name=='time' else len(d))
        for name,v in src.variables.items():
            attrs=v.__dict__.copy();fill=attrs.pop('_FillValue',None)
            w=dst.createVariable(name,v.dtype,v.dimensions,zlib=True,**({'fill_value':fill} if fill is not None else {}));w.setncatts(attrs)
            sl=tuple(idx if d=='time' else slice(None) for d in v.dimensions);w[:]=v[sl]
        for name,values in [('source_v02_frame_index',idx),('matching_v01_frame_index',labels)]:dst.createVariable(name,'i4',('time',))[:]=values
x,y=old['longitude'],old['latitude'];xx,yy=np.meshgrid(x,y);cube=[]
for j,i in enumerate(idx):
    good=np.isfinite(lat[j])&np.isfinite(lon[j])
    values=rad[j].copy();values[(values<0)|(values>65535)]=np.nan
    # Keep missing-radiance vertices: interpolation must not bridge v02 masked boxes.
    interp=LinearNDInterpolator(np.column_stack((lon[j][good],lat[j][good])),values[good])
    cube.append(interp(xx,yy).astype('float32'))
    print('Projected',int(i),flush=True)
cube=np.array(cube)
np.savez_compressed(out/'projected_frames.npz',radiance=cube,latitude=y,longitude=x,epoch=np.array([m[2] for m in matches]),frame_index=idx,matching_v01_frame_index=labels,common_coverage=np.isfinite(cube).all(axis=0))
with (out/'epoch_mapping.csv').open('w',newline='') as f:
    w=csv.writer(f);w.writerow(['v01_frame_index','v02_frame_index','epoch']);w.writerows(matches)
metadata={'source':str(source),'reference_altitude':altitude,'selection':'Exact float64 Epoch equality; no nearest-time substitution or inserted frames','matched_frames':len(idx),'missing':missing,'grid':'Same latitude/longitude grid as earlier demonstration, 0.02 degrees; v02 coordinates use 85 km','interpolation':'Linear interpolation with NaN radiance vertices retained to preserve missing-data masks','reference_v01_frame':955,'reference_v02_frame':int(idx[labels==955][0])}
(out/'metadata.json').write_text(json.dumps(metadata,indent=2))
nb=json.loads((ROOT/'AWE_actual_data/AWE_actual_data.ipynb').read_text())
for c in nb['cells']:
    text=''.join(c['source'])
    text=text.replace('50 actual frames','49 actual frames').replace('indices **930–979**','indices **930–979**, excluding the missing epoch documented in metadata.json')
    text=text.replace('awe_l1r_q20_2023326T0547_00004_v01.nc',source.name).replace('original_frames_0930_0979.nc','original_v02_matched_frames.nc')
    text=text.replace('**87 km OH-layer projection**','**85 km OH-layer projection from v02**')
    text=text.replace('times = data[\'epoch\'] - data[\'epoch\'][25]',"reference_index = int(np.flatnonzero(data['matching_v01_frame_index'] == 955)[0])\ntimes = data['epoch'] - data['epoch'][reference_index]")
    text=text.replace("assert I.shape == (50, 309, 317)","assert I.shape == (49, 309, 317)")
    text=text.replace('assert frames.tolist() == list(range(930, 980))','assert np.all(np.diff(times) > 0)')
    text=text.replace('[0, 25, 49]','[0, reference_index, len(I)-1]').replace('all 50','all 49')
    text=text.replace('previous actual-data fit','v02 exact-epoch fit')
    text=text.replace('FRAME_TO_VIEW = 955',f"FRAME_TO_VIEW = {metadata['reference_v02_frame']}")
    text=text.replace('Original zero-based frame index, 930 through 979.','Zero-based v02 source frame index; see epoch_mapping.csv.')
    text=text.replace('`time[0]` here is original frame 930; `time[25]` is original frame 955.','Subset time indices map to v02 and v01 indices in epoch_mapping.csv. The reference time is the exact Epoch of v01 frame 955.')
    text=text.replace('`source_frame_index`','`source_v02_frame_index`')
    c['source']=text.splitlines(keepends=True)
nb['cells'][0]['source']+=['\nVersion-specific notes: the original v02 file uses an 85 km layer and contains radiance masks. One requested epoch is missing; the gap is retained. Frame labels in plots are v02 source indices. The same geographic grid as before is retained for comparison.\n']
(out/'AWE_v02_actual_data.ipynb').write_text(json.dumps(nb,indent=2),encoding='utf-8')
shutil.copy(ROOT/'AWE_actual_data/requirements.txt',out/'requirements.txt')
print(json.dumps(metadata),flush=True)
