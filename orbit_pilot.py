"""85 km orbit inventory and overlapping-window pilot on the native Master grid.

Run with Python 3.12 and workspace .packages. Inputs are read-only.
"""
from pathlib import Path
import os, sys, json, csv, argparse, warnings, time
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/'.packages'))
os.environ['MPLCONFIGDIR'] = str(ROOT/'.mplcache')
import numpy as np
import netCDF4
from scipy.spatial import cKDTree
from scipy.ndimage import gaussian_filter, maximum_filter, map_coordinates
from scipy.signal import fftconvolve
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = ROOT/'orbit85_pilot'
BASE = Path('E:/soc')
STEM = 'q20_2023326T0542_00004_v02.nc'
L1 = BASE/'l1r/2023/11'/('awe_l1r_'+STEM)
L3 = BASE/'l3c/2023/11'/('awe_l3c_'+STEM)
MASTER = BASE/'l1r/2023/11/Master00004.mgf'
ORBIT = 4
PREFIX = 'orbit85'
SOURCE_LABEL = 'L1R'
if os.environ.get('AWE_ORBIT_CONFIG'):
    settings=json.loads(Path(os.environ['AWE_ORBIT_CONFIG']).read_text())
    ORBIT=settings['orbit'];PREFIX=settings['prefix'];SOURCE_LABEL=settings['source_label']
    L1=Path(settings['input']);L3=Path(settings['l3c']);MASTER=Path(settings['master'])
    OUT=ROOT/(PREFIX+'_pilot')

def arr(v):
    return np.asarray(np.ma.filled(np.ma.asarray(v, dtype=float), np.nan))

def savejson(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False), encoding='utf-8')

def xyz(lat, lon):
    a,b=np.deg2rad(lat),np.deg2rad(lon)
    return np.column_stack((np.cos(a).ravel()*np.cos(b).ravel(), np.cos(a).ravel()*np.sin(b).ravel(), np.sin(a).ravel()))

def master():
    nx,ny=map(int,np.fromfile(MASTER,dtype='<i4',count=2))
    assert MASTER.stat().st_size == 8+16*nx*ny
    la=np.memmap(MASTER,dtype='<f8',mode='r',offset=8,shape=(nx,ny)).T
    lo=np.memmap(MASTER,dtype='<f8',mode='r',offset=8+8*nx*ny,shape=(nx,ny)).T
    return la,lo

def inventory():
    OUT.mkdir(exist_ok=True)
    la,lo=master(); ny,nx=la.shape
    sy,sx=np.mgrid[0:ny:4,0:nx:4]
    tree=cKDTree(xyz(la[::4,::4],lo[::4,::4]))
    rows=[]; maxerr=0.
    with netCDF4.Dataset(L1) as d:
        epochs=arr(d['Epoch'][:]); assert np.all(np.diff(epochs)>0)
        for i in range(len(epochs)):
            lat=arr(d['Latitude'][i])[::-1]; lon=arr(d['Longitude'][i])[::-1]
            rad=arr(d['Radiance'][i])[::-1]
            valid=np.isfinite(lat)&np.isfinite(lon)
            yy,xx=np.where(valid)
            if not len(yy):
                rows.append(dict(frame=i,epoch=float(epochs[i]),x0=-1,y0=-1,valid_fraction=0.,median_kR=None,p90_kR=None,geometry_error_deg=None)); continue
            k=np.argmin((yy-150)**2+(xx-150)**2); y,x=int(yy[k]),int(xx[k])
            _,n=tree.query(xyz(lat[y,x],lon[y,x])[0]); gy,gx=int(sy.ravel()[n]),int(sx.ravel()[n])
            ya,yb=max(0,gy-5),min(ny,gy+6); xa,xb=max(0,gx-5),min(nx,gx+6)
            delta=(la[ya:yb,xa:xb]-lat[y,x])**2+(((lo[ya:yb,xa:xb]-lon[y,x]+180)%360-180)*np.cos(np.deg2rad(lat[y,x])))**2
            iy,ix=np.unravel_index(np.argmin(delta),delta.shape)
            y0,x0=ya+iy-y,xa+ix-x
            assert 0<=y0<=ny-300 and 0<=x0<=nx-300
            err=max(float(np.max(np.abs(lat[valid]-la[y0:y0+300,x0:x0+300][valid]))),float(np.max(np.abs((lon[valid]-lo[y0:y0+300,x0:x0+300][valid]+180)%360-180))))
            assert err<2e-5, (i,err)
            maxerr=max(maxerr,err)
            good=np.isfinite(rad)&(rad>=0)&(rad<=65535)
            rows.append(dict(frame=i,epoch=float(epochs[i]),x0=int(x0),y0=int(y0),valid_fraction=float(good.mean()),median_kR=float(np.median(rad[good])) if good.any() else None,p90_kR=float(np.percentile(rad[good],90)) if good.any() else None,geometry_error_deg=err))
            if i%200==0: print('Inventory frame',i,flush=True)
        attrs={k:str(getattr(d,k)) for k in ['OH_layer_altitude','Date_Start','Date_End','Horizontal_Resolution']}
    with (OUT/'frame_inventory.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    with netCDF4.Dataset(L3) as d:
        l3la=arr(d['Latitude'][0])[::-1]; l3lo=arr(d['Longitude'][0])[::-1]
        good=np.isfinite(l3la)&np.isfinite(l3lo)
        e=max(float(np.max(np.abs(l3la[good]-la[good]))),float(np.max(np.abs((l3lo[good]-lo[good]+180)%360-180))))
        assert e<2e-5
        rad=arr(d['Radiance'][0])[::-1]; cloud=arr(d['MLCloud'][0])[::-1]; ep=arr(d['Pixel_Epoch'][0])[::-1]
        valid=np.isfinite(rad)&(rad>0)&np.isfinite(ep)
        rad=np.where(valid,rad,np.nan);cloud=np.where(valid,cloud,np.nan)
        np.savez_compressed(OUT/'l3c_reference.npz',radiance=rad.astype('f4'),mlcloud=cloud.astype('f4'),epoch=ep)
        overview(rad,cloud)
    gaps=np.flatnonzero(np.diff(epochs)>1.5*np.median(np.diff(epochs)))
    report=dict(l1r=str(L1),l3c=str(L3),master=str(MASTER),frames=len(rows),duration_seconds=float(epochs[-1]-epochs[0]),metadata=attrs,master_shape_y_x=list(la.shape),master_max_l1r_coordinate_error_deg=maxerr,master_max_l3c_coordinate_error_deg=e,gaps=[dict(after_frame=int(i),seconds=float(epochs[i+1]-epochs[i])) for i in gaps],source_label=SOURCE_LABEL,note='Inventory covers the available orbit product (nighttime acquisition), not an assumed 90-minute continuous sequence. Native Master cells are nominally 2 km. No further spatial reinterpolation is needed after input preparation: all valid input coordinate pairs were checked against exact grid placements.')
    savejson(OUT/'inventory.json',report)
    fig,ax=plt.subplots(3,1,figsize=(14,8),sharex=True,layout='constrained')
    ax[0].plot([r['frame'] for r in rows],[r['valid_fraction'] for r in rows]);ax[0].set(ylabel='Valid fraction',ylim=(0,1),title=f'Orbit {ORBIT} | {len(rows):,} frames | 85 km')
    ax[1].plot([r['median_kR'] for r in rows]);ax[1].set_ylabel('Median radiance (kR)')
    ax[2].plot(np.arange(1,len(epochs)),np.diff(epochs));ax[2].set(ylabel='Frame interval (s)',xlabel='Input frame index (zero-based)')
    fig.savefig(OUT/'inventory.png',dpi=160);plt.close(fig)
    return rows

def overview(rad,cloud,windows=None):
    fig,ax=plt.subplots(2,1,figsize=(18,5.8),sharex=True,layout='constrained')
    for a,data,title,cmap,lim in zip(ax,[rad,cloud],['L3C Q-line radiance (existing coadded swath)','Existing L3C MLCloud score (not motion-derived cloud radiance)'],['gray','viridis'],[np.nanpercentile(rad,[2,98]),[0,100]]):
        im=a.imshow(data,origin='lower',aspect='auto',cmap=cmap,vmin=lim[0],vmax=lim[1],interpolation='nearest')
        a.set(title=title,ylabel='Master cross-track cell');a.set_facecolor('#ddd')
        fig.colorbar(im,ax=a,label='kR' if a is ax[0] else 'MLCloud (%)',pad=.01)
        if windows:
            for label,x in windows: a.axvline(x,color='tomato',lw=.8);a.text(x,520,label,color='tomato',fontsize=8)
    ax[-1].set_xlabel('Master along-track cell (nominal 2 km/cell; vertical axis stretched)')
    fig.suptitle(f'Orbit {ORBIT:05d} | 85 km | {SOURCE_LABEL}',fontsize=15)
    fig.savefig(OUT/'l3c_and_cloud_orbit.png',dpi=180);plt.close(fig)

def highpass(a):
    ok=np.isfinite(a);w=gaussian_filter(ok.astype(float),5)
    return np.where(w>.98,a-gaussian_filter(np.nan_to_num(a),5)/np.maximum(w,1e-9),np.nan)

def ncc(t,s):
    ok=np.isfinite(s);b=np.nan_to_num(s);n=t.size;k=np.ones(t.shape)
    count=fftconvolve(ok.astype(float),k,mode='valid');sm=fftconvolve(b,k,mode='valid');ss=fftconvolve(b*b,k,mode='valid')
    tt=t-t.mean();den=np.sqrt(np.maximum(ss-sm*sm/n,0)*np.sum(tt*tt))
    c=fftconvolve(b,tt[::-1,::-1],mode='valid')/np.maximum(den,1e-12);c[count<n-.1]=np.nan
    return c

def motions(cube,times,diagnostic_folder=None,smooth_sigma=0):
    h={};rows=[]
    for i in range(0,50,2):
        a=cube[i]
        if smooth_sigma:
            ok=np.isfinite(a);w=gaussian_filter(ok.astype(float),smooth_sigma)
            a=np.where(w>.98,gaussian_filter(np.nan_to_num(a),smooth_sigma)/np.maximum(w,1e-9),np.nan)
        h[i]=highpass(a)
    # Only even frames: held-out odd frames do not inform velocities or templates.
    for start in (0,12,24):
        for lag in (8,16,24):
            end=start+lag
            if end>=50:continue
            dt=times[end]-times[start];radius=int(np.ceil(dt*1.1));half=12
            for y in range(half+radius,cube.shape[1]-half-radius,24):
                for x in range(half+radius,cube.shape[2]-half-radius,24):
                    t=h[start][y-half:y+half+1,x-half:x+half+1]
                    if not np.isfinite(t).all() or t.std()<.08:continue
                    c=ncc(t,h[end][y-half-radius:y+half+radius+1,x-half-radius:x+half+radius+1])
                    c0=np.nan_to_num(c,nan=-2);ys,xs=np.where((c0==maximum_filter(c0,5))&(c0>.35))
                    for k in np.argsort(c0[ys,xs])[::-1][:4]:
                        py,px=int(ys[k]),int(xs[k])
                        if min(py,px)<=0 or max(py,px)>=2*radius:continue
                        dy,dx=float(py-radius),float(px-radius)
                        for axis in (0,1):
                            z=c[py-1:py+2,px] if axis==0 else c[py,px-1:px+2];den=z[0]-2*z[1]+z[2]
                            v=float(np.clip(.5*(z[0]-z[2])/den,-.5,.5)) if np.isfinite(z).all() and den<0 else 0.
                            if axis==0:dy+=v
                            else:dx+=v
                        rows.append([dx/dt,dy/dt,float(c[py,px]),start,end,x,y])
    a=np.array(rows)
    if diagnostic_folder is not None:
        np.savetxt(diagnostic_folder/'motion_candidates.csv',a,delimiter=',',header='vx,vy,ncc,start,end,x,y',comments='')
    if len(a)<20:raise ValueError('Too few motion candidates')
    weights=(a[:,2]-.35)**2
    hist,edges,_=np.histogram2d(a[:,0],a[:,1],bins=np.linspace(-1.1,1.1,111),weights=weights)
    hist=gaussian_filter(hist,1)
    if diagnostic_folder is not None:
        fig,ax=plt.subplots(figsize=(6,5),layout='constrained')
        im=ax.imshow(hist.T,origin='lower',extent=[edges[0],edges[-1],edges[0],edges[-1]],aspect='equal')
        ax.set(xlabel='Along-track cells/s',ylabel='Cross-track cells/s',title='Motion evidence (two modes may be unresolved)')
        fig.colorbar(im,ax=ax,label='Weighted candidate support');fig.savefig(diagnostic_folder/'motion.png',dpi=150);plt.close(fig)
    peaks=np.argwhere((hist==maximum_filter(hist,11))&(hist>.15*hist.max()))
    peaks=sorted(peaks,key=lambda p:hist[tuple(p)],reverse=True)
    modes=[]
    for p in peaks:
        center=np.array([edges[p[0]]+.01,edges[p[1]]+.01]);near=np.linalg.norm(a[:,:2]-center,axis=1)<.085
        if near.sum()<15:continue
        v=np.average(a[near,:2],axis=0,weights=weights[near])
        if any(np.linalg.norm(v-m['velocity'])<.12 for m in modes):continue
        modes.append(dict(velocity=v.tolist(),support=int(near.sum()),median_ncc=float(np.median(a[near,2])),pair_count=len(set(map(tuple,a[near,3:5])))))
    if not modes:raise ValueError('No supported motion mode was resolved')
    # Strongest two modes, labeled by speed; labels do not assert physical identity.
    modes=sorted(modes[:2],key=lambda m:np.linalg.norm(m['velocity']))
    if len(modes)>1 and np.linalg.norm(np.subtract(modes[0]['velocity'],modes[1]['velocity']))<.15:modes=modes[:1]
    return modes,a,hist,edges

def run_window(start,rows,max_iterations=64,fill_small_gaps=False,motion_smoothing=0):
    folder=OUT/f'window_{start:04d}';folder.mkdir(exist_ok=True)
    if (folder/'result.json').exists():
        cached=json.loads((folder/'result.json').read_text())
        if cached.get('fill_small_gaps',False)==fill_small_gaps and cached.get('motion_smoothing',0)==motion_smoothing and cached['status']=='fit' and (cached.get('model_components',2)==1 or cached.get('checkpoint_iteration',64)==max_iterations):return cached
    idx=np.arange(start,start+50);mid=start+25;ref=rows[mid]
    la,lo=master();ny,nx=la.shape
    x0=max(0,ref['x0']-20);x1=min(nx,ref['x0']+320);y0=max(0,ref['y0']);y1=min(ny,ref['y0']+300)
    cube=np.full((50,y1-y0,x1-x0),np.nan,dtype='f4')
    geometry=np.zeros(cube.shape,bool)
    with netCDF4.Dataset(L1) as d:
        for j,i in enumerate(idx):
            r=rows[i];ax,bx=max(x0,r['x0']),min(x1,r['x0']+300);ay,by=max(y0,r['y0']),min(y1,r['y0']+300)
            if ax>=bx or ay>=by:continue
            p=arr(d['Radiance'][int(i)])[::-1];p[(p<0)|(p>65535)]=np.nan
            cube[j,ay-y0:by-y0,ax-x0:bx-x0]=p[ay-r['y0']:by-r['y0'],ax-r['x0']:bx-r['x0']]
            if fill_small_gaps:
                g=np.isfinite(arr(d['Latitude'][int(i)])[::-1])&np.isfinite(arr(d['Longitude'][int(i)])[::-1])
                geometry[j,ay-y0:by-y0,ax-x0:bx-x0]=g[ay-r['y0']:by-r['y0'],ax-r['x0']:bx-r['x0']]
    observed=cube.copy();interpolated=np.zeros(cube.shape,bool)
    if fill_small_gaps:
        from small_gap_fill import fill_small
        for j in range(50):cube[j],interpolated[j],_=fill_small(cube[j],geometry[j])
        np.testing.assert_array_equal(cube[np.isfinite(observed)],observed[np.isfinite(observed)])
        assert np.all(geometry[interpolated]) and not np.any(np.isfinite(observed)&interpolated)
        np.savez_compressed(folder/'interpolation_provenance.npz',interpolated=interpolated,frame_indices=idx,bounds=[x0,x1,y0,y1])
    times=np.array([rows[i]['epoch']-ref['epoch'] for i in idx]);train=np.arange(0,50,2);test=np.arange(1,50,2)
    result=dict(start=int(start),end=int(start+49),reference_frame=int(mid),reference_epoch=ref['epoch'],bounds=[x0,x1,y0,y1],status='unresolved',duration_seconds=float(times[-1]-times[0]))
    result.update(fill_small_gaps=fill_small_gaps,motion_smoothing=motion_smoothing,interpolated_frame_pixels=int(interpolated.sum()),interpolated_per_frame=interpolated.sum(axis=(1,2)).tolist())
    try:
        if motion_smoothing:modes,candidates,hist,edges=motions(cube,times,folder,smooth_sigma=motion_smoothing)
        else:modes,candidates,hist,edges=motions(cube,times,folder)
    except ValueError as e:
        result['reason']=str(e);savejson(folder/'result.json',result);print(result,flush=True);return result
    va=np.array(modes[0]['velocity']);vc=np.array(modes[1]['velocity']) if len(modes)>1 else va.copy();result['modes']=modes
    result['model_components']=len(modes)
    np.savetxt(folder/'motion_candidates.csv',candidates,delimiter=',',header='vx,vy,ncc,start,end,x,y',comments='')
    fig,ax=plt.subplots(figsize=(6,5),layout='constrained');im=ax.imshow(hist.T,origin='lower',extent=[edges[0],edges[-1],edges[0],edges[-1]],aspect='equal');ax.scatter([m['velocity'][0] for m in modes],[m['velocity'][1] for m in modes],c=['cyan','red'][:len(modes)]);ax.set(xlabel='Along-track cells/s',ylabel='Cross-track cells/s',title=f'Motion modes | frames {start}–{start+49}');fig.colorbar(im,ax=ax,label='Weighted candidate support');fig.savefig(folder/'motion.png',dpi=150);plt.close(fig)
    yy,xx=np.indices(cube.shape[1:],dtype=float)
    def warp(a,v,t):
        coords=[yy-v[1]*t,xx-v[0]*t];w=map_coordinates(np.isfinite(a).astype('f4'),coords,order=1,mode='constant',cval=0)
        b=map_coordinates(np.nan_to_num(a),coords,order=1,mode='constant',cval=0)
        return np.where(w>.999,b/np.maximum(w,1e-9),np.nan)
    def predict(a,v,ids):return np.array([warp(a,v,times[j]) for j in ids])
    genuine_cache={}
    def stack(res,v):
        aligned=np.array([warp(r,v,-times[j]) for r,j in zip(res,train)]);count=np.isfinite(aligned).sum(0)
        if fill_small_gaps:
            key=tuple(v)
            if key not in genuine_cache:genuine_cache[key]=np.array([np.isfinite(warp(observed[j],v,-times[j])) for j in train])
            genuine=genuine_cache[key]
            count=(np.isfinite(aligned)&genuine).sum(0)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',RuntimeWarning);a=np.nanmedian(aligned,axis=0)
        a[count<8]=np.nan
        return a,count
    slow,sc=stack(cube[train],va);fast=np.where(np.isfinite(slow),0.,np.nan).astype('f4');slowseq=predict(slow,va,train)
    baseline=predict(slow,va,test);history=[];fixed=None;masktrain=None
    np.savez_compressed(folder/'checkpoint_000.npz',slow=slow,fast=fast,slow_support=sc,fast_support=sc,slow_velocity=va,fast_velocity=vc,bounds=np.array([x0,x1,y0,y1]),epoch=ref['epoch'])
    if len(modes)==1:
        residual=observed[test]-baseline;fixed=np.isfinite(residual)
        if fixed.sum()<10000:raise RuntimeError('Too little held-out support for single layer')
        tr=observed[train]-slowseq;rmse=float(np.sqrt(np.mean(residual[fixed]**2)))
        history=[dict(iteration=1,heldout_rmse_kR=rmse,training_rmse_kR=float(np.sqrt(np.nanmean(tr**2))),heldout_pixels=int(fixed.sum()))]
        quality='single_slow_motion' if np.linalg.norm(va)<.18 and modes[0]['pair_count']>=3 else 'single_motion_ambiguous'
        result.update(status='fit',quality=quality,history=history,slow_only_heldout_rmse_kR=rmse,heldout_improvement_fraction=0.,training_frames=idx[train].tolist(),heldout_frames=idx[test].tolist(),minimum_template_support=8,checkpoint_iteration=1,interpretation='One supported motion mode; no fast component estimated. This is not proof of cloud absence.')
        np.savez_compressed(folder/'checkpoint_001.npz',slow=slow,fast=fast,slow_support=sc,fast_support=sc,slow_velocity=va,fast_velocity=vc,bounds=np.array([x0,x1,y0,y1]),epoch=ref['epoch'])
        np.savez_compressed(folder/'reference_maps.npz',observed=observed[25],filled_observed=cube[25],interpolated=interpolated[25],slow=slow,fast=fast,residual=observed[25]-slow,latitude=la[y0:y1,x0:x1],longitude=lo[y0:y1,x0:x1],evaluation_mask=fixed,frame_indices=idx,relative_seconds=times)
        savejson(folder/'result.json',result);print('Single-layer window',start,'velocity',va,'heldout',rmse,flush=True);return result
    for it in range(1,max_iterations+1):
        c,fc=stack(cube[train]-slowseq,vc);fast=np.where(np.isfinite(c),c,fast)
        off=np.nanmedian(fast);fast-=off;slow+=off
        fastseq=predict(fast,vc,train);c,sc=stack(cube[train]-fastseq,va);slow=np.where(np.isfinite(c),c,slow);slowseq=predict(slow,va,train)
        if it in (1,2,4,8,16,32,64) or it==max_iterations:
            pred=predict(slow,va,test)+predict(fast,vc,test);res=observed[test]-pred;tr=observed[train]-slowseq-fastseq
            if fixed is None:
                fixed=np.isfinite(res)&np.isfinite(baseline);masktrain=np.isfinite(tr)
                if fixed.sum()<10000:raise RuntimeError('Too little held-out support')
                base_rmse=float(np.sqrt(np.mean((observed[test]-baseline)[fixed]**2)))
            assert np.isfinite(res[fixed]).all() and np.isfinite(tr[masktrain]).all()
            history.append(dict(iteration=it,heldout_rmse_kR=float(np.sqrt(np.mean(res[fixed]**2))),training_rmse_kR=float(np.sqrt(np.mean(tr[masktrain]**2))),heldout_pixels=int(fixed.sum())))
            print('Window',start,'iteration',it,'heldout',history[-1]['heldout_rmse_kR'],flush=True)
            if it in (16,32,64) or it==max_iterations:np.savez_compressed(folder/f'checkpoint_{it:03d}.npz',slow=slow,fast=fast,slow_support=sc,fast_support=fc,slow_velocity=va,fast_velocity=vc,bounds=np.array([x0,x1,y0,y1]),epoch=ref['epoch'])
    result.update(status='fit',history=history,slow_only_heldout_rmse_kR=base_rmse,heldout_improvement_fraction=1-history[-1]['heldout_rmse_kR']/base_rmse,training_frames=idx[train].tolist(),heldout_frames=idx[test].tolist(),minimum_template_support=8)
    result['quality']='pilot_supported' if result['heldout_improvement_fraction']>.05 and all(m['pair_count']>=3 for m in modes) else 'ambiguous'
    result['checkpoint_iteration']=max_iterations
    np.savez_compressed(folder/'reference_maps.npz',observed=observed[25],filled_observed=cube[25],interpolated=interpolated[25],slow=slow,fast=fast,residual=observed[25]-slow-fast,latitude=la[y0:y1,x0:x1],longitude=lo[y0:y1,x0:x1],evaluation_mask=fixed,frame_indices=idx,relative_seconds=times)
    savejson(folder/'result.json',result)
    return result

def pilot():
    rows=list(csv.DictReader((OUT/'frame_inventory.csv').open()))
    for r in rows:
        for k in r:r[k]=None if r[k]=='' else (int(r[k]) if k in ('frame','x0','y0') else float(r[k]))
    # Three orbit sections, each with a 25-frame-offset overlap; includes the previous test area.
    starts=[200,225,1154,1179,1650,1675]
    results=[]
    for start in starts:
        print('Starting pilot window',start,flush=True);results.append(run_window(start,rows))
    savejson(OUT/'pilot_summary.json',results)
    z=np.load(OUT/'l3c_reference.npz');overview(z['radiance'],z['mlcloud'],[(str(r['start']),sum(r['bounds'][:2])/2) for r in results])
    return results

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['inventory','pilot','all']);args=p.parse_args()
    if args.mode in ('inventory','all'):inventory()
    if args.mode in ('pilot','all'):pilot()
