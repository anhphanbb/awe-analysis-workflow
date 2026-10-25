"""Project native L1F onto a Master grid using inverse bilinear geometry.

Only adjacent camera cells are used. All four radiances must be finite;
no extrapolation, hole filling, temporal interpolation or integer rounding.
"""
import orbit_pilot as p
from orbit_pilot import np, netCDF4, arr, Path, savejson
from scipy.spatial import cKDTree

def resample(lat, lon, radiance, target_lat, target_lon):
    h,w=lat.shape
    valid=np.isfinite(lat)&np.isfinite(lon)
    ids=np.flatnonzero(valid)
    tree=cKDTree(p.xyz(lat[valid],lon[valid]))
    _,nearest=tree.query(p.xyz(target_lat,target_lon))
    yy,xx=np.unravel_index(ids[nearest],lat.shape)
    x=xx.astype(float);y=yy.astype(float)
    ta=target_lat.ravel();tb=target_lon.ravel()
    def corners(a,iy,ix):
        return np.array([a[iy,ix],a[iy,ix+1],a[iy+1,ix],a[iy+1,ix+1]])
    def evaluate(c,u,v):
        return c[0]*(1-u)*(1-v)+c[1]*u*(1-v)+c[2]*(1-u)*v+c[3]*u*v
    for _ in range(8):
        ix=np.clip(np.floor(x),0,w-2).astype(int);iy=np.clip(np.floor(y),0,h-2).astype(int)
        u=x-ix;v=y-iy
        a=corners(lat,iy,ix)-ta
        b=(corners(lon,iy,ix)-tb+180)%360-180
        da=evaluate(a,u,v);db=evaluate(b,u,v)
        ax=(a[1]-a[0])*(1-v)+(a[3]-a[2])*v
        ay=(a[2]-a[0])*(1-u)+(a[3]-a[1])*u
        bx=(b[1]-b[0])*(1-v)+(b[3]-b[2])*v
        by=(b[2]-b[0])*(1-u)+(b[3]-b[1])*u
        det=ax*by-ay*bx
        good=np.isfinite(det)&(np.abs(det)>1e-12)
        dx=np.zeros_like(x);dy=np.zeros_like(y)
        np.divide(da*by-db*ay,det,out=dx,where=good)
        np.divide(db*ax-da*bx,det,out=dy,where=good)
        x=np.clip(x-dx,-1,w);y=np.clip(y-dy,-1,h)
    ix=np.clip(np.floor(x),0,w-2).astype(int);iy=np.clip(np.floor(y),0,h-2).astype(int)
    u=x-ix;v=y-iy
    a=corners(lat,iy,ix)-ta;b=(corners(lon,iy,ix)-tb+180)%360-180
    geometry=(u>=0)&(u<=1)&(v>=0)&(v<=1)&(np.abs(evaluate(a,u,v))<1e-6)&(np.abs(evaluate(b,u,v))<1e-6)
    values=corners(radiance,iy,ix)
    supported=geometry&np.all(np.isfinite(values),axis=0)
    result=np.where(supported,evaluate(values,u,v),np.nan)
    return result.reshape(target_lat.shape).astype('f4'),geometry.reshape(target_lat.shape)

def prepare(source, target):
    source=Path(source);target=Path(target)
    signature=dict(path=str(source),size=source.stat().st_size,mtime_ns=source.stat().st_mtime_ns,
                   master=str(p.MASTER),master_mtime_ns=p.MASTER.stat().st_mtime_ns,adapter_version=1)
    if target.exists():
        with netCDF4.Dataset(target) as d:
            if p.json.loads(d.input_signature)!=signature:raise ValueError('Cached projection belongs to different inputs')
        return
    target.parent.mkdir(exist_ok=True)
    la,lo=p.master();ny,nx=la.shape
    sy,sx=np.mgrid[0:ny:4,0:nx:4];tree=cKDTree(p.xyz(la[::4,::4],lo[::4,::4]))
    temporary=target.with_suffix('.partial.nc')
    with netCDF4.Dataset(source) as src, netCDF4.Dataset(temporary,'w') as dst:
        epochs=arr(src['Epoch'][:]);assert np.all(np.diff(epochs)>0)
        for key in src.ncattrs():dst.setncattr(key,src.getncattr(key))
        dst.Data_level='L1F_PROJECTED';dst.Horizontal_Resolution='2 km nominal Master grid'
        dst.input_signature=p.json.dumps(signature)
        dst.projection_policy=__doc__
        dst.createDimension('time',len(epochs));dst.createDimension('y',300);dst.createDimension('x',300)
        dst.createVariable('Epoch','f8',('time',))[:]=epochs
        for key in ['Radiance','Latitude','Longitude']:
            dst.createVariable(key,'f4',('time','y','x'),zlib=True,complevel=1,fill_value=np.nan,chunksizes=(1,300,300))
        counts=[]
        for i in range(len(epochs)):
            ca=float(src['Center_Latitude'][i]);cb=float(src['Center_Longitude'][i])
            _,k=tree.query(p.xyz(ca,cb)[0]);gy=int(sy.ravel()[k]);gx=int(sx.ravel()[k])
            ya,yb=max(0,gy-5),min(ny,gy+6);xa,xb=max(0,gx-5),min(nx,gx+6)
            _,k=cKDTree(p.xyz(la[ya:yb,xa:xb],lo[ya:yb,xa:xb])).query(p.xyz(ca,cb)[0])
            dy,dx=np.unravel_index(k,(yb-ya,xb-xa));y0=ya+dy-150;x0=xa+dx-150
            assert 0<=y0<=ny-300 and 0<=x0<=nx-300
            ta=la[y0:y0+300,x0:x0+300];tb=lo[y0:y0+300,x0:x0+300]
            r=arr(src['Radiance'][i]);r[(r<0)|(r>65535)]=np.nan
            out,geometry=resample(arr(src['Latitude'][i]),arr(src['Longitude'][i]),r,ta,tb)
            dst['Radiance'][i]=out[::-1]
            dst['Latitude'][i]=np.where(geometry,ta,np.nan)[::-1]
            dst['Longitude'][i]=np.where(geometry,tb,np.nan)[::-1]
            counts.append(int(np.isfinite(out).sum()))
            if i%100==0:print('Projected L1F',i,'/',len(epochs),flush=True)
    temporary.replace(target)
    savejson(target.with_suffix('.json'),dict(input=signature,frames=len(epochs),finite_pixels_per_frame=counts,policy=__doc__))
