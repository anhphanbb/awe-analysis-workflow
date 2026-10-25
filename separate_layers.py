"""Exploratory motion-based separation on the fixed 87 km map grid."""
import os, sys, json, warnings
from pathlib import Path
ROOT=Path(__file__).parent
sys.path.insert(0,str(ROOT/'.packages'))
os.environ['MPLCONFIGDIR']=str(ROOT/'.mplcache')
import numpy as np
from scipy.ndimage import gaussian_filter, map_coordinates, maximum_filter
from scipy.signal import fftconvolve
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image

OUT=ROOT/'layer_separation'
OUT.mkdir(exist_ok=True)
D=np.load(ROOT/'parallax_demo/projected_frames.npz')
I=D['radiance'].astype(float)
lat,lon=D['latitude'],D['longitude']
times=D['epoch']-D['epoch'][25]
frames=D['frame_index']
valid=np.isfinite(I)
ny,nx=I.shape[1:]
yy,xx=np.indices((ny,nx),dtype=float)

def highpass(a,sigma=5):
    ok=np.isfinite(a)
    w=gaussian_filter(ok.astype(float),sigma)
    bg=gaussian_filter(np.nan_to_num(a),sigma)/np.maximum(w,1e-9)
    return np.where(ok & (w>.98),a-bg,np.nan)

def ncc(template, search):
    """Normalized cross-correlation; reject any window containing missing data."""
    ok=np.isfinite(search)
    b=np.nan_to_num(search)
    kernel=np.ones(template.shape)
    n=template.size
    count=fftconvolve(ok.astype(float),kernel,mode='valid')
    sm=fftconvolve(b,kernel,mode='valid')
    ss=fftconvolve(b*b,kernel,mode='valid')
    t=template-template.mean()
    denom=np.sqrt(np.maximum(ss-sm*sm/n,0)*np.sum(t*t))
    c=fftconvolve(b,t[::-1,::-1],mode='valid')/np.maximum(denom,1e-12)
    c[count<n-.1]=np.nan
    return c

def motions():
    H=np.array([highpass(a) for a in I])
    rows=[]
    half=15
    for lag in (8,16,24):
        radius=int(np.ceil(lag*1.1))
        for start in (0,12,25):
            end=start+lag
            if end>=len(I):continue
            dt=times[end]-times[start]
            for cy in range(half+radius,ny-half-radius,18):
                for cx in range(half+radius,nx-half-radius,18):
                    t=H[start,cy-half:cy+half+1,cx-half:cx+half+1]
                    if not np.isfinite(t).all() or t.std()<.08:continue
                    s=H[end,cy-half-radius:cy+half+radius+1,cx-half-radius:cx+half+radius+1]
                    c=ncc(t,s)
                    if not np.isfinite(c).any():continue
                    # Keep distinct local peaks, not only the winning displacement.
                    peak=(np.nan_to_num(c,nan=-2)==maximum_filter(np.nan_to_num(c,nan=-2),size=5)) & (c>.35)
                    py,px=np.where(peak)
                    order=np.argsort(c[py,px])[::-1][:4]
                    for rank,k in enumerate(order):
                        dy,dx=float(py[k]-radius),float(px[k]-radius)
                        if abs(dy)>=radius or abs(dx)>=radius:continue
                        # Subpixel quadratic interpolation about each local maximum.
                        for axis in (0,1):
                            p=int(py[k]); q=int(px[k])
                            z=c[p-1:p+2,q] if axis==0 else c[p,q-1:q+2]
                            if len(z)==3 and np.isfinite(z).all():
                                den=z[0]-2*z[1]+z[2]
                                delta=np.clip(.5*(z[0]-z[2])/den,-.5,.5) if den<0 else 0
                                if axis==0:dy+=delta
                                else:dx+=delta
                        rows.append([cx,cy,start,end,dx/dt,dy/dt,float(c[py[k],px[k]]),rank])
    a=np.array(rows)
    np.savetxt(OUT/'motion_candidates.csv',a,delimiter=',',header='x_pixel,y_pixel,start_index,end_index,vx_pixel_per_s,vy_pixel_per_s,ncc,peak_rank',comments='')
    return a

def mode_velocity(rows, slow=False):
    speed=np.hypot(rows[:,4],rows[:,5])
    choose=(speed<.18) if slow else (speed>.22)
    a=rows[choose]
    bins=np.linspace(-1,1,101)
    weights=np.maximum(a[:,6]-.35,0)**2/(1+a[:,7])
    hist,xe,ye=np.histogram2d(a[:,4],a[:,5],bins=(bins,bins),weights=weights)
    hist=gaussian_filter(hist,1)
    ix,iy=np.unravel_index(np.argmax(hist),hist.shape)
    center=np.array([(xe[ix]+xe[ix+1])/2,(ye[iy]+ye[iy+1])/2])
    close=np.linalg.norm(a[:,4:6]-center,axis=1)<.09
    v=np.average(a[close,4:6],axis=0,weights=weights[close])
    return v, int(close.sum()),float(np.median(a[close,6]))

def warp(a,v,dt):
    """Translate a reference image by (vx,vy)*dt; no wrap or extrapolation."""
    coords=[yy-v[1]*dt,xx-v[0]*dt]
    w=map_coordinates(np.isfinite(a).astype(float),coords,order=1,mode='constant',cval=0)
    b=map_coordinates(np.nan_to_num(a),coords,order=1,mode='constant',cval=0)
    return np.where(w>.999,b/np.maximum(w,1e-12),np.nan)

def stack_layer(residual,v):
    a=np.array([warp(r,v,-t) for r,t in zip(residual,times)])
    count=np.isfinite(a).sum(axis=0)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore',RuntimeWarning)
        template=np.nanmedian(a,axis=0)
    template[count<10]=np.nan
    return template,count

def render(A,C,R,confidence,metrics,rows,va,vc):
    plt.rcParams.update({'font.size':10,'figure.facecolor':'#f5f7fa','axes.facecolor':'#dce0e6'})
    kx=111.32*np.cos(np.deg2rad(lat.mean()))*.02
    ky=111.32*.02
    fig,axes=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
    q=rows[:,7]==0
    axes[0].scatter(rows[q,4]*kx,rows[q,5]*ky,c=rows[q,6],s=12,cmap='viridis',vmin=.35,vmax=1,alpha=.65)
    for v,name,color in [(va,'Slow template','#e54b4b'),(vc,'Fast template','#0059c9')]:
        axes[0].plot(v[0]*kx,v[1]*ky,'x',ms=12,mew=3,color=color,label=name)
    axes[0].set(xlabel='Apparent eastward speed (km/s)',ylabel='Apparent northward speed (km/s)',title='Local best correlation peaks')
    axes[0].legend();axes[0].grid(alpha=.25)
    axes[1].hist(np.hypot(rows[q,4]*kx,rows[q,5]*ky),bins=40,color='#6485ac')
    axes[1].set(xlabel='Apparent speed on 87 km map (km/s)',ylabel='Patch / frame-pair count',title='Measured speed distribution (not physical wind)')
    fig.savefig(OUT/'motion_diagnostics.png',dpi=140);plt.close(fig)
    extent=[lon[0]-.01,lon[-1]+.01,lat[0]-.01,lat[-1]+.01]
    fig,axes=plt.subplots(2,2,figsize=(11,9),layout='constrained')
    lim=np.nanpercentile(np.abs(C),98)
    rlim=np.nanpercentile(np.abs(R),98)
    lo,hi=np.nanpercentile(I,[2,98])
    configs=[(I,'Observed radiance','gray',lo,hi),(A,'Slow component — candidate airglow','gray',lo,hi),(C,'Fast component — signed cloud-associated contrast','RdBu_r',-lim,lim),(R,'Unexplained residual','RdBu_r',-rlim,rlim)]
    ims=[]
    for ax,(cube,title,cmap,vmin,vmax) in zip(axes.flat,configs):
        im=ax.imshow(cube[25],origin='lower',extent=extent,cmap=cmap,vmin=vmin,vmax=vmax,aspect=1/np.cos(np.deg2rad(lat.mean())),interpolation='nearest')
        ims.append(im);ax.set(title=title,xlabel='Longitude (°E)',ylabel='Latitude (°N)');ax.grid(alpha=.15)
        fig.colorbar(im,ax=ax,shrink=.8,label='kRayleigh')
    title=fig.suptitle('')
    for j,f in enumerate(frames):
        for im,(cube,*_) in zip(ims,configs):im.set_data(cube[j])
        title.set_text(f'Motion-based two-layer fit | frame {f} | Δt={times[j]:+.1f} s\nFixed scales • missing support is gray • fast component has zero-median convention')
        fig.savefig(OUT/f'frame_{f:04d}.png',dpi=95)
    plt.close(fig)
    pics=[Image.open(OUT/f'frame_{f:04d}.png').convert('RGB') for f in frames]
    pics[0].save(OUT/'separation.gif',save_all=True,append_images=pics[1:],duration=150,loop=0)
    pics[25].save(OUT/'preview.png')
    fig,axes=plt.subplots(1,2,figsize=(11,5),layout='constrained')
    for ax,a,title in zip(axes,[confidence[25],np.where(confidence[25]>.6,I[25],np.nan)],['Fast-pattern coherence near frame 955','Fast-pattern mask (>0.6 coherence)']):
        im=ax.imshow(a,origin='lower',extent=extent,aspect=1/np.cos(np.deg2rad(lat.mean())),vmin=0 if ax==axes[0] else lo,vmax=1 if ax==axes[0] else hi,cmap='viridis' if ax==axes[0] else 'gray')
        ax.set(title=title,xlabel='Longitude (°E)',ylabel='Latitude (°N)');fig.colorbar(im,ax=ax,shrink=.7)
    fig.savefig(OUT/'confidence.png',dpi=130);plt.close(fig)

def main():
    rows=motions()
    va,na,ca=mode_velocity(rows,True)
    vc,nc,cc=mode_velocity(rows,False)
    print('MOTION',json.dumps({'slow':va.tolist(),'fast':vc.tolist(),'slow_support':na,'fast_support':nc,'slow_ncc':ca,'fast_ncc':cc}),flush=True)
    slow,counta=stack_layer(I,va)
    fast=np.where(counta>=10,0.,np.nan)
    history=[]
    # Alternating robust template fits with two constant translation velocities.
    for it in range(8):
        slowseq=np.array([warp(slow,va,t) for t in times])
        fast,countc=stack_layer(I-slowseq,vc)
        # Fix the additive-offset ambiguity: fast template has median zero.
        offset=np.nanmedian(fast);fast-=offset;slow+=offset
        fastseq=np.array([warp(fast,vc,t) for t in times])
        slow,counta=stack_layer(I-fastseq,va)
        slowseq=np.array([warp(slow,va,t) for t in times])
        r=I-slowseq-fastseq
        history.append(float(np.sqrt(np.nanmean(r*r))))
        print('ITER',it,history[-1],flush=True)
    # Predictive check: fit templates on even frames and assess odd frames.
    # This checks interpolation in time, not independent physical validation.
    train=I.copy();train[1::2]=np.nan
    sa,_=stack_layer(train,va)
    for it in range(5):
        sb=np.array([warp(sa,va,t) for t in times])
        fc,_=stack_layer(train-sb,vc)
        fc-=np.nanmedian(fc)
        fb=np.array([warp(fc,vc,t) for t in times])
        sa,_=stack_layer(train-fb,va)
    pred=np.array([warp(sa,va,t) for t in times])+fb
    baseline,_=stack_layer(train,va)
    basepred=np.array([warp(baseline,va,t) for t in times])
    test=np.isfinite(pred)&np.isfinite(basepred)&valid
    test[::2]=False
    test_rmse=float(np.sqrt(np.mean((I[test]-pred[test])**2)))
    baseline_rmse=float(np.sqrt(np.mean((I[test]-basepred[test])**2)))
    # Local explained-energy score after removing slow background. Not probability.
    confidence=[]
    for obs,a,c in zip(I,slowseq,fastseq):
        ok=np.isfinite(obs)&np.isfinite(a)&np.isfinite(c)
        energy=gaussian_filter(np.where(ok,(obs-a)**2,0),3)
        error=gaussian_filter(np.where(ok,(obs-a-c)**2,0),3)
        coverage=gaussian_filter(ok.astype(float),3)
        score=np.clip(1-error/np.maximum(energy,1e-8),0,1)
        score[(coverage<.95)|(energy<.04)]=np.nan
        confidence.append(score)
    confidence=np.array(confidence)
    k=np.array([111.32*np.cos(np.deg2rad(lat.mean()))*.02,111.32*.02])
    metrics={'slow_velocity_pixels_per_s':va.tolist(),'fast_velocity_pixels_per_s':vc.tolist(),'slow_apparent_velocity_east_north_km_s':(va*k).tolist(),'fast_apparent_velocity_east_north_km_s':(vc*k).tolist(),'slow_peak_support':na,'fast_peak_support':nc,'slow_median_ncc':ca,'fast_median_ncc':cc,'iterations_rmse_kR':history,'heldout_odd_frames_two_layer_rmse_kR':test_rmse,'heldout_odd_frames_slow_only_rmse_kR':baseline_rmse,'heldout_pixels':int(test.sum()),'frame_indices_zero_based':[930,979],'model':'Two constant-velocity translated templates; alternating temporal medians; signed fast contrast with median zero.'}
    (OUT/'metrics.json').write_text(json.dumps(metrics,indent=2))
    np.savez_compressed(OUT/'separated_layers.npz',observed=I.astype('float32'),slow=slowseq.astype('float32'),fast=fastseq.astype('float32'),residual=r.astype('float32'),fast_coherence=confidence.astype('float32'),fast_mask=confidence>.6,latitude=lat,longitude=lon,frame_index=frames,relative_seconds=times,slow_template=slow,fast_template=fast,slow_support=counta,fast_support=countc)
    render(slowseq,fastseq,r,confidence,metrics,rows,va,vc)
    print(json.dumps(metrics),flush=True)

if __name__=='__main__':main()
