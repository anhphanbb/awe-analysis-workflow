"""Reproduce the v02 notebook near orbit 19 05:23:40 UTC; then run juwavelet.



Uses the notebook's matching, warping and alternating-median fit code unchanged.

"""

import os

for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'

from pathlib import Path

import sys,json,warnings,hashlib

ROOT=Path(__file__).resolve().parent

sys.path.insert(0,str(ROOT/'.packages'))

os.environ['MPLCONFIGDIR']=str(ROOT/'.mplcache')

import numpy as np

import netCDF4

from scipy.interpolate import LinearNDInterpolator

from scipy.ndimage import gaussian_filter,map_coordinates,maximum_filter,distance_transform_edt

from scipy.signal import fftconvolve

import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt

from datetime import datetime,timedelta,timezone



OUT=ROOT/'AWE_v02_orbit19_052340';OUT.mkdir(exist_ok=True)

SOURCE=Path('E:/soc/l1r/2023/11/awe_l1r_q20_2023327T0456_00019_v02.nc')

NOTEBOOK=ROOT/'AWE_v02_actual_data/AWE_v02_actual_data.ipynb'

JUDIR=Path('C:/Users/domin/juwavelet')

ANCHOR=datetime(2023,1,1,tzinfo=timezone.utc)

def utc(e):return (ANCHOR+timedelta(seconds=float(e)-1356566418)).isoformat(timespec='milliseconds').replace('+00:00','Z')

def arr(a):return np.asarray(np.ma.filled(np.ma.asarray(a,dtype=float),np.nan))

def save(name,obj):(OUT/name).write_text(json.dumps(obj,indent=2,allow_nan=False),encoding='utf-8')



def fit():

    if (OUT/'fit_summary.json').exists():return

    target=1356566418+(datetime(2023,11,23,5,23,40,tzinfo=timezone.utc)-ANCHOR).total_seconds()

    with netCDF4.Dataset(SOURCE) as d:

        epochs=arr(d['Epoch'][:]);reference=int(np.argmin(abs(epochs-target)))

        frames=np.arange(reference-25,reference+25);ep=epochs[frames]

        raw=arr(d['Radiance'][frames]);la=arr(d['Latitude'][frames]);lo=arr(d['Longitude'][frames])

    raw[(raw<0)|(raw>65535)]=np.nan

    good=np.isfinite(la[25])&np.isfinite(lo[25]);step=.02

    lat=np.arange(np.floor(la[25][good].min()/step),np.ceil(la[25][good].max()/step)+1)*step

    lon=np.arange(np.floor(lo[25][good].min()/step),np.ceil(lo[25][good].max()/step)+1)*step

    xx,yy=np.meshgrid(lon,lat);cube=[]

    for j,frame in enumerate(frames):

        valid=np.isfinite(la[j])&np.isfinite(lo[j])

        interp=LinearNDInterpolator(np.column_stack((lo[j][valid],la[j][valid])),raw[j][valid])

        cube.append(interp(xx,yy).astype('f4'))

        if j%10==0:print('Projected source frame',frame,flush=True)

    I=np.asarray(cube,dtype=float);times=ep-ep[25];ny,nx=I.shape[1:];yy,xx=np.indices((ny,nx),dtype=float)

    assert np.all(np.diff(ep)>0) and reference==1467 and abs(ep[25]-target)<.6

    data=dict(radiance=I.astype('f4'),latitude=lat,longitude=lon,epoch=ep,frame_index=frames,common_coverage=np.isfinite(I).all(axis=0))

    np.savez_compressed(OUT/'projected_frames.npz',**data)

    nb=json.loads(NOTEBOOK.read_text());OUTPUT=OUT/'actual_results_64';OUTPUT.mkdir(exist_ok=True)

    ns=dict(np=np,gaussian_filter=gaussian_filter,map_coordinates=map_coordinates,maximum_filter=maximum_filter,

            fftconvolve=fftconvolve,warnings=warnings,json=json,I=I,times=times,ny=ny,nx=nx,yy=yy,xx=xx,

            OUTPUT=OUTPUT,data=data,lat=lat,lon=lon,frames=frames,ITERATIONS=64)

    exec(compile(''.join(nb['cells'][6]['source']),str(NOTEBOOK)+'#matching','exec'),ns)

    candidates=ns['motions']()

    va,na,ca=ns['mode_velocity'](candidates,slow=True);vc,nc,cc=ns['mode_velocity'](candidates,slow=False)

    assert np.isfinite(va).all() and np.isfinite(vc).all() and min(na,nc)>=3

    ns.update(va=va,vc=vc)

    print('Fresh velocities:',va,vc,'support',na,nc,flush=True)

    exec(compile(''.join(nb['cells'][9]['source']),str(NOTEBOOK)+'#fit','exec'),ns)

    np.savez_compressed(OUT/'fitted_sequence.npz',slow=ns['slowseq'].astype('f4'),fast=ns['fastseq'].astype('f4'),residual=ns['residual'].astype('f4'))

    summary=dict(orbit=19,requested_utc='2023-11-23T05:23:40.000Z',reference_frame=reference,reference_utc=utc(ep[25]),

        first_frame=int(frames[0]),last_frame=int(frames[-1]),first_utc=utc(ep[0]),last_utc=utc(ep[-1]),

        frames=len(frames),duration_seconds=float(times[-1]-times[0]),epochs=ep.tolist(),utc=[utc(e) for e in ep],

        source=str(SOURCE),source_size=SOURCE.stat().st_size,source_mtime_ns=SOURCE.stat().st_mtime_ns,

        reference_altitude_km=85,grid_shape=list(I.shape[1:]),grid_step_degrees=step,

        slow_velocity_pixels_per_s=va.tolist(),fast_velocity_pixels_per_s=vc.tolist(),

        slow_peak_support=na,fast_peak_support=nc,slow_median_ncc=ca,fast_median_ncc=cc,

        iterations=64,initial_rmse_kR=ns['history'][0]['reconstruction_rmse_kR'],final_rmse_kR=ns['history'][-1]['reconstruction_rmse_kR'],

        evaluation_pixels_per_frame=int(ns['fixed_mask'].sum()),notebook_sha256=hashlib.sha256(NOTEBOOK.read_bytes()).hexdigest(),

        method='Same matching and 64-iteration alternating-median fit as the original v02 notebook. Fresh velocities. All 50 frames fit; no heldout validation. Minimum 10-frame template support; retain prior template where updates lack support; fast median fixed to zero.',

        projection='0.02-degree latitude/longitude grid covering the reference frame. Linear interpolation retains NaN radiance vertices. No temporal interpolation.',

        utc_conversion='AWE processing-code anchor: GPS Epoch 1356566418 = 2023-01-01T00:00:00 UTC; GPS-UTC = 18 s for this acquisition.')

    save('fit_summary.json',summary)



def wavelets(max_wavelength=200.,cluster_threshold=.1,angle_step=5,iterations=64):

    if (OUT/'wavelet_summary.json').exists():
        previous=json.loads((OUT/'wavelet_summary.json').read_text())
        settings=previous.get('settings',{})
        if previous.get('footprint')=='all_finite_slow' and settings.get('max_wavelength_exclusive_km')==max_wavelength and settings.get('cluster_threshold')==cluster_threshold and settings.get('angle_step_degrees')==angle_step and settings.get('separation_iterations',64)==iterations:return

    sys.path.insert(0,str(JUDIR));from juwavelet import transform

    from wave_cluster_limited import identify_cluster2d

    data=np.load(OUT/'projected_frames.npz');fit=np.load(OUT/f'actual_results_{iterations}/checkpoint_{iterations:04d}.npz')

    slow=fit['slow'];obs=data['radiance'][25];support=fit['slow_support']

    valid=np.isfinite(slow)  # Match the displayed slow field; no extra observation/support mask.

    yy,xx=np.indices(slow.shape);A=np.column_stack((np.ones(valid.sum()),xx[valid],yy[valid]))

    plane=np.linalg.lstsq(A,slow[valid],rcond=None)[0]

    detrended=np.where(valid,slow-(plane[0]+plane[1]*xx+plane[2]*yy),0)

    dy=np.deg2rad(.02)*6371;dx=dy*np.cos(np.deg2rad(data['latitude'].mean()))

    distance=distance_transform_edt(np.pad(valid,1),sampling=(dy,dx))[1:-1,1:-1]

    taper=.5-.5*np.cos(np.pi*np.clip(distance/40,0,1));scale=float(np.std(detrended[distance>=40]))

    assert scale>0 and (distance>=40).sum()>100

    proc=detrended*taper/scale;minimum=12.;dj=.125;js=int(np.ceil(np.log2(max_wavelength/minimum)/dj));jt=int(180/angle_step)
    assert angle_step>0 and 180%angle_step==0

    # Independent synthetic test of x/y convention, kilometre scale and axis angle.

    y,x=np.indices((81,81));lam=48.;theta=np.deg2rad(30)

    synth=np.cos(2*np.pi/lam*(dx*x*np.cos(theta)+dy*y*np.sin(theta)))

    with np.errstate(divide='ignore',invalid='ignore'):

        test=transform.decompose2d(synth.T,dx=dx,dy=dy,s0=24,dj=.125,js=17,jt=36,mode='scaled',dtype=np.complex64)

        power=np.mean(abs(test['decomposition'][:,:,30:51,30:51])**2,axis=(2,3));si,ti=np.unravel_index(np.argmax(power),power.shape)

        assert abs(test['period'][si]-lam)<5 and abs(test['theta'][ti]-theta)<np.deg2rad(6)

        del test

        print('juwavelet synthetic wavelength and orientation check PASS; transforming slow field',flush=True)

        bank=transform.decompose2d(proc.T,dx=dx,dy=dy,s0=minimum,dj=dj,js=js,jt=jt,mode='scaled',dtype=np.complex64)

    coeff=bank['decomposition'];spectrum=np.mean(abs(coeff)**2,axis=(2,3))

    print('Identifying strongest four wavelet clusters',flush=True)

    baseline=None
    if cluster_threshold<.1:
        ba,bi,bl=identify_cluster2d(bank,min_amp=.25,thr=.1,max_clusters=4)
        baseline=dict(threshold=.1,clusters=len(ba),selected_coefficients=int((bl>=0).sum()),spatial_pixels=int(np.any(bl>=0,axis=(0,1)).sum()))
        del bl
    amps,indices,labels=identify_cluster2d(bank,min_amp=.25,thr=cluster_threshold,max_clusters=4)
    expansion=dict(baseline=baseline,selected_coefficients=int((labels>=0).sum()),spatial_pixels=int(np.any(labels>=0,axis=(0,1)).sum()))

    maps=[];info=[]

    for k,amp in enumerate(amps):

        weights=np.zeros(coeff.shape[:2]);partial=np.zeros_like(coeff)

        for si in range(js):

            for ti in range(jt):

                mask=labels[si,ti]==k;partial[si,ti][mask]=coeff[si,ti][mask]

                weights[si,ti]=abs(coeff[si,ti][mask]).sum()

        wavelength=float(np.sum(weights*bank['period'][:,None])/weights.sum())

        axis=float(np.degrees(np.angle(np.sum(weights*np.exp(2j*bank['theta'][None,:])))/2)%180)

        rec=transform.reconstruct2d(dict(bank,decomposition=partial)).real.T*scale

        maps.append(np.where(valid,rec,np.nan));si,ti,px,py=indices[k]

        info.append(dict(rank=k+1,peak_normalized_amplitude=float(amp),weighted_wavelength_km=wavelength,

            selected_coefficients=int((labels==k).sum()),spatial_pixels=int(np.any(labels==k,axis=(0,1)).sum()),wave_normal_axis_degrees=axis,peak_wavelength_km=float(bank['period'][si]),peak_axis_degrees=float(np.degrees(bank['theta'][ti])),

            peak_longitude=float(data['longitude'][px]),peak_latitude=float(data['latitude'][py]),

            peak_edge_distance_km=float(distance[py,px]),search_edge=bool(si in (0,js-1))))

        del partial

        print('Cluster',k+1,info[-1],flush=True)

    np.savez_compressed(OUT/'wavelet_results.npz',slow=np.where(valid,slow,np.nan),detrended=np.where(valid,detrended,np.nan),taper=taper,

        valid=valid,clusters=np.asarray(maps),period_km=bank['period'],axis_degrees=np.degrees(bank['theta']),spectrum=spectrum)

    save('wavelet_summary.json',dict(source=f'Fitted slow reference template, frame 1467 at iteration {iterations}',footprint='all_finite_slow',finite_pixels=int(valid.sum()),missing_pixels=int((~valid).sum()),clusters=info,cluster_expansion=expansion,

        settings=dict(separation_iterations=iterations,min_wavelength_km=minimum,max_wavelength_km=float(bank['period'][-1]),scale_step_octaves=dj,angle_step_degrees=angle_step,

        dx_km=float(dx),dy_km=float(dy),taper_km=40,min_amplitude=.25,cluster_threshold=cluster_threshold,max_wavelength_exclusive_km=max_wavelength,wavelength_steps=js,angle_steps=jt,mode='scaled',normalization_std_kR=scale),

        checks=['Synthetic 48 km / 30 degree wave recovered within scale and angle resolution'],

        transform_sha256=hashlib.sha256((JUDIR/'juwavelet/transform.py').read_bytes()).hexdigest(),

        notes='Plane removed; cosine taper over 40 km from missing-data edges; zero outside valid footprint for FFT. Up to four strongest coefficient clusters, not significance-tested independent waves. Wavelength uses local geographic spacing. Axis is wave normal counterclockwise from east modulo 180 degrees; no travel direction or speed inferred.'))



if __name__=='__main__':

    fit();wavelets();print('ANALYSIS COMPLETE',OUT,flush=True)

