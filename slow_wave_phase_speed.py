"""Local juwavelet phase speed from L1R minus estimated fast radiance.

python slow_wave_phase_speed.py --frame 955
Uses a single accepted separation window spanning the sequence to avoid seams.
No prior wave angle or wavelength. Actual Epoch intervals, numerical radiance,
and a shifted-template control; no PNG quantization or assumed fixed cadence.
"""
import orbit_full as f
from orbit_full import np,arr,netCDF4,Path,json,plt
import sys,argparse,csv,hashlib,textwrap
from scipy.signal.windows import tukey
from scipy.ndimage import map_coordinates

def wavelet_module(path):
    sys.path.insert(0,str(path))
    from juwavelet import transform
    return transform

def prepare(a):
    y,x=np.indices(a.shape);ok=np.isfinite(a)
    if ok.mean()<.95:raise ValueError(f'Only {ok.mean():.1%} valid crop; choose another point or smaller crop')
    A=np.column_stack([np.ones(ok.sum()),x[ok],y[ok]])
    b=np.linalg.lstsq(A,a[ok],rcond=None)[0]
    residual=np.where(ok,a-(b[0]+b[1]*x+b[2]*y),0.)
    return residual*np.outer(tukey(a.shape[0],.3),tukey(a.shape[1],.3))

def decompose(transform,a,minimum,js,jt=36,dj=.125,theta=None):
    filt={} if theta is None else dict(min_theta=theta-1e-7,max_theta=theta+1e-7)
    with np.errstate(divide='ignore',invalid='ignore'):
        return transform.decompose2d(prepare(a).T,dx=2.,dy=2.,s0=minimum,dj=dj,js=js,jt=jt,aspect=1,mode='scaled',filt=filt,dtype=np.complex64)

def phase_fit(times,coeff):
    # Each fixed point has its own spatial phase; relative complex products
    # align these offsets before combining the seven local measurements.
    reference=len(times)//2
    cross=np.sum(coeff*np.conj(coeff[reference]),axis=1)
    if np.any(np.abs(cross)<=1e-12):raise ValueError('Vanishing wavelet phase amplitude')
    phase=np.unwrap(np.angle(cross));t=times-times[reference]
    slope,intercept=np.polyfit(t,phase,1);residual=phase-(slope*t+intercept)
    rms=float(np.sqrt(np.mean(residual**2)))
    return dict(slope_rad_s=float(slope),intercept_rad=float(intercept),phase=phase,fit=slope*t+intercept,phase_scatter_rad=rms,total_fitted_phase_rad=float(abs(slope)*(t[-1]-t[0])),max_adjacent_phase_rad=float(np.max(np.abs(np.diff(phase)))))

def track(transform,cube,times,wavelength,theta,angle_index,jt=36):
    center=(cube.shape[1]-1)/2;offsets=np.linspace(-wavelength/4,wavelength/4,7)/2
    xs=center+offsets*np.cos(theta);ys=center+offsets*np.sin(theta)
    values=[];field=None
    for j,a in enumerate(cube):
        z=decompose(transform,a,wavelength,1,jt=jt,theta=theta)['decomposition'][0,angle_index].T
        if j==len(cube)//2:field=z.real
        values.append(map_coordinates(z.real,[ys,xs],order=1)+1j*map_coordinates(z.imag,[ys,xs],order=1))
    values=np.asarray(values);return phase_fit(times,values),values,field,np.column_stack([xs,ys])

def synthetic_check(transform):
    y,x=np.indices((81,81));theta=np.deg2rad(30);lam=32.;speed=.08
    times=np.arange(25)*1.1;times[13:]+=1.1 # preserve a real gap
    for sign in [1,-1]:
        cube=np.asarray([np.cos(2*np.pi/lam*(2*x*np.cos(theta)+2*y*np.sin(theta)-sign*speed*t)) for t in times])
        fit,_,_,_=track(transform,cube,times,lam,theta,6)
        recovered=-fit['slope_rad_s']*lam/(2*np.pi)
        assert abs(recovered-sign*speed)<.005,(recovered,sign*speed)
    bank=decompose(transform,cube[12],16,17)
    power=np.mean(np.abs(bank['decomposition'][:,:,38:43,38:43])**2,axis=(2,3));si,ti=np.unravel_index(np.argmax(power),power.shape)
    assert abs(bank['period'][si]-lam)<3 and abs(bank['theta'][ti]-theta)<np.deg2rad(6)
    return 'PASS: automatic wavelength/direction search; positive/negative 80 m/s waves; nonuniform timestamps; phase sign and units'

def build_sequence(args):
    rows=f.rows_read();indices=np.arange(args.frame-args.half_frames,args.frame+args.half_frames+1)
    if indices.min()<0 or indices.max()>=len(rows):raise ValueError('Frame sequence outside source range')
    la,lo=f.p.master();row=rows[args.frame]
    cy=row['y0']+args.pixel_y;cx=row['x0']+args.pixel_x
    radius=args.radius;yy,xx=np.mgrid[cy-radius:cy+radius+1,cx-radius:cx+radius+1]
    if xx.min()<0 or yy.min()<0 or xx.max()>=la.shape[1] or yy.max()>=la.shape[0]:raise ValueError('Crop outside Mastergrid')
    results=json.loads((args.orbit_dir/'window_summary.json').read_text())
    candidates=[r for r in results if r['selected_model'] in ('single_slow','two_motion') and r['start']<=indices.min() and r['end']>=indices.max() and r['bounds'][0]<=xx.min() and r['bounds'][1]>xx.max() and r['bounds'][2]<=yy.min() and r['bounds'][3]>yy.max()]
    if not candidates:raise ValueError('No single accepted window spans this crop and time range; reduce radius or half-frames')
    r=min(candidates,key=lambda r:abs(r['reference_frame']-args.frame));f.OUT=args.orbit_dir;c=f.load_model(r)
    offsets=json.loads((args.orbit_dir/'offset_alignment.json').read_text());offset=next((v for group in offsets['components'] for s,v in zip(group['window_starts'],group['offsets_kR']) if s==r['start']),0.)
    times=np.array([rows[i]['epoch'] for i in indices]);observations=[];slows=[];fasts=[];supports=[]
    with netCDF4.Dataset(f.p.L1) as d:
        for i,e in zip(indices,times):
            row=rows[i];a=arr(d['Radiance'][int(i)])[::-1];a[(a<0)|(a>65535)]=np.nan
            sy=yy-row['y0'];sx=xx-row['x0'];inside=(sy>=0)&(sx>=0)&(sy<300)&(sx<300)
            obs=np.full(xx.shape,np.nan);obs[inside]=a[sy[inside],sx[inside]]
            dt=e-r['reference_epoch'];x=xx-r['bounds'][0];y=yy-r['bounds'][2]
            slow=f.sample(c['slow'],c['slow_velocity'],dt,x,y)-offset
            fast=f.sample(c['fast'],c['fast_velocity'],dt,x,y)+offset if r['selected_model']=='two_motion' else np.zeros(xx.shape)
            sc=f.sample(c['slow_support'].astype(float),c['slow_velocity'],dt,x,y);fc=f.sample(c['fast_support'].astype(float),c['fast_velocity'],dt,x,y) if r['selected_model']=='two_motion' else sc
            ok=(sc>=8)&(fc>=8)&np.isfinite(slow)&np.isfinite(fast)
            observations.append(obs);slows.append(np.where(ok,slow,np.nan));fasts.append(np.where(ok,fast,np.nan));supports.append(ok)
    # Local geographic derivative: columns are the East/North displacement
    # (km) for +1 native x/y cell. Convert wavevector using inverse transpose.
    lat0=la[cy,cx];lon0=lo[cy,cx];J=np.zeros((2,2))
    for k,(dy,dx) in enumerate([(0,1),(1,0)]):
        dlon=((lo[cy+dy,cx+dx]-lo[cy-dy,cx-dx]+180)%360)-180
        J[:,k]=[np.deg2rad(dlon)*6371*np.cos(np.deg2rad(lat0))/2,np.deg2rad(la[cy+dy,cx+dx]-la[cy-dy,cx-dx])*6371/2]
    return indices,times,np.asarray(observations),np.asarray(slows),np.asarray(fasts),r,c,J,dict(latitude=float(lat0),longitude=float(lon0),master_x=int(cx),master_y=int(cy),frame=args.frame,pixel_x=args.pixel_x,pixel_y=args.pixel_y,pixel_convention='zero-based native cutout after reversing disk y, consistent with Mastergrid')

def main(args):
    args.output.mkdir(exist_ok=True);transform=wavelet_module(args.juwavelet_dir)
    checks=synthetic_check(transform);print(checks,flush=True)
    indices,times,obs,template,fast,r,c,J,location=build_sequence(args)
    cleaned=obs-fast;cube=cleaned if args.slow_source=='observation-minus-fast' else template
    # Use a fixed common footprint to avoid time-dependent missing-data edges.
    common=np.isfinite(cube).all(axis=0)&np.isfinite(template).all(axis=0)
    if common.mean()<.95:raise ValueError(f'Common supported crop is only {common.mean():.1%}; reduce --radius or choose another pixel')
    cube=np.where(common,cube,np.nan);template=np.where(common,template,np.nan)
    ref=args.half_frames;n=cube.shape[1];mid=n//2
    if not common[mid-3:mid+4,mid-3:mid+4].all():raise ValueError('Target neighborhood contains missing observations; choose another point')
    dj=.125;js=int(np.floor(np.log2(args.max_wavelength/args.min_wavelength)/dj))+1
    # Select the wave from the denoised spatial slow pattern, without using its
    # imposed velocity to choose temporal phase. Measure phase on cube below.
    bank=decompose(transform,template[ref],args.min_wavelength,js,dj=dj)
    power=np.mean(np.abs(bank['decomposition'][:,:,mid-3:mid+4,mid-3:mid+4])**2,axis=(2,3))
    si,ti=np.unravel_index(np.argmax(power),power.shape);lam=float(bank['period'][si]);theta=float(bank['theta'][ti])
    fit,coeff,field,points=track(transform,cube,times,lam,theta,ti)
    control,_,_,_=track(transform,template,times,lam,theta,ti)
    k_grid=(2*np.pi/lam)*2*np.array([np.cos(theta),np.sin(theta)])
    k_en=np.linalg.solve(J.T,k_grid);knorm=np.linalg.norm(k_en);axis=k_en/knorm
    signed_speed=-fit['slope_rad_s']/knorm*1000;control_speed=-control['slope_rad_s']/knorm*1000
    split_vector=J@c['slow_velocity']*1000;split_projected=float(split_vector@axis)
    bearing=float(np.degrees(np.arctan2(axis[0],axis[1]))%360)
    flags=[]
    if si in (0,js-1):flags.append('Dominant wavelength is at a search boundary; extend the range before interpretation.')
    if fit['total_fitted_phase_rad']<3*fit['phase_scatter_rad']:flags.append('Fitted phase change is small relative to scatter; wavelet speed/direction is poorly constrained.')
    if fit['max_adjacent_phase_rad']>np.pi/2:flags.append('Large adjacent phase changes; possible phase wrapping or changing wave structure.')
    summary=dict(location=location,source_l1r=str(f.p.L1),orbit_results=str(args.orbit_dir),separation_window=[r['start'],r['end']],selected_model=r['selected_model'],frames=indices.tolist(),duration_seconds=float(times[-1]-times[0]),source=args.slow_source,common_valid_fraction=float(common.mean()),wavelength_nominal_km=lam,wavelength_local_geographic_km=float(2*np.pi/knorm),wave_axis_native_degrees=float(np.degrees(theta)),positive_axis_bearing_deg=bearing,wavelet_signed_phase_speed_m_s=float(signed_speed),wavelet_propagation_bearing_deg=bearing if signed_speed>=0 else (bearing+180)%360,layer_split_axis_projected_speed_m_s=split_projected,layer_split_full_speed_m_s=float(np.linalg.norm(split_vector)),layer_split_east_north_m_s=split_vector.tolist(),shifted_slow_template_control_speed_m_s=float(control_speed),phase_change_deg=float(np.degrees(fit['total_fitted_phase_rad'])),phase_scatter_deg=float(np.degrees(fit['phase_scatter_rad'])),flags=flags,checks=checks,method='Automatic local wavelet power peak in the saved slow spatial template over wavelengths and 0..175-degree axes; lock that wave filter and seven fixed points across frames; unwrap relative complex phase; linear fit versus actual Epoch seconds. This differs from the previous PNG/cosine-amplitude workflow.',limitations='Observation-minus-fast retains temporal variability but depends on the layer-splitting estimate and includes its fitting frames. Not independent validation or a pure-airglow retrieval. A short interval and weak phase change can limit speed accuracy. Template-control speed largely reflects imposed translation. Speeds are signed normal-to-wave phase speeds, not necessarily a full flow velocity.',search=dict(min_km=args.min_wavelength,max_km=args.max_wavelength,scale_step_octaves=dj,angle_step_deg=5,crop_radius_cells=args.radius),juwavelet_transform_sha256=hashlib.sha256((args.juwavelet_dir/'juwavelet/transform.py').read_bytes()).hexdigest())
    f.savejson(args.output/'summary.json',summary)
    np.savez_compressed(args.output/'slow_sequence.npz',frame_index=indices,epoch=times,observed=obs,estimated_fast=fast,observation_minus_fast=cleaned,shifted_slow_template=template,common_mask=common,wavelet_coefficients=coeff,phase_rad=fit['phase'],sample_points_xy=points,local_grid_jacobian_km=J)
    with (args.output/'phase_series.csv').open('w',newline='') as file:
        w=csv.writer(file);w.writerow(['frame','epoch','relative_seconds','phase_rad','fitted_phase_rad','template_phase_rad']);w.writerows(zip(indices,times,times-times[ref],fit['phase'],fit['fit'],control['phase']))
    fig,axes=plt.subplots(2,3,figsize=(16,9),layout='constrained');limits=np.nanpercentile(obs[ref],[2,98]);extent=[-args.radius*2,args.radius*2,-args.radius*2,args.radius*2]
    for ax,image,title in zip(axes[0],[obs[ref],cube[ref],field],['Original L1R radiance','Slow estimate: '+args.slow_source,'Selected wavelet component (real coefficient)']):
        im=ax.imshow(image,origin='lower',extent=extent,cmap='gray');ax.plot(0,0,'r+');ax.set(title=title,xlabel='Native x distance (nominal km)',ylabel='Native y distance (nominal km)');fig.colorbar(im,ax=ax)
    im=axes[1,0].imshow(power,origin='lower',aspect='auto',extent=[0,180,0,js],cmap='viridis');axes[1,0].plot(np.degrees(theta)+2.5,si+.5,'rx');axes[1,0].set(xlabel='Wave axis (degrees in native grid)',ylabel='Wavelength scale index',title=f'Search in slow template: {lam:.1f} km, {np.degrees(theta):.0f} degrees');fig.colorbar(im,ax=axes[1,0])
    t=times-times[ref];ax=axes[1,1];ax.plot(t,np.degrees(fit['phase']),'o',label='Original minus fast: phase');ax.plot(t,np.degrees(fit['fit']),label='Linear phase fit');ax.plot(t,np.degrees(control['phase']),'--',label='Shifted-template control');ax.set(xlabel='Seconds from selected frame',ylabel='Relative phase (degrees)',title='Same wave filter and fixed points through time');ax.legend(fontsize=8)
    ax=axes[1,2];ax.axis('off');ax.text(0,.95,f'Wavelet phase speed: {signed_speed:+.1f} m/s\nLayer split, along same axis: {split_projected:+.1f} m/s\nLayer split, full magnitude: {np.linalg.norm(split_vector):.1f} m/s\nShifted-template control: {control_speed:+.1f} m/s\n\nPositive-axis bearing: {bearing:.1f} deg from north\nPhase change: {np.degrees(fit["total_fitted_phase_rad"]):.2f} deg\nPhase scatter: {np.degrees(fit["phase_scatter_rad"]):.2f} deg\n\n'+('\n\n'.join(textwrap.fill(flag,48) for flag in flags) if flags else 'Exploratory comparison;\nno independent validation.'),va='top',fontsize=11,wrap=True)
    fig.suptitle(f'L1R frame {args.frame} | {location["latitude"]:.3f} N, {location["longitude"]:.3f} E | {len(indices)} frames / {times[-1]-times[0]:.1f} s');fig.savefig(args.output/'speed_comparison.png',dpi=140);plt.close(fig)
    (args.output/'README.md').write_text('# Slow-component phase speed\n\n'+json.dumps(summary,indent=2)+'\n\nFiles: slow_sequence.npz contains numerical radiance arrays and exact timestamps; phase_series.csv contains fitted phase progression; speed_comparison.png compares both methods. Pixel coordinates are zero-based in the native Mastergrid-oriented cutout, not raw disk-y coordinates. Positive signed speeds use the axis bearing reported in summary.json.\n',encoding='utf-8')
    print(json.dumps(summary,indent=2),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--frame',type=int,default=955);ap.add_argument('--half-frames',type=int,default=12);ap.add_argument('--pixel-x',type=int,default=150);ap.add_argument('--pixel-y',type=int,default=150);ap.add_argument('--radius',type=int,default=40);ap.add_argument('--min-wavelength',type=float,default=12);ap.add_argument('--max-wavelength',type=float,default=80);ap.add_argument('--slow-source',choices=['observation-minus-fast','shifted-template'],default='observation-minus-fast');ap.add_argument('--orbit-dir',type=Path,default=f.p.ROOT/'orbit85_full_filled');ap.add_argument('--juwavelet-dir',type=Path,default=Path('C:/Users/domin/juwavelet'));ap.add_argument('--output',type=Path,default=f.p.ROOT/'slow_phase_frame0955');args=ap.parse_args()
    if not 0<=args.pixel_x<300 or not 0<=args.pixel_y<300:ap.error('Pixel coordinates must be in 0..299')
    if args.half_frames<2 or args.radius<12 or not 0<args.min_wavelength<args.max_wavelength:ap.error('Invalid sequence length, crop, or wavelength range')
    main(args)
