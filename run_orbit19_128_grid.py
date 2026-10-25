"""Continue the same separation through iteration 128; crop-only threshold/angle grid."""
import shutil
import build_orbit19_actual as b
np=b.np
parent=b.OUT

def extend_fit():
    summary=b.json.loads((parent/'fit_summary.json').read_text())
    if summary['iterations']==128:return
    for name in ['fit_summary.json','fitted_sequence.npz','orbit19_052340.html']:
        backup=parent/('iteration64_'+name)
        if not backup.exists():shutil.copy2(parent/name,backup)
    data=np.load(parent/'projected_frames.npz');checkpoint=np.load(parent/'actual_results_64/checkpoint_0064.npz')
    I=data['radiance'].astype(float);times=data['epoch']-data['epoch'][25];ny,nx=I.shape[1:];yy,xx=np.indices((ny,nx),dtype=float)
    output=parent/'actual_results_128';output.mkdir(exist_ok=True)
    history=b.json.loads((parent/'actual_results_64/metrics.json').read_text())
    ns=dict(np=np,gaussian_filter=b.gaussian_filter,map_coordinates=b.map_coordinates,maximum_filter=b.maximum_filter,
        fftconvolve=b.fftconvolve,warnings=b.warnings,json=b.json,I=I,times=times,ny=ny,nx=nx,yy=yy,xx=xx,
        OUTPUT=output,data=data,lat=data['latitude'],lon=data['longitude'],frames=data['frame_index'],ITERATIONS=128,
        va=checkpoint['slow_velocity'],vc=checkpoint['fast_velocity'],checkpoint=checkpoint,old_history=history)
    nb=b.json.loads(b.NOTEBOOK.read_text())
    exec(''.join(nb['cells'][6]['source']),ns)
    code=''.join(nb['cells'][9]['source'])
    code=code.replace('slow, slow_support = stack_layer(I, va)',"slow=checkpoint['slow'].copy();slow_support=checkpoint['slow_support'].copy()")
    code=code.replace('fast = np.zeros_like(slow)',"fast=checkpoint['fast'].copy()")
    code=code.replace('history = []','history = list(old_history)').replace('fixed_mask = None',"fixed_mask = checkpoint['evaluation_mask'].copy()")
    code=code.replace('range(1, ITERATIONS+1)','range(65, ITERATIONS+1)')
    exec(compile(code,'notebook_fit_continued_to_128','exec'),ns)
    assert len(ns['history'])==128
    np.testing.assert_array_equal(ns['fixed_mask'],checkpoint['evaluation_mask'])
    np.savez_compressed(parent/'fitted_sequence.npz',slow=ns['slowseq'].astype('f4'),fast=ns['fastseq'].astype('f4'),residual=ns['residual'].astype('f4'))
    summary.update(iterations=128,final_rmse_kR=ns['history'][-1]['reconstruction_rmse_kR'],rmse_at_64_kR=history[-1]['reconstruction_rmse_kR'],method=summary['method'].replace('64-iteration','128-iteration')+' Continued exactly from saved iteration 64 with the same velocities and evaluation mask.')
    b.save('fit_summary.json',summary)

def grid():
    data=np.load(parent/'projected_frames.npz');fit=np.load(parent/'actual_results_128/checkpoint_0128.npz')
    ys=(data['latitude']>=-4.8)&(data['latitude']<=-1.8);xs=(data['longitude']>=-59.4)&(data['longitude']<=-56.4)
    crop=fit['slow'][np.ix_(ys,xs)]
    np.savez_compressed(parent/'concentric_waves_crop.npz',slow=crop,latitude=data['latitude'][ys],longitude=data['longitude'][xs],frame=1467)
    results=[]
    for step in [5,10,20]:
        for threshold in [.02,.05]:
            out=parent/'grid128'/f'step{step:02d}_thr{int(threshold*1000):03d}'
            out.mkdir(parents=True,exist_ok=True);(out/'actual_results_128').mkdir(exist_ok=True)
            np.savez_compressed(out/'projected_frames.npz',latitude=data['latitude'][ys],longitude=data['longitude'][xs],radiance=data['radiance'][:,ys][:,:,xs])
            np.savez_compressed(out/'actual_results_128/checkpoint_0128.npz',slow=crop,slow_support=fit['slow_support'][np.ix_(ys,xs)])
            b.OUT=out;b.wavelets(max_wavelength=80.,cluster_threshold=threshold,angle_step=step,iterations=128)
            summary=b.json.loads((out/'wavelet_summary.json').read_text())
            summary.update(source='Concentric-wave crop of fitted slow component, frame 1467 at iteration 128',input_shape=list(crop.shape),analysis_scope='Crop only; plane, taper and normalization recomputed within crop.',crop_bounds=dict(longitude=[-59.4,-56.4],latitude=[-4.8,-1.8]))
            b.save('wavelet_summary.json',summary)
            z=np.load(out/'wavelet_results.npz');np.testing.assert_allclose(z['slow'],crop,equal_nan=True)
            assert np.all(z['period_km']<80)
            results.append(dict(angle_step=step,threshold=threshold,folder=str(out),clusters=summary['clusters'],spatial_pixels=summary['cluster_expansion']['spatial_pixels']))
            print('GRID COMPLETE',step,threshold,flush=True)
    b.OUT=parent;b.save('grid128_summary.json',dict(iterations=128,runs=results,selection='5 degree / 0.02 is the requested default, not an optimized or independently validated winner.'))

if __name__=='__main__':extend_fit();grid()
