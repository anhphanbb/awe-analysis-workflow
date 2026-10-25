"""Apply juwavelet to only the saved concentric-wave numerical crop."""
import build_orbit19_actual as b
np=b.np
parent=b.OUT
crop=np.load(parent/'concentric_waves_crop.npz')
data=np.load(parent/'projected_frames.npz')
fit=np.load(parent/'actual_results_64/checkpoint_0064.npz')
ys=np.isin(data['latitude'],crop['latitude']);xs=np.isin(data['longitude'],crop['longitude'])
out=parent/'cropped_wavelet';out.mkdir(exist_ok=True)
(out/'actual_results_64').mkdir(exist_ok=True)
np.savez_compressed(out/'projected_frames.npz',latitude=crop['latitude'],longitude=crop['longitude'],radiance=data['radiance'][:,ys][:,:,xs])
np.savez_compressed(out/'actual_results_64/checkpoint_0064.npz',slow=crop['slow'],slow_support=fit['slow_support'][np.ix_(ys,xs)])
b.OUT=out
b.wavelets(max_wavelength=80.,cluster_threshold=.05)
summary=b.json.loads((out/'wavelet_summary.json').read_text())
summary.update(source='Concentric-wave crop of fitted slow component, frame 1467 at iteration 64',crop_bounds=dict(longitude=[float(crop['longitude'][0]),float(crop['longitude'][-1])],latitude=[float(crop['latitude'][0]),float(crop['latitude'][-1])]),input_shape=list(crop['slow'].shape),analysis_scope='Crop only; plane removal, edge taper and normalization computed on the crop before juwavelet.')
b.save('wavelet_summary.json',summary)
w=np.load(out/'wavelet_results.npz');np.testing.assert_allclose(w['slow'],crop['slow'],equal_nan=True)
print('PASS: juwavelet input equals the saved crop',crop['slow'].shape,flush=True)
