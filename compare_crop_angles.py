"""Compare 10 and 20 degree wavelet angle grids against the 5 degree baseline."""
import shutil
import build_orbit19_actual as b
parent=b.OUT
source=parent/'cropped_wavelet'
for step in [10,20]:
    out=parent/f'cropped_wavelet_step{step:02d}';out.mkdir(exist_ok=True)
    (out/'actual_results_64').mkdir(exist_ok=True)
    shutil.copy2(source/'projected_frames.npz',out/'projected_frames.npz')
    shutil.copy2(source/'actual_results_64/checkpoint_0064.npz',out/'actual_results_64/checkpoint_0064.npz')
    b.OUT=out
    b.wavelets(max_wavelength=80.,cluster_threshold=.05,angle_step=step)
    summary=b.json.loads((out/'wavelet_summary.json').read_text())
    original=b.json.loads((source/'wavelet_summary.json').read_text())
    for k in ['source','crop_bounds','input_shape','analysis_scope']:summary[k]=original[k]
    b.save('wavelet_summary.json',summary)
    z=b.np.load(out/'wavelet_results.npz');baseline=b.np.load(source/'wavelet_results.npz')
    b.np.testing.assert_array_equal(z['slow'],baseline['slow'])
    b.np.testing.assert_array_equal(z['taper'],baseline['taper'])
    b.np.testing.assert_array_equal(z['period_km'],baseline['period_km'])
    b.np.testing.assert_allclose(z['axis_degrees'],b.np.arange(0,180,step),atol=1e-10)
    assert (z['period_km']<80).all()
    print('PASS: identical crop, taper and wavelengths; angle step',step,flush=True)
