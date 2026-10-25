"""Numerical and product-integrity checks for the orbit pilot."""
from plot_orbit_pilot import *

def main():
    rng=np.random.default_rng(42);a=gaussian_filter(rng.normal(size=(80,90)),1)
    moved=sample(a,np.array([.5,-.25]),12.)
    np.testing.assert_allclose(moved[:-3,6:],a[3:,:-6],atol=1e-12)
    assert np.isnan(moved[:,:6]).all() and np.isnan(moved[-3:]).all()
    c=ncc(a[25:46,30:51],moved[15:56,20:61]);py,px=np.unravel_index(np.nanargmax(c),c.shape)
    assert (px-10,py-10)==(6,-3)
    base=8+gaussian_filter(rng.normal(size=(300,340)),2)*4
    times=(np.arange(50)-25)*1.1
    cube=np.array([sample(base,np.array([.04,.02]),t) for t in times])
    modes,_,_,_=motions(cube,times)
    assert len(modes)==1, modes
    np.testing.assert_allclose(modes[0]['velocity'],[.04,.02],atol=.02)
    results=json.loads((OUT/'pilot_summary.json').read_text());assert len(results)==6
    verified=0
    for r in results:
        if r['status']!='fit':continue
        assert not set(r['training_frames'])&set(r['heldout_frames'])
        assert len(r['training_frames'])==len(r['heldout_frames'])==25
        folder=OUT/f"window_{r['start']:04d}";z=np.load(folder/'reference_maps.npz');c=np.load(folder/f"checkpoint_{r.get('checkpoint_iteration',64):03d}.npz")
        ok=np.isfinite(z['residual']);np.testing.assert_allclose((z['observed']-z['slow']-z['fast'])[ok],z['residual'][ok],atol=2e-6)
        for layer in ('slow','fast'):np.testing.assert_array_equal(z[layer],c[layer])
        assert len(set(h['heldout_pixels'] for h in r['history']))==1
        matched=np.load(folder/'l3c_time_matched.npz');keep=np.isfinite(matched['fast'])
        assert keep.sum()==r['l3c_matched_pixels']
        t=matched['epoch'][keep]-r['reference_epoch'];assert np.all(t>=z['relative_seconds'][12]) and np.all(t<=z['relative_seconds'][37])
        verified+=1
    with netCDF4.Dataset(OUT/'pilot_on_mastergrid.nc') as d:
        f=arr(d['Pilot_Fast_Contrast'][:]);s=arr(d['Pilot_Slow'][:]);e=arr(d['Pixel_Epoch'][:]);ids=np.asarray(d['Pilot_Window_Start'][:].filled(-1))
        assert f.shape==(551,7801)
        assert np.array_equal(np.isfinite(f),ids>=0)
        assert np.isfinite(s[np.isfinite(f)]).all() and np.isfinite(e[np.isfinite(f)]).all()
        assert set(np.unique(ids[ids>=0]))<=set(r['start'] for r in results)
    report=dict(status='PASS',fitted_windows=verified,checks=['Known translation and no edge wrapping','NCC displacement direction','Single-motion synthetic scene yields one mode without a forced fast layer','Disjoint training and held-out frames','Reconstruction identity','Saved templates match plotted reference maps','Fixed held-out evaluation masks','L3C pixel times within central-window limits','Master output shape and per-pixel provenance'])
    savejson(OUT/'verification.json',report);print(json.dumps(report,indent=2))

if __name__=='__main__':main()
