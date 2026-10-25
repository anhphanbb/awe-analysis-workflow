"""Analytic energy-ratio tests and provenance checks for whole-orbit q."""
import orbit_full as f
from orbit_full import np,netCDF4,arr,json
import argparse

def numerical():
    shape=(61,61);obs=np.full(shape,2.);slow=np.zeros(shape);support=np.ones(shape,bool)
    for fast,expected in [(2.,1.),(1.,.75),(0.,0.),(-2.,-3.)]:
        raw,q,energy,error,cov,ok=f.local_fast_score(obs,slow,np.full(shape,fast),support)
        np.testing.assert_allclose(raw[30,30],expected,atol=1e-12)
        np.testing.assert_allclose(q[30,30],np.clip(expected,0,1),atol=1e-12)
        assert np.isnan(q[0,0]) # no reflected padding inventing edge support
    weak=f.local_fast_score(obs*.01,slow,obs*.01,support)[1]
    assert not np.isfinite(weak).any()
    obs[30,30]=np.nan
    assert np.isnan(f.local_fast_score(obs,slow,np.ones(shape),support)[1][30,30])
    obs[20:40,20:40]=np.nan
    assert np.isnan(f.local_fast_score(obs,slow,np.ones(shape),support)[1][19,30])
    return ['Analytic perfect, partial, zero and worse predictions; clipping retains negative raw diagnostic','Weak energy, missing centers, missing neighborhoods and image edges remain unavailable']

def products():
    with netCDF4.Dataset(f.OUT/'orbit85_separation.nc') as d:
        q=arr(d['Q_Score'][:]);raw=arr(d['Q_Raw'][:]);energy=arr(d['Q_Residual_Energy'][:]);error=arr(d['Q_Error_Energy'][:]);cov=arr(d['Q_Local_Coverage'][:])
        mask=np.asarray(d['Q_Candidate_Mask'][:]);reason=np.asarray(d['Q_Status'][:]);status=np.asarray(d['Status'][:]);choice=np.asarray(d['Dominant_Window_Start'][:]);epoch=arr(d['Pixel_Epoch'][:]);settings=json.loads(d.q_settings_json)
    ok=np.isfinite(q);assert ok.any();assert np.array_equal(ok,reason==5)
    np.testing.assert_allclose(raw[ok],1-error[ok]/energy[ok],atol=2e-6)
    np.testing.assert_allclose(q[ok],np.clip(raw[ok],0,1),atol=1e-7)
    assert np.all(mask[~ok]==-1) and np.array_equal(mask[ok],(q[ok]>settings['threshold']).astype('i1'))
    assert np.all(status[ok]==3) and np.all(cov[ok]>=settings['min_local_coverage']-1e-6)
    assert np.all(energy[ok]>=settings['min_weighted_residual_energy_kR2']-1e-6)
    results=json.loads((f.OUT/'window_summary.json').read_text());rows=f.rows_read()
    with netCDF4.Dataset(f.p.L1) as source:
        for r in results:
            yy,xx=np.where(ok&(choice==r['start']))
            if not len(yy):continue
            assert r['selected_model']=='two_motion'
            times=epoch[yy,xx];allowed=np.array([rows[i]['epoch'] for i in r['heldout_frames']])
            assert np.isin(times,allowed).all()
            for i in r['heldout_frames']:
                take=times==rows[i]['epoch']
                if not take.any():continue
                a=arr(source['Radiance'][i])[::-1];sy=yy[take]-rows[i]['y0'];sx=xx[take]-rows[i]['x0']
                assert ((sy>=0)&(sy<300)&(sx>=0)&(sx<300)).all()
                assert np.isfinite(a[sy,sx]).all()
    summary=json.loads((f.OUT/'q_summary.json').read_text())
    assert summary['eligible_pixels']==int(ok.sum()) and summary['candidate_pixels']==int((mask==1).sum())
    return ['Saved q recomputes from energies; threshold, support, and unavailable masks consistent','Every scored pixel uses an exact heldout epoch of its dominant two-motion model and a genuine finite source observation']

def display_products():
    with netCDF4.Dataset(f.OUT/'orbit85_separation.nc') as d:
        if 'Q_AllFrames_Score' not in d.variables:return []
        q=arr(d['Q_AllFrames_Score'][:]);held=arr(d['Q_Score'][:]);raw=arr(d['Q_AllFrames_Raw'][:]);energy=arr(d['Q_AllFrames_Residual_Energy'][:]);error=arr(d['Q_AllFrames_Error_Energy'][:]);obs=arr(d['Q_AllFrames_Observation'][:]);mask=np.asarray(d['Q_AllFrames_Candidate_Mask'][:]);reason=np.asarray(d['Q_AllFrames_Status'][:]);threshold=json.loads(d.q_all_frames_settings_json)['threshold']
    ok=np.isfinite(q);old=np.isfinite(held)
    np.testing.assert_allclose(q[old],held[old],atol=1e-7)
    assert ok.sum()>old.sum() and not np.any(reason==2)
    assert np.array_equal(ok,reason==5) and np.isfinite(obs[ok]).all()
    np.testing.assert_allclose(raw[ok],1-error[ok]/energy[ok],atol=2e-6)
    assert np.array_equal(mask[ok],(q[ok]>threshold).astype('i1')) and np.all(mask[~ok]==-1)
    for name in ['q_radiance_swaths.png']+[f'q_radiance_segment_{i}.png' for i in range(1,5)]:assert (f.OUT/name).is_file()
    return ['All-frame diagnostic preserves every heldout score, adds directly evaluated timestamps, and retains observation/threshold provenance','Combined radiance, q, masked observed radiance and four enlarged swath plots exist']

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--fill-small-gaps',action='store_true');ap.add_argument('--numerical-only',action='store_true');ap.add_argument('--start-offset',type=int,default=0);args=ap.parse_args()
    if args.fill_small_gaps:f.OUT=f.p.ROOT/'orbit85_full_filled'
    if args.start_offset:f.OUT=f.OUT.with_name(f.OUT.name+f'_shift{args.start_offset:02d}')
    checks=numerical()
    if not args.numerical_only:checks+=products()+display_products()
    report=dict(status='PASS',checks=checks);f.savejson(f.OUT/'q_verification.json',report);print(json.dumps(report,indent=2))
