import combine_orbit_passes as c
from combine_orbit_passes import np,f,json,netCDF4,arr
import argparse

def numerical():
    qa=np.array([np.nan,np.nan,.2,.2,.8,np.nan,.8,.8,.6]);qb=np.array([np.nan,np.nan,np.nan,.2,np.nan,.8,.8,.2,.6]);valid=np.ones(9,bool);valid[0]=False
    s,p=c.agreement(qa,qb,valid)
    np.testing.assert_array_equal(s,[0,1,2,3,4,5,6,7,3])
    reverse,_=c.agreement(qb,qa,valid);np.testing.assert_array_equal(reverse,[0,1,2,3,5,4,6,7,3])
    for shift in [0,12,24]:
        starts=f.window_starts(1908,shift);coverage=np.zeros(1908,int)
        for start in starts:coverage[start:start+50]+=1
        assert coverage.min()>0 and starts[0]==0 and starts[-1]==1858
    assert set(f.window_starts(1908))&set(f.window_starts(1908,12))=={0,1858}
    return ['All agreement categories, strict threshold boundary and pass-swap symmetry','Shifted layouts cover every frame and share only endpoint starts']

def products():
    out=f.p.ROOT/'orbit85_shift_comparison';a=c.read_pass(f.p.ROOT/'orbit85_full_filled');b=c.read_pass(f.p.ROOT/'orbit85_full_filled_shift12')
    summary=json.loads((out/'summary.json').read_text());valid=np.isfinite(a['L3C_Radiance'])&np.isfinite(a['Pixel_Epoch'])
    with netCDF4.Dataset(out/'orbit85_shift_comparison.nc') as d:
        np.testing.assert_array_equal(arr(d['Pixel_Epoch'][:]),a['Pixel_Epoch'])
        for mode,key,suffix in [('withheld','Q_Score','Heldout'),('all_frames','Q_AllFrames_Score','AllFrames')]:
            s,p=c.agreement(a[key],b[key],valid)
            np.testing.assert_array_equal(d[mode+'_Agreement'][:],s);np.testing.assert_array_equal(d[mode+'_Eligible_Passes'][:],p)
            np.testing.assert_array_equal(arr(d['Q_A_'+suffix][:]),a[key]);np.testing.assert_array_equal(arr(d['Q_B_'+suffix][:]),b[key])
            assert sum(summary['modes'][mode]['categories'].values())==valid.size
            assert summary['modes'][mode]['eligible_either']==int((p>0).sum())
        same=np.asarray(d['Same_Dominant_Window'][:]).astype(bool)
        assert np.all(a['Dominant_Window_Start'][same]==b['Dominant_Window_Start'][same])
        assert set(np.unique(a['Dominant_Window_Start'][same]))<={0,1858}
    for name in ['withheld_whole_orbit.png','all_frames_whole_orbit.png']:assert (out/name).is_file()
    return ['Combined q fields preserved exactly; categories and eligibility recompute from both inputs','Same-window endpoint agreement flagged separately; summary and output plots verified']

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--numerical-only',action='store_true');args=ap.parse_args();checks=numerical()
    if not args.numerical_only:checks+=products()
    report=dict(status='PASS',checks=checks)
    if not args.numerical_only:f.savejson(f.p.ROOT/'orbit85_shift_comparison/verification.json',report)
    print(json.dumps(report,indent=2))
