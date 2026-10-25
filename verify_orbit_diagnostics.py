"""Controls and reproducibility checks for the rejected-window diagnosis."""
from diagnose_orbit_windows import *

def main():
    times=(np.arange(50)-25)*1.1;yy,xx=np.indices((300,340),dtype=float)
    rng=np.random.default_rng(91);base=8+gaussian_filter(rng.normal(size=(300,340)),2)*4
    cube=np.array([full.sample(base,np.array([.04,.02]),t,xx,yy) for t in times])
    s,_=match(cube,times,'smooth_1px');assert s['reason']=='single_slow_mode',s['modes'];np.testing.assert_allclose(s['modes'][0]['velocity'],[.04,.02],atol=.02)
    controls=[]
    for seed in (12,35,81):
        rng=np.random.default_rng(seed);noise=8+rng.normal(scale=.4,size=(50,300,340)).astype('f4');r,_=match(noise,times,'smooth_1px');controls.append(dict(seed=seed,reason=r['reason'],retained_candidates=r['counts']['retained_candidates'],modes=r['modes']))
        assert not r['modes'],r
    d=json.loads((OUT/'matching_summary.json').read_text());assert len(d)==76 and sum(x['rejected'] for x in d)==36
    for r in d:
        c=r['variants']['baseline']['counts'];parts=['reference_missing','texture_below_threshold','no_valid_search_position','correlation_below_threshold','boundary_only_patch','patches_with_matches'];assert sum(c.get(k,0) for k in parts)==c['patches_considered']
        if r['rejected']:
            a=np.load(OUT/f"window_{r['start']:04d}/baseline_candidates.npy");b=np.load(OUT/f"window_{r['start']:04d}/lower_texture_candidates.npy");np.testing.assert_array_equal(a,b)
    for r in json.loads((OUT/'refit_summary.json').read_text()):
        original=json.loads((OUT/f"window_{r['start']:04d}/matching.json").read_text())['variants'][r['variant']]
        assert r['modes']==original['modes'];assert not set(r['training_frames'])&set(r['heldout_frames'])
    comparisons=json.loads((OUT/'heldout_comparison.json').read_text());assert len(comparisons)==33
    for r in comparisons:
        assert r['common_heldout_samples']>10000
        np.testing.assert_allclose(r['improvement_over_original_or_stationary'],1-r['rmse_kR'][r['variant']]/r['rmse_kR']['original'])
    report=dict(status='PASS',checks=['Known slow translation recovered after smoothing without forced second mode','Three independent-noise controls produce no supported mode','Rejection accounting exhaustively partitions considered patches','Lowering texture gate leaves rejected-window candidate arrays identical','Refits use training-only diagnostic velocities and disjoint heldout frames','Common-sample comparison arithmetic and coverage'],noise_controls=controls,limitation='Synthetic controls are limited tests and do not establish physical layer correctness.')
    savejson(OUT/'verification.json',report);print(json.dumps(report,indent=2))

if __name__=='__main__':main()
