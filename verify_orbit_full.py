"""Check model selection, offset signs, and completed orbit output integrity."""
import orbit_full as f
from orbit_full import np, json, netCDF4, arr, savejson
from pathlib import Path
import argparse

def numerical():
    from small_gap_fill import fill_small
    y,x=np.indices((20,20));original=(x+2*y).astype(float);a=original.copy()
    a[5:7,5:7]=np.nan;a[0,3]=np.nan;a[10:15,10:15]=np.nan;a[3,10]=np.nan
    geometry=np.ones(a.shape,bool);geometry[3,10]=False
    b,mask,sizes=fill_small(a,geometry)
    assert mask.sum()==4 and sizes==[4]
    np.testing.assert_allclose(b[mask],original[mask])
    np.testing.assert_array_equal(b[np.isfinite(a)],a[np.isfinite(a)])
    assert np.array_equal(np.isfinite(b),np.isfinite(a)|mask)
    slow=dict(velocity=[.04,.02],pair_count=5)
    fast=dict(velocity=[.65,0.],pair_count=5)
    assert f.classify(dict(status='fit',modes=[slow],model_components=1))=='single_slow'
    assert f.classify(dict(status='fit',modes=[fast],model_components=1))=='ambiguous_motion'
    assert f.classify(dict(status='fit',modes=[slow,fast],model_components=2,heldout_improvement_fraction=.02))=='single_slow'
    assert f.classify(dict(status='fit',modes=[slow,fast],model_components=2,heldout_improvement_fraction=.2))=='two_motion'
    assert f.classify(dict(status='unresolved'))=='insufficient_evidence'
    rng=np.random.default_rng(8);texture=rng.normal(size=(100,100))
    r=[dict(start=s,reference_epoch=0.,bounds=[0,100,0,100],selected_model='two_motion') for s in [0,25]]
    models={i:dict(fast=texture+v,slow=np.full((100,100),8-v),fast_velocity=np.zeros(2),fast_support=np.full((100,100),25)) for i,v in enumerate([2.,-3.])}
    old=f.OUT;f.OUT=old/'_checks';f.OUT.mkdir(exist_ok=True)
    try:offsets,components,groups=f.offset_graph(r,models)
    finally:f.OUT=old
    assert len(groups)==1
    np.testing.assert_allclose([offsets[0],offsets[1]],[-2.5,2.5],atol=1e-10)
    np.testing.assert_allclose(models[0]['fast']+offsets[0],models[1]['fast']+offsets[1],atol=1e-10)
    for i in models:
        np.testing.assert_allclose(models[i]['slow']-offsets[i]+models[i]['fast']+offsets[i],8+texture,atol=1e-10)
    # Unconnected windows must keep distinct gauges.
    r[1]['start']=100
    f.OUT=old/'_checks'
    try:_,_,groups=f.offset_graph(r,models)
    finally:f.OUT=old
    assert len(groups)==2
    return ['Single-slow selection, ambiguous single-fast rejection, and predictive-benefit selection','Offset alignment sign and exact reconstruction invariance','Independent gauges for disconnected overlap groups']

def products():
    out=f.OUT;cfg=json.loads((out/'run_config.json').read_text());results=json.loads((out/'window_summary.json').read_text());s=json.loads((out/'summary.json').read_text())
    assert [r['start'] for r in results]==cfg['starts'];assert cfg['starts']==f.window_starts(cfg['source_frames'],cfg.get('start_offset',0))
    coverage=np.zeros(cfg['source_frames'],int)
    for r in results:
        assert r['status']!='error'
        coverage[r['start']:r['end']+1]+=1
        if cfg.get('fill_small_gaps'):
            z=np.load(out/f"window_{r['start']:04d}/interpolation_provenance.npz")
            assert int(z['interpolated'].sum())==r['interpolated_frame_pixels']
            if r['status']=='fit':
                ref=np.load(out/f"window_{r['start']:04d}/reference_maps.npz")
                assert not np.any(ref['evaluation_mask']&z['interpolated'][1::2])
                mask=np.isfinite(ref['observed'])
                np.testing.assert_array_equal(ref['observed'][mask],ref['filled_observed'][mask])
        if r['status']!='fit':continue
        assert not set(r['training_frames'])&set(r['heldout_frames'])
        assert len(r['training_frames'])==len(r['heldout_frames'])==25
        assert len({h['heldout_pixels'] for h in r['history']})==1
    assert (coverage>0).all() and len(coverage)==cfg['source_frames']
    with netCDF4.Dataset(out/'orbit85_separation.nc') as d:
        status=np.asarray(d['Status'][:]);a=arr(d['Slow_Component'][:]);b=arr(d['Fast_Contrast'][:]);e=arr(d['Pixel_Epoch'][:]);ids=np.asarray(d['Dominant_Window_Start'][:]);count=np.asarray(d['Contributing_Windows'][:]);support=arr(d['Mean_Training_Support'][:]);group=np.asarray(d['Offset_Component'][:])
        good=status>=2
        assert status.shape==(551,7801)
        assert np.array_equal(np.isfinite(a),good) and np.array_equal(np.isfinite(b),good)
        assert np.isfinite(e[good]).all() and (support[good]>=8-1e-5).all() and (count[good]>=1).all()
        assert np.all(b[status==2]==0) and (group[status==3]>=0).all()
        assert set(np.unique(ids[good]))<=set(cfg['starts'])
        assert int(good.sum())==sum(s['pixels_by_status'][str(k)] for k in [2,3])
        assert set(np.unique(status))<={0,1,2,3}
    assert s['maximum_l3c_epoch_match_error_seconds']<1e-6
    alignment=json.loads((out/'offset_alignment.json').read_text())
    for c in alignment['components']:assert abs(np.mean(c['offsets_kR']))<1e-6
    return [f"All {cfg['source_frames']} frames included in {len(results)} completed windows",'Training/selection separation and fixed evaluation masks','Map dimensions, support, finite values, status and window provenance','Single-layer modeled zero distinguished from missing estimates','Exact L3C source Epoch matching','Zero-mean gauge within each connected alignment group']

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--numerical-only',action='store_true');ap.add_argument('--fill-small-gaps',action='store_true');ap.add_argument('--start-offset',type=int,default=0);args=ap.parse_args()
    if args.fill_small_gaps:f.OUT=f.p.ROOT/'orbit85_full_filled'
    if args.start_offset:f.OUT=f.OUT.with_name(f.OUT.name+f'_shift{args.start_offset:02d}')
    checks=numerical()+['Small-hole interpolation recovers a planar field; large, edge and invalid-geometry gaps stay missing; observations unchanged']
    if args.fill_small_gaps:checks+=['Per-window interpolation provenance checked; interpolated holdout pixels excluded from evaluation'] if not args.numerical_only else []
    if not args.numerical_only:checks+=products()
    report=dict(status='PASS',checks=checks);savejson(f.OUT/('numerical_checks.json' if args.numerical_only else 'verification.json'),report);print(json.dumps(report,indent=2))
