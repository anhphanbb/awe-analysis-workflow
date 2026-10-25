"""Run an isolated whole-orbit analysis from L1R or projected L1F, including verification.

python run_orbit.py --orbit 19 --workers 6 --iterations 32
"""
import os,sys,json,argparse
from pathlib import Path

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--orbit',type=int,required=True)
    ap.add_argument('--source',choices=['l1r','l1f'],default='l1r')
    ap.add_argument('--data-root',type=Path,default=Path('E:/soc'))
    ap.add_argument('--workers',type=int,default=6)
    ap.add_argument('--iterations',type=int,default=32)
    ap.add_argument('--prepare-only',action='store_true')
    args=ap.parse_args()
    root=Path(__file__).resolve().parent
    matches=list((args.data_root/args.source).rglob(f'awe_{args.source}_q20_*_{args.orbit:05d}_v02.nc'))
    if len(matches)!=1:raise ValueError(f'Expected one {args.source.upper()} Q20 product, found {matches}')
    source=matches[0];relative=source.parent.relative_to(args.data_root/args.source)
    l3=args.data_root/'l3c'/relative/source.name.replace(args.source,'l3c')
    master=args.data_root/'l1r'/relative/f'Master{args.orbit:05d}.mgf'
    for path in (l3,master):
        if not path.exists():raise FileNotFoundError(path)
    prefix=f'orbit{args.orbit:05d}_85km_{args.source}'
    cache=root/(prefix+'_input');cache.mkdir(exist_ok=True)
    config=cache/'config.json'
    settings=dict(orbit=args.orbit,prefix=prefix,source_label='projected L1F' if args.source=='l1f' else 'L1R',input=str(cache/'projected_l1f.nc') if args.source=='l1f' else str(source),original_source=str(source),l3c=str(l3),master=str(master))
    if config.exists() and json.loads(config.read_text())!=settings:raise ValueError('Configuration changed; use a new output directory')
    config.write_text(json.dumps(settings,indent=2))
    os.environ['AWE_ORBIT_CONFIG']=str(config)
    import orbit_full as f
    adapter_checks=[]
    if args.source=='l1f':
        import l1f_adapter
        import verify_l1f_adapter
        adapter_checks=verify_l1f_adapter.check()
        f.savejson(cache/'adapter_verification.json',dict(status='PASS',checks=adapter_checks))
        l1f_adapter.prepare(source,settings['input'])
    if not (f.PILOT/'inventory.json').exists():f.p.inventory()
    if args.prepare_only:return
    f.fit(args.workers,args.iterations)
    f.assemble()
    f.score_orbit()
    f.score_orbit(all_frames=True)
    import verify_orbit_full as verify
    checks=adapter_checks+verify.numerical()+verify.products()
    f.savejson(f.OUT/'verification.json',dict(status='PASS',checks=checks))
    print('COMPLETE:',f.OUT,flush=True)

if __name__=='__main__':main()
