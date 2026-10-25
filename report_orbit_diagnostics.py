"""Compare diagnostic refits on identical held-out samples and local overlaps."""
from diagnose_orbit_windows import *

def load_fit(r,root):
    it=r.get('checkpoint_iteration',32)
    if r.get('selected_model')=='single_slow' and r.get('model_components',2)>1:it=0
    with np.load(root/f"window_{r['start']:04d}/checkpoint_{it:03d}.npz") as z:return {k:z[k] for k in z.files}

def prediction(c,times,ids):
    yy,xx=np.indices(c['slow'].shape,dtype=float);pred=[]
    for j in ids:
        a=full.sample(c['slow'],c['slow_velocity'],times[j],xx,yy);b=full.sample(c['fast'],c['fast_velocity'],times[j],xx,yy)
        sa=full.sample(c['slow_support'].astype(float),c['slow_velocity'],times[j],xx,yy);sb=full.sample(c['fast_support'].astype(float),c['fast_velocity'],times[j],xx,yy)
        pred.append(np.where((sa>=8)&(sb>=8),a+b,np.nan))
    return np.array(pred)

def compare_pair(a,ca,b,cb,layer):
    x0=max(a['bounds'][0],b['bounds'][0]);x1=min(a['bounds'][1],b['bounds'][1]);y0=max(a['bounds'][2],b['bounds'][2]);y1=min(a['bounds'][3],b['bounds'][3])
    if x1<=x0 or y1<=y0:return None
    yy,xx=np.mgrid[y0:y1,x0:x1];epoch=(a['reference_epoch']+b['reference_epoch'])/2
    aa=full.sample(ca[layer],ca[layer+'_velocity'],epoch-a['reference_epoch'],xx-a['bounds'][0],yy-a['bounds'][2]);bb=full.sample(cb[layer],cb[layer+'_velocity'],epoch-b['reference_epoch'],xx-b['bounds'][0],yy-b['bounds'][2])
    sa=full.sample(ca[layer+'_support'].astype(float),ca[layer+'_velocity'],epoch-a['reference_epoch'],xx-a['bounds'][0],yy-a['bounds'][2]);sb=full.sample(cb[layer+'_support'].astype(float),cb[layer+'_velocity'],epoch-b['reference_epoch'],xx-b['bounds'][0],yy-b['bounds'][2]);good=np.isfinite(aa)&np.isfinite(bb)&(sa>=8)&(sb>=8)
    if good.sum()<1000 or min(np.std(aa[good]),np.std(bb[good]))<.03:return None
    delta=aa[good]-bb[good];offset=float(np.median(delta));corr=float(np.corrcoef(aa[good],bb[good])[0,1]);return dict(neighbor=b['start'],layer=layer,pixels=int(good.sum()),correlation=corr,offset_removed_rmse_kR=float(np.sqrt(np.mean((delta-offset)**2))),velocity_difference_cells_s=float(np.linalg.norm(ca[layer+'_velocity']-cb[layer+'_velocity'])))

def evaluate():
    orig=json.loads((full.OUT/'window_summary.json').read_text());refs=json.loads((OUT/'refit_summary.json').read_text());diags=json.loads((OUT/'matching_summary.json').read_text());rows=full.rows_read();origmap={r['start']:r for r in orig};records=[]
    accepted_orig={r['start']:r for r in orig if r['selected_model'] in ('single_slow','two_motion')}
    allmodels={('original',s):load_fit(r,full.OUT) for s,r in accepted_orig.items()}
    for r in refs:allmodels[(r['variant'],r['start'])]=load_fit(r,OUT/f"fits_{r['variant']}")
    for d in diags:
        if not d['rejected']:continue
        start=d['start'];fits=[r for r in refs if r['start']==start]
        if not fits:continue
        cube,times,_=cube_for(start,rows);test=np.arange(1,50,2);train=np.arange(0,50,2)
        with warnings.catch_warnings():warnings.simplefilter('ignore',RuntimeWarning);stationary=np.nanmedian(cube[train],axis=0)
        stationary[np.isfinite(cube[train]).sum(0)<8]=np.nan
        pred={'stationary':np.broadcast_to(stationary,cube[test].shape)}
        old=origmap[start]
        if old['status']=='fit':pred['original']=prediction(load_fit(old,full.OUT),times,test)
        else:pred['original']=pred['stationary']
        for r in fits:pred[r['variant']]=prediction(allmodels[(r['variant'],start)],times,test)
        mask=np.isfinite(cube[test])
        for v in pred.values():mask&=np.isfinite(v)
        assert mask.sum()>10000
        scores={k:float(np.sqrt(np.mean((cube[test]-v)[mask]**2))) for k,v in pred.items()}
        for r in fits:
            v=r['variant'];c=allmodels[(v,start)];neighbors=[]
            for ns in sorted(set(accepted_orig)|{q['start'] for q in refs if q['variant']==v and q['selected_model'] in ('single_slow','two_motion')}):
                if not 0<abs(ns-start)<=50:continue
                nr=next((q for q in refs if q['variant']==v and q['start']==ns and q['selected_model'] in ('single_slow','two_motion')),None)
                if nr is not None:nc=allmodels[(v,ns)];source='diagnostic'
                else:nr=accepted_orig[ns];nc=allmodels[('original',ns)];source='original'
                layers=['slow']+(['fast'] if r['selected_model']=='two_motion' and nr['selected_model']=='two_motion' else [])
                for layer in layers:
                    match=compare_pair(r,c,nr,nc,layer)
                    if match:match['neighbor_source']=source;neighbors.append(match)
            target='fast' if r['selected_model']=='two_motion' else 'slow'
            target_checks=[n for n in neighbors if n['layer']==target]
            consistent=any(n['correlation']>=.5 and n['velocity_difference_cells_s']<=.15 for n in target_checks)
            gain=1-scores[v]/scores['original'];stationary_gain=1-scores[v]/scores['stationary']
            outcome='promising' if r['selected_model'] in ('single_slow','two_motion') and gain>.05 and consistent else 'not_confirmed'
            if r['selected_model']=='single_slow' and consistent and 0<=gain<=.05:outcome='consistent_single_slow'
            rec=dict(start=start,end=start+49,original_class=old['selected_model'],variant=v,selected_model=r['selected_model'],comparison_baseline='original_single_motion_fit' if old['status']=='fit' else 'stationary_even_frame_median',common_heldout_samples=int(mask.sum()),common_fraction_of_observed_holdout=float(mask.sum()/np.isfinite(cube[test]).sum()),rmse_kR=scores,improvement_over_original_or_stationary=gain,improvement_over_stationary=stationary_gain,overlap_checks=neighbors,has_consistent_neighbor=consistent,diagnostic_outcome=outcome)
            records.append(rec)
        print('Evaluated',start,flush=True)
    savejson(OUT/'heldout_comparison.json',records)
    summarize(diags,orig,refs,records)

def summarize(diags,original,refs,records):
    rejected=[d for d in diags if d['rejected']];baseline=Counter();accepted=Counter()
    for d in diags:(baseline if d['rejected'] else accepted).update(d['variants']['baseline']['counts'])
    variants={v:dict(Counter(d['variants'][v]['reason'] for d in rejected)) for v in VARIANTS}
    outcomes=dict(Counter(r['diagnostic_outcome'] for r in records));gains=[100*r['improvement_over_original_or_stationary'] for r in records];promising=[r for r in records if r['diagnostic_outcome']=='promising']
    summary=dict(original_rejected_windows=len(rejected),original_accepted_windows=len(diags)-len(rejected),rejected_baseline_patch_counts=dict(baseline),accepted_baseline_patch_counts=dict(accepted),radiance_out_of_range_count=sum(d['source']['pixels_rejected_by_radiance_range'] for d in diags),valid_patch_minimum_texture_std_kR=min(d['variants']['baseline']['texture_std_percentiles_kR'][0] for d in rejected),variant_motion_outcomes=variants,refitted_windows=len(refs),refit_model_classes=dict(Counter(r['selected_model'] for r in refs)),comparison_outcomes=outcomes,common_sample_gain_percent_range=[min(gains),max(gains)] if gains else [],common_sample_gain_percent_median=float(np.median(gains)) if gains else None,promising_starts=[r['start'] for r in promising],consistent_single_slow_starts=[r['start'] for r in records if r['diagnostic_outcome']=='consistent_single_slow'],not_confirmed_starts=[r['start'] for r in records if r['diagnostic_outcome']=='not_confirmed'],not_refitted_starts=sorted(set(d['start'] for d in rejected)-set(r['start'] for r in refs)),notes='No original orbit products changed. Smoothing is used only for even-frame motion estimation; original unsmoothed radiance is fitted and tested. Multiple variants and neighboring overlapping windows make this exploratory evidence, not independent validation. Counts across overlapping windows are not unique pixel counts.')
    savejson(OUT/'summary.json',summary)
    fields=['start','original_class','variant','selected_model','comparison_baseline','common_heldout_samples','common_fraction_of_observed_holdout','improvement_over_original_or_stationary','has_consistent_neighbor','diagnostic_outcome']
    with (OUT/'window_comparisons.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows({k:r[k] for k in fields} for r in records)
    fig,ax=plt.subplots(2,2,figsize=(15,10),layout='constrained')
    reasons=['reference_missing','texture_below_threshold','no_valid_search_position','correlation_below_threshold','boundary_only_patch','patches_with_matches'];labels=['Missing reference pixels','Texture <0.08 kR','No valid search position','NCC ≤0.35','Boundary-only matches','Matched patches']
    vals=[baseline[k]/baseline['patches_considered']*100 for k in reasons];ax[0,0].barh(labels,vals,color=['#a9a9a9','#e6a04b','#999','#df7a65','#aaa','#388a75']);ax[0,0].invert_yaxis();ax[0,0].set(xlabel='% of all considered patches',title='Why matching failed in the 36 rejected windows')
    reasons2=['two_modes','single_slow_mode','only_non_slow_mode','no_cluster_with_15_candidates','fewer_than_20_candidates'];colors=['#6c63bd','#3b9c75','#cd9b4d','#b5b5b5','#dedede'];bottom=np.zeros(4)
    for reason,color in zip(reasons2,colors):
        v=[variants[x].get(reason,0) for x in VARIANTS];ax[0,1].bar(np.arange(4),v,bottom=bottom,color=color,label=reason.replace('_',' '));bottom+=v
    ax[0,1].set(xticks=np.arange(4),xticklabels=['Original','Lower texture\ncutoff','1-pixel\nsmoothing','41×41\npatches'],ylabel='Windows',title='Controlled matching tests (same NCC and cluster gates)');ax[0,1].legend(fontsize=8,loc='lower left')
    if records:
        for r in records:
            col='#3b9c75' if r['diagnostic_outcome']=='promising' else '#6c63bd' if r['diagnostic_outcome']=='consistent_single_slow' else '#cd9b4d';ax[1,0].scatter(r['start'],100*r['improvement_over_original_or_stationary'],color=col)
        ax[1,0].axhline(5,color='#555',ls='--');ax[1,0].axhline(0,color='#999',lw=.7);ax[1,0].set(xlabel='Window start frame',ylabel='Held-out RMSE improvement (%)',title='Green: improved two-motion | purple: consistent slow-only')
        x=[r['rmse_kR']['original'] for r in records];y=[r['rmse_kR'][r['variant']] for r in records];ax[1,1].scatter(x,y,c=['#3b9c75' if r['diagnostic_outcome']=='promising' else '#6c63bd' if r['diagnostic_outcome']=='consistent_single_slow' else '#cd9b4d' for r in records]);bounds=[min(x+y)*.9,max(x+y)*1.1];ax[1,1].plot(bounds,bounds,'--',color='#777');ax[1,1].set(xlabel='Original model or stationary median RMSE (kR)',ylabel='Diagnostic model RMSE (kR)',title='Below diagonal = improved withheld-frame prediction')
    fig.suptitle('Orbit 4 at 85 km | diagnosis of 36 rejected windows',fontsize=17);fig.savefig(OUT/'rejection_diagnosis.png',dpi=160);plt.close(fig)
    lines=['# Rejected-window diagnosis','',f"Audited all 76 windows, including the 36 originally rejected. Original orbit outputs were not changed.",'','## Main findings','',f"No finite radiance samples were excluded by the allowed radiance range. No fully valid baseline patch failed the 0.08 kR texture cutoff; the smallest valid patch standard deviation was {summary['valid_patch_minimum_texture_std_kR']:.3f} kR.",'',f"Of {baseline['patches_considered']:,} considered patches in rejected windows, {baseline['reference_missing']:,} lacked a complete valid reference. Of {baseline['patches_searched']:,} searched patches, {baseline['correlation_below_threshold']:,} failed NCC >0.35 and {baseline['no_valid_search_position']:,} had no fully valid search position. These are repeated patch/frame-pair counts, not unique pixels.",'', 'Fifteen windows had one supported non-slow motion (about 0.63–0.75 cells/s), rather than no motion. Eleven had no 15-candidate cluster; ten had fewer than 20 retained candidates.','', 'Lowering the texture threshold from 0.08 to 0.04 kR did not change the motion outcomes. Larger 41×41 patches generally reduced matching evidence. Gaussian smoothing with sigma=1 pixel, before the existing high-pass filter, recovered a single slow or two-mode candidate in 33 of 36 windows, without lowering the correlation or cluster-support gates.','',f"Refitted {len(refs)} diagnostic windows using ORIGINAL UNSMOOTHED radiance, at 32 iterations for two-layer models. Outcome counts: {summary['refit_model_classes']}.",'',f"On identical supported odd-frame samples, median improvement over the original available model (or a stationary even-frame median where no model existed) was {summary['common_sample_gain_percent_median']:.1f}%. {len(promising)} windows passed the exploratory confirmation rule: accepted single/two model, >5% improvement, and at least one neighboring same-component overlap correlation >=0.5 with velocity difference <=0.15 cells/s.",'','## Limits and recommendation','','- Do not widen the radiance acceptance range or lower the texture gate globally: these were not the bottleneck.','- Mild smoothing for MOTION ESTIMATION is the useful candidate change. Preserve original radiance for layer fitting and residual/mask calculations.','- The neighbor checks involve overlapping observations and are not independent validation. The odd-frame comparisons also inform parameter evaluation, so the results are exploratory.','- A stationary median is a weak comparator in formerly unfit windows. Improvements against it show predictive value but do not by themselves prove correct physical layer separation.','- Unresolved and weakly confirmed cases remain flagged. No cloud-free classification is inferred from a failure to find a second mode.','- The previously discussed local explained-energy cloud mask was not part of this diagnostic; it remains a separate missing product.','','## Per-window comparisons','','| Start | Original outcome | New selected model | Common-sample RMSE gain | Consistent neighbor | Diagnostic result |','|---:|---|---|---:|---|---|']
    for r in sorted(records,key=lambda r:r['start']):lines.append(f"| {r['start']} | {r['original_class']} | {r['selected_model']} | {r['improvement_over_original_or_stationary']:.1%} | {r['has_consistent_neighbor']} | {r['diagnostic_outcome']} |")
    lines+=['',f"No acceptable candidate for refitting at starts: {summary['not_refitted_starts']}.",'', 'Files: `matching_summary.json` contains all rejection counters and contrast/noise statistics; `heldout_comparison.json` contains same-sample errors and overlap checks; `window_comparisons.csv` is the compact table. `fits_smooth_1px/` holds separate experimental fits.','', 'Reproduce: `python diagnose_orbit_windows.py scan`, `python diagnose_orbit_windows.py refit`, then `python report_orbit_diagnostics.py`.']
    lines+=['','## Slow-only interpretation and controls','','Fourteen diagnostic fits are classified as consistent single-slow scenes: their neighboring slow templates agree and their error is slightly below the stationary median (0.4–1.7% improvement). This is distinct from the >5% improvement test used to identify a substantial predictive gain. A nearly stationary single layer need not beat a stationary comparator by 5%; no fast layer is forced. This descriptive classification is not proof of cloud absence.','', 'The numerical controls recover a known single slow motion after smoothing and find no supported mode in three independent-noise sequences. See `verification.json`. These limited controls do not establish physical layer identity.']
    (OUT/'README.md').write_text('\n'.join(lines),encoding='utf-8');print(json.dumps(summary,indent=2),flush=True)

if __name__=='__main__':evaluate()
