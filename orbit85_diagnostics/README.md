# Rejected-window diagnosis

Audited all 76 windows, including the 36 originally rejected. Original orbit outputs were not changed.

## Main findings

No finite radiance samples were excluded by the allowed radiance range. No fully valid baseline patch failed the 0.08 kR texture cutoff; the smallest valid patch standard deviation was 0.225 kR.

Of 39,096 considered patches in rejected windows, 17,160 lacked a complete valid reference. Of 21,936 searched patches, 16,910 failed NCC >0.35 and 2,594 had no fully valid search position. These are repeated patch/frame-pair counts, not unique pixels.

Fifteen windows had one supported non-slow motion (about 0.63–0.75 cells/s), rather than no motion. Eleven had no 15-candidate cluster; ten had fewer than 20 retained candidates.

Lowering the texture threshold from 0.08 to 0.04 kR did not change the motion outcomes. Larger 41×41 patches generally reduced matching evidence. Gaussian smoothing with sigma=1 pixel, before the existing high-pass filter, recovered a single slow or two-mode candidate in 33 of 36 windows, without lowering the correlation or cluster-support gates.

Refitted 33 diagnostic windows using ORIGINAL UNSMOOTHED radiance, at 32 iterations for two-layer models. Outcome counts: {'two_motion': 19, 'single_slow': 14}.

On identical supported odd-frame samples, median improvement over the original available model (or a stationary even-frame median where no model existed) was 12.6%. 19 windows passed the exploratory confirmation rule: accepted single/two model, >5% improvement, and at least one neighboring same-component overlap correlation >=0.5 with velocity difference <=0.15 cells/s.

## Limits and recommendation

- Do not widen the radiance acceptance range or lower the texture gate globally: these were not the bottleneck.
- Mild smoothing for MOTION ESTIMATION is the useful candidate change. Preserve original radiance for layer fitting and residual/mask calculations.
- The neighbor checks involve overlapping observations and are not independent validation. The odd-frame comparisons also inform parameter evaluation, so the results are exploratory.
- A stationary median is a weak comparator in formerly unfit windows. Improvements against it show predictive value but do not by themselves prove correct physical layer separation.
- Unresolved and weakly confirmed cases remain flagged. No cloud-free classification is inferred from a failure to find a second mode.
- The previously discussed local explained-energy cloud mask was not part of this diagnostic; it remains a separate missing product.

## Per-window comparisons

| Start | Original outcome | New selected model | Common-sample RMSE gain | Consistent neighbor | Diagnostic result |
|---:|---|---|---:|---|---|
| 0 | ambiguous_motion | two_motion | 18.6% | True | promising |
| 25 | ambiguous_motion | two_motion | 19.2% | True | promising |
| 50 | ambiguous_motion | two_motion | 17.6% | True | promising |
| 100 | insufficient_evidence | two_motion | 18.4% | True | promising |
| 125 | insufficient_evidence | two_motion | 14.3% | True | promising |
| 150 | insufficient_evidence | two_motion | 10.4% | True | promising |
| 575 | ambiguous_motion | two_motion | 32.0% | True | promising |
| 625 | ambiguous_motion | two_motion | 23.4% | True | promising |
| 650 | ambiguous_motion | two_motion | 24.5% | True | promising |
| 675 | insufficient_evidence | single_slow | 1.3% | True | consistent_single_slow |
| 825 | ambiguous_motion | two_motion | 20.5% | True | promising |
| 850 | ambiguous_motion | two_motion | 18.3% | True | promising |
| 875 | ambiguous_motion | two_motion | 20.8% | True | promising |
| 925 | ambiguous_motion | two_motion | 21.7% | True | promising |
| 950 | ambiguous_motion | two_motion | 26.2% | True | promising |
| 1000 | insufficient_evidence | single_slow | 1.4% | True | consistent_single_slow |
| 1025 | insufficient_evidence | single_slow | 1.7% | True | consistent_single_slow |
| 1275 | ambiguous_motion | two_motion | 30.5% | True | promising |
| 1300 | insufficient_evidence | single_slow | 0.4% | True | consistent_single_slow |
| 1325 | insufficient_evidence | single_slow | 1.3% | True | consistent_single_slow |
| 1350 | insufficient_evidence | two_motion | 8.4% | True | promising |
| 1375 | insufficient_evidence | two_motion | 15.5% | True | promising |
| 1400 | ambiguous_motion | two_motion | 15.7% | True | promising |
| 1450 | insufficient_evidence | single_slow | 1.6% | True | consistent_single_slow |
| 1475 | insufficient_evidence | single_slow | 1.1% | True | consistent_single_slow |
| 1500 | insufficient_evidence | single_slow | 1.1% | True | consistent_single_slow |
| 1525 | insufficient_evidence | single_slow | 1.4% | True | consistent_single_slow |
| 1550 | insufficient_evidence | single_slow | 1.1% | True | consistent_single_slow |
| 1575 | insufficient_evidence | single_slow | 1.0% | True | consistent_single_slow |
| 1600 | insufficient_evidence | single_slow | 1.3% | True | consistent_single_slow |
| 1625 | insufficient_evidence | single_slow | 0.4% | True | consistent_single_slow |
| 1650 | insufficient_evidence | single_slow | 0.8% | True | consistent_single_slow |
| 1675 | insufficient_evidence | two_motion | 12.6% | True | promising |

No acceptable candidate for refitting at starts: [600, 900, 1425].

Files: `matching_summary.json` contains all rejection counters and contrast/noise statistics; `heldout_comparison.json` contains same-sample errors and overlap checks; `window_comparisons.csv` is the compact table. `fits_smooth_1px/` holds separate experimental fits.

Reproduce: `python diagnose_orbit_windows.py scan`, `python diagnose_orbit_windows.py refit`, then `python report_orbit_diagnostics.py`.

## Slow-only interpretation and controls

Fourteen diagnostic fits are classified as consistent single-slow scenes: their neighboring slow templates agree and their error is slightly below the stationary median (0.4–1.7% improvement). This is distinct from the >5% improvement test used to identify a substantial predictive gain. A nearly stationary single layer need not beat a stationary comparator by 5%; no fast layer is forced. This descriptive classification is not proof of cloud absence.

The numerical controls recover a known single slow motion after smoothing and find no supported mode in three independent-noise sequences. See `verification.json`. These limited controls do not establish physical layer identity.