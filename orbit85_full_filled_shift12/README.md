# Whole-orbit separation at 85 km

Processed all 1908 source frames through 76 overlapping windows. Supported estimates cover 98.1% of valid L3C pixels. The remainder were processed but lack accepted motion fits or adequate spatial/temporal support.

Window classifications: {'ambiguous_motion': 4, 'single_slow': 38, 'two_motion': 34}

## Interpretation

- One slow motion is a valid outcome. Fast=0 there means the selected model contains no fast layer; it does not prove cloud absence.
- Two-motion estimates require >5% held-out improvement and evidence from at least three frame pairs per mode. Motion labels alone do not establish physical cloud/airglow identity.
- Every 50-frame window uses even local frames for motion estimation and fitting and odd frames for model selection. Adjacent frames are correlated; these are exploratory selection diagnostics, not independent validation.
- Fits use 32 iterations for this first orbit pass; the pilot showed small additional error reduction from 32 to 64. Saved histories permit targeted refinement.
- Actual timestamps and frame gaps are preserved. Small enclosed radiance holes are interpolated for fitting; originals are retained for evaluation. No spatial regridding beyond the source L1R Mastergrid is performed.
- Maps are evaluated at L3C Pixel_Epoch. Existing L3C is a coadded product, so these are aligned component estimates, not an official reprocessing of the L3C coadd.
- Overlapping windows are weighted toward their reference times and require at least eight training-frame contributions. Different model classes and disconnected offset groups are not blended together.
- Brightness offsets were aligned in 8 connected groups using same-time overlaps. Fast gains an offset and slow loses it, preserving their sum. Absolute offsets remain arbitrary between disconnected groups.
- Gray missing regions are insufficient evidence/support, not cloud-free areas.

## Outputs

- `l3c_cloud_whole_orbit.png`: orbit-wide L3C, cloud score, slow/fast estimates, and coverage status.
- `segment_1.png` through `segment_4.png`: enlarged views.
- `orbit85_separation.nc`: geolocated estimates, status, timestamps and window provenance.
- `window_summary.json`: all window outcomes and held-out histories.
- `offset_alignment.json`: overlap evidence, solved offsets and connected groups.
- `summary.json`: frame and pixel coverage counts.

## Window results

| Start | End | Selected model | Holdout RMSE, final candidate | Improvement vs slow baseline |
|---:|---:|---|---:|---:|
| 0 | 49 | two_motion | 0.3210 | 11.1% |
| 12 | 61 | two_motion | 0.3273 | 12.6% |
| 37 | 86 | two_motion | 0.3277 | 19.8% |
| 62 | 111 | two_motion | 0.3300 | 21.8% |
| 87 | 136 | two_motion | 0.3320 | 17.7% |
| 112 | 161 | two_motion | 0.3230 | 13.3% |
| 137 | 186 | two_motion | 0.3157 | 10.8% |
| 162 | 211 | single_slow | 0.3359 | 0.0% |
| 187 | 236 | single_slow | 0.3383 | 0.0% |
| 212 | 261 | two_motion | 0.3159 | 10.9% |
| 237 | 286 | two_motion | 0.3118 | 10.1% |
| 262 | 311 | single_slow | 0.3438 | 0.0% |
| 287 | 336 | single_slow | 0.3213 | 0.0% |
| 312 | 361 | single_slow | 0.5651 | 0.0% |
| 337 | 386 | single_slow | 0.3274 | 0.0% |
| 362 | 411 | single_slow | 0.3518 | 0.0% |
| 387 | 436 | single_slow | 0.3876 | 0.0% |
| 412 | 461 | single_slow | 0.4362 | 0.0% |
| 437 | 486 | single_slow | 0.5126 | 0.0% |
| 462 | 511 | single_slow | 0.5296 | 0.0% |
| 487 | 536 | single_slow | 0.4854 | 0.0% |
| 512 | 561 | single_slow | 0.4230 | 0.0% |
| 537 | 586 | two_motion | 0.3103 | 22.5% |
| 562 | 611 | two_motion | 0.3186 | 27.6% |
| 587 | 636 | ambiguous_motion | 0.4567 | 0.0% |
| 612 | 661 | two_motion | 0.3286 | 21.5% |
| 637 | 686 | two_motion | 0.3273 | 12.2% |
| 662 | 711 | single_slow | 0.3095 | 5.0% |
| 687 | 736 | single_slow | 0.3103 | 5.0% |
| 712 | 761 | single_slow | 0.3540 | 0.0% |
| 737 | 786 | two_motion | 0.3244 | 14.9% |
| 762 | 811 | two_motion | 0.3329 | 17.0% |
| 787 | 836 | two_motion | 0.3545 | 24.0% |
| 812 | 861 | two_motion | 0.3707 | 33.8% |
| 837 | 886 | ambiguous_motion | 0.4544 | 0.0% |
| 862 | 911 | ambiguous_motion | 0.4375 | 0.0% |
| 887 | 936 | two_motion | 0.3230 | 20.8% |
| 912 | 961 | ambiguous_motion | 0.4032 | 0.0% |
| 937 | 986 | two_motion | 0.3116 | 16.0% |
| 962 | 1011 | single_slow | 0.3221 | 0.0% |
| 987 | 1036 | single_slow | 0.3028 | 0.0% |
| 1012 | 1061 | single_slow | 0.3038 | 0.0% |
| 1037 | 1086 | single_slow | 0.3040 | 0.0% |
| 1062 | 1111 | single_slow | 0.3022 | 0.0% |
| 1087 | 1136 | single_slow | 0.3128 | 0.0% |
| 1112 | 1161 | two_motion | 0.3101 | 20.8% |
| 1137 | 1186 | two_motion | 0.3172 | 41.7% |
| 1162 | 1211 | two_motion | 0.3150 | 45.7% |
| 1187 | 1236 | two_motion | 0.3128 | 35.7% |
| 1212 | 1261 | two_motion | 0.3107 | 24.0% |
| 1237 | 1286 | two_motion | 0.3151 | 27.0% |
| 1262 | 1311 | two_motion | 0.3092 | 23.3% |
| 1287 | 1336 | two_motion | 0.3114 | 10.2% |
| 1312 | 1361 | single_slow | 0.3085 | 4.5% |
| 1337 | 1386 | two_motion | 0.3109 | 5.3% |
| 1362 | 1411 | two_motion | 0.3104 | 10.6% |
| 1387 | 1436 | two_motion | 0.3100 | 17.4% |
| 1412 | 1461 | two_motion | 0.3068 | 17.5% |
| 1437 | 1486 | two_motion | 0.3053 | 7.4% |
| 1462 | 1511 | single_slow | 0.3107 | 0.3% |
| 1487 | 1536 | single_slow | 0.3157 | 0.0% |
| 1512 | 1561 | single_slow | 0.3137 | 0.0% |
| 1537 | 1586 | single_slow | 0.3186 | 0.0% |
| 1562 | 1611 | single_slow | 0.3274 | 0.0% |
| 1587 | 1636 | single_slow | 0.3283 | 2.0% |
| 1612 | 1661 | single_slow | 0.3675 | 4.3% |
| 1637 | 1686 | single_slow | 0.3778 | 0.0% |
| 1662 | 1711 | single_slow | 0.3970 | 3.0% |
| 1687 | 1736 | two_motion | 0.3779 | 14.2% |
| 1712 | 1761 | two_motion | 0.4230 | 9.8% |
| 1737 | 1786 | single_slow | 0.4780 | 0.0% |
| 1762 | 1811 | single_slow | 0.5606 | 0.0% |
| 1787 | 1836 | single_slow | 0.6412 | 0.0% |
| 1812 | 1861 | single_slow | 0.7310 | 0.0% |
| 1837 | 1886 | single_slow | 0.7903 | 0.0% |
| 1858 | 1907 | single_slow | 0.6335 | 0.0% |

## Acceptance criteria and checks

A supported motion mode requires at least 15 matching candidates within 0.085 native cells/s of its peak. A single mode is labeled slow only below 0.18 cells/s with at least three supporting frame pairs. These thresholds are heuristic. A well-supported lone faster mode is classified as ambiguous motion, rather than assuming it is cloud-free.

Run `python verify_orbit_full.py` for numerical and product-integrity checks; the completed results are in `verification.json`.

Reproduce: `python orbit_full.py fit --workers 6 --iterations 32`, then `python orbit_full.py assemble`. Completed windows with matching iteration settings are reused. Source files and pilot outputs are unchanged.

## Small-hole processing

Only enclosed missing components of at most 16 pixels and at most 5x5 cells, with valid geometry, are filled by harmonic interpolation. Original finite radiance is unchanged. Larger holes and footprint edges remain missing. Motion matching uses Gaussian smoothing with sigma 1 pixels; templates use the filled, unsmoothed radiance. Evaluation and minimum support use genuine observations. Each window saves an interpolation mask.

Reproduce this run: `python orbit_full.py fit --fill-small-gaps --start-offset 12 --workers 6 --iterations 32`, then `python orbit_full.py assemble --fill-small-gaps --start-offset 12`. Verify with `python verify_orbit_full.py --fill-small-gaps --start-offset 12`.

Regular window starts are shifted by 12 frames. Endpoint starts 0 and 1858 are retained; they are shared with the original pass and are not independent comparisons.