# Whole-orbit separation at 85 km

Processed all 1908 source frames through 76 overlapping windows. Supported estimates cover 99.4% of valid L3C pixels. The remainder were processed but lack accepted motion fits or adequate spatial/temporal support.

Window classifications: {'ambiguous_motion': 2, 'single_slow': 38, 'two_motion': 36}

## Interpretation

- One slow motion is a valid outcome. Fast=0 there means the selected model contains no fast layer; it does not prove cloud absence.
- Two-motion estimates require >5% held-out improvement and evidence from at least three frame pairs per mode. Motion labels alone do not establish physical cloud/airglow identity.
- Every 50-frame window uses even local frames for motion estimation and fitting and odd frames for model selection. Adjacent frames are correlated; these are exploratory selection diagnostics, not independent validation.
- Fits use 32 iterations for this first orbit pass; the pilot showed small additional error reduction from 32 to 64. Saved histories permit targeted refinement.
- Actual timestamps and frame gaps are preserved. Small enclosed radiance holes are interpolated for fitting; originals are retained for evaluation. No spatial regridding beyond the source L1R Mastergrid is performed.
- Maps are evaluated at L3C Pixel_Epoch. Existing L3C is a coadded product, so these are aligned component estimates, not an official reprocessing of the L3C coadd.
- Overlapping windows are weighted toward their reference times and require at least eight training-frame contributions. Different model classes and disconnected offset groups are not blended together.
- Brightness offsets were aligned in 9 connected groups using same-time overlaps. Fast gains an offset and slow loses it, preserving their sum. Absolute offsets remain arbitrary between disconnected groups.
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
| 25 | 74 | two_motion | 0.3287 | 15.3% |
| 50 | 99 | two_motion | 0.3285 | 21.7% |
| 75 | 124 | two_motion | 0.3307 | 20.8% |
| 100 | 149 | two_motion | 0.3269 | 15.9% |
| 125 | 174 | two_motion | 0.3170 | 12.2% |
| 150 | 199 | two_motion | 0.3128 | 8.3% |
| 175 | 224 | single_slow | 0.3381 | 0.0% |
| 200 | 249 | two_motion | 0.3200 | 8.4% |
| 225 | 274 | two_motion | 0.3129 | 11.4% |
| 250 | 299 | single_slow | 0.3507 | 0.0% |
| 275 | 324 | single_slow | 0.3264 | 0.0% |
| 300 | 349 | single_slow | 0.5657 | 0.0% |
| 325 | 374 | single_slow | 0.3193 | 0.0% |
| 350 | 399 | single_slow | 0.3381 | 0.0% |
| 375 | 424 | single_slow | 0.3722 | 0.0% |
| 400 | 449 | single_slow | 0.4104 | 0.0% |
| 425 | 474 | single_slow | 0.4813 | 0.0% |
| 450 | 499 | two_motion | 0.3483 | 34.4% |
| 475 | 524 | single_slow | 0.5093 | 0.0% |
| 500 | 549 | single_slow | 0.4533 | 0.0% |
| 525 | 574 | single_slow | 0.3979 | 0.0% |
| 550 | 599 | two_motion | 0.3131 | 26.3% |
| 575 | 624 | two_motion | 0.3218 | 28.2% |
| 600 | 649 | ambiguous_motion | 0.4466 | 0.0% |
| 625 | 674 | two_motion | 0.3341 | 16.6% |
| 650 | 699 | two_motion | 0.3108 | 7.2% |
| 675 | 724 | single_slow | 0.3107 | 4.8% |
| 700 | 749 | two_motion | 0.3130 | 5.9% |
| 725 | 774 | single_slow | 0.3692 | 0.0% |
| 750 | 799 | single_slow | 0.3914 | 0.0% |
| 775 | 824 | two_motion | 0.3386 | 17.7% |
| 800 | 849 | two_motion | 0.3654 | 30.2% |
| 825 | 874 | two_motion | 0.3719 | 35.5% |
| 850 | 899 | two_motion | 0.3662 | 33.2% |
| 875 | 924 | two_motion | 0.3371 | 25.1% |
| 900 | 949 | ambiguous_motion | 0.4046 | 0.0% |
| 925 | 974 | two_motion | 0.3126 | 19.0% |
| 950 | 999 | two_motion | 0.3060 | 10.3% |
| 975 | 1024 | single_slow | 0.3054 | 0.0% |
| 1000 | 1049 | single_slow | 0.3019 | 0.0% |
| 1025 | 1074 | single_slow | 0.3037 | 0.0% |
| 1050 | 1099 | single_slow | 0.3011 | 0.0% |
| 1075 | 1124 | single_slow | 0.3097 | 0.0% |
| 1100 | 1149 | single_slow | 0.3307 | 0.0% |
| 1125 | 1174 | two_motion | 0.3130 | 34.5% |
| 1150 | 1199 | two_motion | 0.3145 | 46.5% |
| 1175 | 1224 | two_motion | 0.3147 | 41.6% |
| 1200 | 1249 | two_motion | 0.3129 | 28.9% |
| 1225 | 1274 | two_motion | 0.3150 | 24.4% |
| 1250 | 1299 | two_motion | 0.3092 | 27.5% |
| 1275 | 1324 | two_motion | 0.3116 | 16.7% |
| 1300 | 1349 | two_motion | 0.3103 | 6.7% |
| 1325 | 1374 | single_slow | 0.3097 | 4.9% |
| 1350 | 1399 | two_motion | 0.3064 | 7.1% |
| 1375 | 1424 | two_motion | 0.3114 | 14.9% |
| 1400 | 1449 | two_motion | 0.3094 | 18.4% |
| 1425 | 1474 | two_motion | 0.3064 | 13.7% |
| 1450 | 1499 | single_slow | 0.3154 | 0.0% |
| 1475 | 1524 | single_slow | 0.3141 | 0.0% |
| 1500 | 1549 | single_slow | 0.3130 | 0.3% |
| 1525 | 1574 | single_slow | 0.3140 | 0.0% |
| 1550 | 1599 | single_slow | 0.3251 | 0.0% |
| 1575 | 1624 | single_slow | 0.3306 | 0.0% |
| 1600 | 1649 | single_slow | 0.3542 | 0.0% |
| 1625 | 1674 | single_slow | 0.3662 | 0.0% |
| 1650 | 1699 | single_slow | 0.4086 | 2.2% |
| 1675 | 1724 | two_motion | 0.3769 | 10.8% |
| 1700 | 1749 | two_motion | 0.4021 | 12.2% |
| 1725 | 1774 | single_slow | 0.4625 | 0.0% |
| 1750 | 1799 | single_slow | 0.5239 | 0.0% |
| 1775 | 1824 | single_slow | 0.5948 | 0.0% |
| 1800 | 1849 | single_slow | 0.7173 | 0.0% |
| 1825 | 1874 | single_slow | 0.7832 | 0.0% |
| 1850 | 1899 | single_slow | 0.6518 | 0.0% |
| 1858 | 1907 | single_slow | 0.6335 | 0.0% |

## Acceptance criteria and checks

A supported motion mode requires at least 15 matching candidates within 0.085 native cells/s of its peak. A single mode is labeled slow only below 0.18 cells/s with at least three supporting frame pairs. These thresholds are heuristic. A well-supported lone faster mode is classified as ambiguous motion, rather than assuming it is cloud-free.

Run `python verify_orbit_full.py` for numerical and product-integrity checks; the completed results are in `verification.json`.

Reproduce: `python orbit_full.py fit --workers 6 --iterations 32`, then `python orbit_full.py assemble`. Completed windows with matching iteration settings are reused. Source files and pilot outputs are unchanged.

## Small-hole processing

Only enclosed missing components of at most 16 pixels and at most 5x5 cells, with valid geometry, are filled by harmonic interpolation. Original finite radiance is unchanged. Larger holes and footprint edges remain missing. Motion matching uses Gaussian smoothing with sigma 1 pixels; templates use the filled, unsmoothed radiance. Evaluation and minimum support use genuine observations. Each window saves an interpolation mask.

Reproduce this run: `python orbit_full.py fit --fill-small-gaps --workers 6 --iterations 32`, then `python orbit_full.py assemble --fill-small-gaps`. Verify with `python verify_orbit_full.py --fill-small-gaps`.