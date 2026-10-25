# Whole-orbit separation at 85 km

Processed all 1908 source frames through 76 overlapping windows. Supported estimates cover 58.9% of valid L3C pixels. The remainder were processed but lack accepted motion fits or adequate spatial/temporal support.

Window classifications: {'ambiguous_motion': 15, 'insufficient_evidence': 21, 'single_slow': 22, 'two_motion': 18}

## Interpretation

- One slow motion is a valid outcome. Fast=0 there means the selected model contains no fast layer; it does not prove cloud absence.
- Two-motion estimates require >5% held-out improvement and evidence from at least three frame pairs per mode. Motion labels alone do not establish physical cloud/airglow identity.
- Every 50-frame window uses even local frames for motion estimation and fitting and odd frames for model selection. Adjacent frames are correlated; these are exploratory selection diagnostics, not independent validation.
- Fits use 32 iterations for this first orbit pass; the pilot showed small additional error reduction from 32 to 64. Saved histories permit targeted refinement.
- Actual timestamps, missing radiance, and frame gaps are preserved. No spatial regridding beyond the source L1R Mastergrid is performed.
- Maps are evaluated at L3C Pixel_Epoch. Existing L3C is a coadded product, so these are aligned component estimates, not an official reprocessing of the L3C coadd.
- Overlapping windows are weighted toward their reference times and require at least eight training-frame contributions. Different model classes and disconnected offset groups are not blended together.
- Brightness offsets were aligned in 7 connected groups using same-time overlaps. Fast gains an offset and slow loses it, preserving their sum. Absolute offsets remain arbitrary between disconnected groups.
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
| 0 | 49 | ambiguous_motion | 0.3991 | 0.0% |
| 25 | 74 | ambiguous_motion | 0.4015 | 0.0% |
| 50 | 99 | ambiguous_motion | 0.4005 | 0.0% |
| 75 | 124 | two_motion | 0.3289 | 21.5% |
| 100 | 149 | insufficient_evidence | — | — |
| 125 | 174 | insufficient_evidence | — | — |
| 150 | 199 | insufficient_evidence | — | — |
| 175 | 224 | single_slow | 0.3387 | 0.0% |
| 200 | 249 | two_motion | 0.3202 | 8.3% |
| 225 | 274 | two_motion | 0.3128 | 11.4% |
| 250 | 299 | single_slow | 0.3521 | 0.0% |
| 275 | 324 | single_slow | 0.3284 | 0.0% |
| 300 | 349 | single_slow | 0.5665 | 0.0% |
| 325 | 374 | single_slow | 0.3195 | 0.0% |
| 350 | 399 | single_slow | 0.3380 | 0.0% |
| 375 | 424 | single_slow | 0.3720 | 0.0% |
| 400 | 449 | single_slow | 0.4100 | 0.0% |
| 425 | 474 | single_slow | 0.4810 | 0.0% |
| 450 | 499 | two_motion | 0.3483 | 34.4% |
| 475 | 524 | two_motion | 0.3485 | 31.6% |
| 500 | 549 | single_slow | 0.4534 | 0.0% |
| 525 | 574 | single_slow | 0.3982 | 0.0% |
| 550 | 599 | two_motion | 0.3132 | 26.3% |
| 575 | 624 | ambiguous_motion | 0.4716 | 0.0% |
| 600 | 649 | ambiguous_motion | 0.4474 | 0.0% |
| 625 | 674 | ambiguous_motion | 0.4355 | 0.0% |
| 650 | 699 | ambiguous_motion | 0.4100 | 0.0% |
| 675 | 724 | insufficient_evidence | — | — |
| 700 | 749 | two_motion | 0.3129 | 5.9% |
| 725 | 774 | two_motion | 0.3214 | 13.0% |
| 750 | 799 | two_motion | 0.3307 | 15.4% |
| 775 | 824 | two_motion | 0.3384 | 17.7% |
| 800 | 849 | two_motion | 0.3661 | 30.0% |
| 825 | 874 | ambiguous_motion | 0.4669 | 0.0% |
| 850 | 899 | ambiguous_motion | 0.4477 | 0.0% |
| 875 | 924 | ambiguous_motion | 0.4249 | 0.0% |
| 900 | 949 | ambiguous_motion | 0.4066 | 0.0% |
| 925 | 974 | ambiguous_motion | 0.3984 | 0.0% |
| 950 | 999 | ambiguous_motion | 0.4138 | 0.0% |
| 975 | 1024 | single_slow | 0.3081 | 0.0% |
| 1000 | 1049 | insufficient_evidence | — | — |
| 1025 | 1074 | insufficient_evidence | — | — |
| 1050 | 1099 | single_slow | 0.3019 | 0.0% |
| 1075 | 1124 | single_slow | 0.3108 | 0.0% |
| 1100 | 1149 | single_slow | 0.3308 | 0.0% |
| 1125 | 1174 | two_motion | 0.3130 | 34.4% |
| 1150 | 1199 | two_motion | 0.3145 | 46.5% |
| 1175 | 1224 | two_motion | 0.3137 | 42.1% |
| 1200 | 1249 | two_motion | 0.3141 | 28.5% |
| 1225 | 1274 | two_motion | 0.3149 | 24.4% |
| 1250 | 1299 | two_motion | 0.3088 | 27.5% |
| 1275 | 1324 | ambiguous_motion | 0.4458 | 0.0% |
| 1300 | 1349 | insufficient_evidence | — | — |
| 1325 | 1374 | insufficient_evidence | — | — |
| 1350 | 1399 | insufficient_evidence | — | — |
| 1375 | 1424 | insufficient_evidence | — | — |
| 1400 | 1449 | ambiguous_motion | 0.3678 | 0.0% |
| 1425 | 1474 | insufficient_evidence | — | — |
| 1450 | 1499 | insufficient_evidence | — | — |
| 1475 | 1524 | insufficient_evidence | — | — |
| 1500 | 1549 | insufficient_evidence | — | — |
| 1525 | 1574 | insufficient_evidence | — | — |
| 1550 | 1599 | insufficient_evidence | — | — |
| 1575 | 1624 | insufficient_evidence | — | — |
| 1600 | 1649 | insufficient_evidence | — | — |
| 1625 | 1674 | insufficient_evidence | — | — |
| 1650 | 1699 | insufficient_evidence | — | — |
| 1675 | 1724 | insufficient_evidence | — | — |
| 1700 | 1749 | two_motion | 0.4015 | 11.8% |
| 1725 | 1774 | single_slow | 0.4622 | 0.0% |
| 1750 | 1799 | single_slow | 0.5243 | 0.0% |
| 1775 | 1824 | single_slow | 0.5951 | 0.0% |
| 1800 | 1849 | single_slow | 0.7174 | 0.0% |
| 1825 | 1874 | single_slow | 0.7836 | 0.0% |
| 1850 | 1899 | single_slow | 0.6535 | 0.0% |
| 1858 | 1907 | single_slow | 0.6326 | 0.0% |

## Acceptance criteria and checks

A supported motion mode requires at least 15 matching candidates within 0.085 native cells/s of its peak. A single mode is labeled slow only below 0.18 cells/s with at least three supporting frame pairs. These thresholds are heuristic. A well-supported lone faster mode is classified as ambiguous motion, rather than assuming it is cloud-free.

Run `python verify_orbit_full.py` for numerical and product-integrity checks; the completed results are in `verification.json`.

Reproduce: `python orbit_full.py fit --workers 6 --iterations 32`, then `python orbit_full.py assemble`. Completed windows with matching iteration settings are reused. Source files and pilot outputs are unchanged.