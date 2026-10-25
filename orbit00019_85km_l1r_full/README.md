# Whole-orbit separation at 85 km

Processed all 1923 source frames through 76 overlapping windows. Supported estimates cover 59.4% of valid L3C pixels. The remainder were processed but lack accepted motion fits or adequate spatial/temporal support.

Window classifications: {'ambiguous_motion': 19, 'insufficient_evidence': 22, 'single_slow': 22, 'two_motion': 13}

## Interpretation

- One slow motion is a valid outcome. Fast=0 there means the selected model contains no fast layer; it does not prove cloud absence.
- Two-motion estimates require >5% held-out improvement and evidence from at least three frame pairs per mode. Motion labels alone do not establish physical cloud/airglow identity.
- Every 50-frame window uses even local frames for motion estimation and fitting and odd frames for model selection. Adjacent frames are correlated; these are exploratory selection diagnostics, not independent validation.
- Fits use 32 iterations for this first orbit pass; the pilot showed small additional error reduction from 32 to 64. Saved histories permit targeted refinement.
- Actual timestamps, missing radiance, and frame gaps are preserved. L1F runs use a floating-point projection onto the Master grid; all diagnostics refer to those resampled observations. No further spatial regridding is performed.
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
| 0 | 49 | single_slow | 0.3567 | 0.0% |
| 25 | 74 | single_slow | 0.3552 | 0.0% |
| 50 | 99 | single_slow | 0.3397 | 2.4% |
| 75 | 124 | single_slow | 0.3500 | 0.0% |
| 100 | 149 | single_slow | 0.3619 | 0.0% |
| 125 | 174 | ambiguous_motion | 0.3706 | 0.0% |
| 150 | 199 | ambiguous_motion | 0.3807 | 0.0% |
| 175 | 224 | ambiguous_motion | 0.3704 | 0.0% |
| 200 | 249 | ambiguous_motion | 0.3645 | 0.0% |
| 225 | 274 | insufficient_evidence | — | — |
| 250 | 299 | single_slow | 0.3332 | 0.0% |
| 275 | 324 | two_motion | 0.3855 | 23.4% |
| 300 | 349 | ambiguous_motion | 0.5316 | 0.0% |
| 325 | 374 | ambiguous_motion | 0.5367 | 0.0% |
| 350 | 399 | ambiguous_motion | 0.5728 | 0.0% |
| 375 | 424 | ambiguous_motion | 0.5057 | 0.0% |
| 400 | 449 | ambiguous_motion | 0.4746 | 0.0% |
| 425 | 474 | single_slow | 0.3381 | 0.0% |
| 450 | 499 | single_slow | 0.3687 | 0.0% |
| 475 | 524 | ambiguous_motion | 0.4361 | 0.0% |
| 500 | 549 | two_motion | 0.3306 | 22.9% |
| 525 | 574 | two_motion | 0.3251 | 22.9% |
| 550 | 599 | single_slow | 0.3976 | 0.0% |
| 575 | 624 | single_slow | 0.3626 | 0.0% |
| 600 | 649 | single_slow | 0.3400 | 0.0% |
| 625 | 674 | ambiguous_motion | 0.4235 | 0.0% |
| 650 | 699 | two_motion | 0.3240 | 12.5% |
| 675 | 724 | two_motion | 0.3371 | 12.8% |
| 700 | 749 | single_slow | 0.3929 | 0.0% |
| 725 | 774 | two_motion | 0.3424 | 17.2% |
| 750 | 799 | ambiguous_motion | 0.4494 | 0.0% |
| 775 | 824 | ambiguous_motion | 0.4922 | 0.0% |
| 800 | 849 | two_motion | 0.3295 | 21.2% |
| 825 | 874 | insufficient_evidence | — | — |
| 850 | 899 | insufficient_evidence | — | — |
| 875 | 924 | insufficient_evidence | — | — |
| 900 | 949 | insufficient_evidence | — | — |
| 925 | 974 | insufficient_evidence | — | — |
| 950 | 999 | insufficient_evidence | — | — |
| 975 | 1024 | insufficient_evidence | — | — |
| 1000 | 1049 | insufficient_evidence | — | — |
| 1025 | 1074 | single_slow | 0.3082 | 0.0% |
| 1050 | 1099 | insufficient_evidence | — | — |
| 1075 | 1124 | insufficient_evidence | — | — |
| 1100 | 1149 | insufficient_evidence | — | — |
| 1125 | 1174 | insufficient_evidence | — | — |
| 1150 | 1199 | insufficient_evidence | — | — |
| 1175 | 1224 | single_slow | 0.3114 | 0.1% |
| 1200 | 1249 | single_slow | 0.3051 | 0.0% |
| 1225 | 1274 | single_slow | 0.3225 | 0.0% |
| 1250 | 1299 | insufficient_evidence | — | — |
| 1275 | 1324 | insufficient_evidence | — | — |
| 1300 | 1349 | insufficient_evidence | — | — |
| 1325 | 1374 | insufficient_evidence | — | — |
| 1350 | 1399 | insufficient_evidence | — | — |
| 1375 | 1424 | insufficient_evidence | — | — |
| 1400 | 1449 | ambiguous_motion | 0.3732 | 0.0% |
| 1425 | 1474 | two_motion | 0.3063 | 30.2% |
| 1450 | 1499 | ambiguous_motion | 0.4094 | 0.0% |
| 1475 | 1524 | two_motion | 0.3132 | 20.4% |
| 1500 | 1549 | ambiguous_motion | 0.4289 | 0.0% |
| 1525 | 1574 | ambiguous_motion | 0.3973 | 0.0% |
| 1550 | 1599 | ambiguous_motion | 0.4571 | 0.0% |
| 1575 | 1624 | ambiguous_motion | 0.4730 | 0.0% |
| 1600 | 1649 | two_motion | 0.3560 | 22.4% |
| 1625 | 1674 | single_slow | 0.4398 | 0.0% |
| 1650 | 1699 | single_slow | 0.4745 | 0.0% |
| 1675 | 1724 | single_slow | 0.4322 | 0.0% |
| 1700 | 1749 | single_slow | 0.5234 | 0.0% |
| 1725 | 1774 | single_slow | 0.5288 | 0.0% |
| 1750 | 1799 | single_slow | 0.6706 | 0.0% |
| 1775 | 1824 | two_motion | 0.5700 | 14.9% |
| 1800 | 1849 | two_motion | 0.5640 | 11.6% |
| 1825 | 1874 | two_motion | 0.6035 | 5.7% |
| 1850 | 1899 | insufficient_evidence | — | — |
| 1873 | 1922 | insufficient_evidence | — | — |

## Acceptance criteria and checks

A supported motion mode requires at least 15 matching candidates within 0.085 native cells/s of its peak. A single mode is labeled slow only below 0.18 cells/s with at least three supporting frame pairs. These thresholds are heuristic. A well-supported lone faster mode is classified as ambiguous motion, rather than assuming it is cloud-free.



Input: L1R; orbit 19.
Reproduce: `python run_orbit.py --orbit 19 --source l1r`.