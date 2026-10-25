# Orbit 4: 85 km pilot

Inventory: 1908 frames over 35.20 minutes; 12 retained timestamp gaps. Master shape: 551 × 7801. All valid L1R coordinates and L3C coordinates matched the Master within 7.6e-06 degrees (float32 rounding).

## Six-window experiment

| Frames | Classification | Slow-only holdout RMSE | Two-layer RMSE (64) | Improvement | L3C-matched pixels |
|---|---|---:|---:|---:|---:|
| 200–249 | pilot_supported | 0.349 | 0.319 | 8.5% | 29310 |
| 225–274 | pilot_supported | 0.353 | 0.311 | 11.8% | 29290 |
| 1154–1203 | pilot_supported | 0.591 | 0.313 | 47.0% | 31080 |
| 1179–1228 | pilot_supported | 0.524 | 0.313 | 40.2% | 29668 |
| 1650–1699 | unresolved: No supported motion mode was resolved | — | — | — | — |
| 1675–1724 | unresolved: No supported motion mode was resolved | — | — | — | — |

## Method and limits

- Native L1R radiance is placed directly on verified Master cells; no new spatial interpolation. Source files are untouched.
- Each window has 50 observations; only even local indices estimate velocities and train templates. Odd indices evaluate predictions. Neighbouring frames remain temporally correlated; this is an exploratory holdout, not independent atmospheric validation.
- Two strongest distinct supported motion modes are selected without fixed slow/fast speed bins. The slower is labeled slow; neither label proves a physical identity. Full candidates and diagnostic histograms are saved.
- Constant translation is expressed in native along/cross-track grid cells, nominally 2 km. This is a local approximation over each window, not an Earth-fixed velocity retrieval.
- Minimum template support is 8 training frames. Unsupported updates retain previous values; displayed cloud estimates require final support of at least 8 for both layers.
- Held-out sample masks are fixed from iteration 1 onward. Baseline is the initial slow-only temporal-median template at the same slow velocity.
- The classification pilot_supported requires >5% held-out improvement and at least 3 frame pairs per mode; it is a heuristic flag, not calibrated confidence.
- L3C is the existing coadded product. Fast contrast is evaluated at its per-pixel representative Epoch, restricted to the central 25 frames of each window. It is not a reconstruction of the upstream 19-frame coadd.
- MLCloud is the existing percentage-valued product. It is not interchangeable with the signed motion-derived fast radiance.
- Each fast template has median zero; absolute radiance partition remains ambiguous. The orbit pilot map chooses the closest reference time without blending or solving between-window offsets.
- Only the pilot regions are processed; gray areas in its map are unprocessed or unsupported, not cloud-free.

## Overlap checks

- Windows [200, 225]: fast-component correlation 0.864; offset-removed RMSE 0.213 kR; 64176 supported common pixels.
- Windows [1154, 1179]: fast-component correlation 0.949; offset-removed RMSE 0.253 kR; 62635 supported common pixels.

## Files

- `l3c_cloud_and_pilot.png`: orbit comparison.
- `window_*/comparison.png`: detailed L3C, MLCloud, fast contrast and reconstruction panels.
- `heldout_convergence.png`, `overlap_consistency.json`: diagnostics.
- `pilot_on_mastergrid.nc`: geolocated L3C reference and sparse pilot fields, with window provenance.
- `inventory.json`, `frame_inventory.csv`: complete source inventory.
- `window_*/checkpoint_*.npz`: fitted templates at 16/32/64 iterations.

Reproduce using Python 3.12 with NumPy, SciPy, matplotlib and netCDF4:
```
python orbit_pilot.py all
python plot_orbit_pilot.py
```
The scripts use the existing workspace `.packages` directory. Completed windows are reused; remove only an explicitly chosen window output to refit it.

## Single-motion scenes and interpretation

A second layer is optional. The estimator now accepts one supported motion mode and fits a single translated template in that case. A mode below 0.18 native cells/s with evidence from at least three frame pairs is provisionally labeled single slow motion; this speed threshold is heuristic. The displayed fast field is left missing when no fast component is estimated, rather than presenting zeros as measured cloud absence.
A known single-motion synthetic scene passes the one-mode check. In these real pilot windows, the four successful fits have two supported modes. The two late windows do not meet the current minimum of 15 candidates within a 0.085-cell/s neighborhood for even one mode. They therefore remain insufficient-evidence cases, not evidence for cloud presence or absence.
The next whole-orbit pass should retain three outcomes: one supported slow motion, two supported motions with predictive benefit, and insufficient evidence. Resolve brightness offsets between adjacent templates before constructing a continuous fast-contrast mosaic. The current pilot has not performed that stitching.

Verification: run `python verify_orbit_pilot.py`; see `verification.json` for the completed numerical and product checks.