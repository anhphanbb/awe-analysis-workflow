# AWE analysis workflow and experiment archive

Research scripts, notebooks, reports, and parameter experiments for AWE parallax, slow/fast layer separation, whole-orbit processing, and juwavelet analysis. This is an archive of exploratory work, not a validated retrieval product or a portable turnkey package.

Open `index.html` to browse all HTML reports, or start with `AWE_v02_orbit19_052340/orbit19_052340.html`. Source and result directories retain their original names so the progression remains traceable. See [EXPERIMENTS.md](EXPERIMENTS.md).

## Latest orbit 19 result

- L1R Q20 v02, orbit 00019; target 2023-11-23 05:23:40 UTC, nearest frame 1467 at 05:23:40.527 UTC (zero-based).
- 50-frame local fit, 0.02-degree geographic grid, 128 separation iterations.
- juwavelet applied to a 151 by 151 pixel crop of the fitted slow component: longitude -59.4 to -56.4 degrees, latitude -4.8 to -1.8 degrees.
- Wavelengths 12 to 74.026 km, 22 logarithmic scales, strictly below 80 km; scale spacing 0.125 octaves.
- Grid search: coefficient clustering threshold 0.02 / 0.05 and angular steps 5 / 10 / 20 degrees. Default 0.02 and 5 degrees; minimum peak amplitude 0.25 in normalized coefficient units.
- Plane removal and 40 km edge taper. No extra observed-frame or support-count gate on the finite slow component. Missing values already present in L1R remain missing; upstream solar-panel masking was not reversed.
- Leading reconstructed wave scale is approximately 37 km. Clusters are algorithmic coefficient groups, not significance-tested independent physical waves.

## Reproduction

Use Python 3.12 with numpy, scipy, matplotlib, netCDF4, Pillow, numba, and juwavelet. The original run used a local juwavelet checkout; the transform SHA256 is recorded in the wavelet summaries. Install that dependency separately and retain its upstream AGPL license. A dependency list is supplied for convenience; versions are not a locked environment.

Obtain the original AWE L1R/L1F/L3C products and Master geometry files separately. The main orbit runner supports a configurable data root:

```sh
python run_orbit.py --orbit 19 --source l1r --data-root /path/to/soc --workers 6 --iterations 32
```

Local orbit 19 scripts retain experiment-specific paths. Set `SOURCE` and `JUDIR` in `build_orbit19_actual.py` for your local AWE product and juwavelet checkout. Review source paths in the older scripts/notebooks before running them. Execute from the repository root.

```sh
python build_orbit19_actual.py
python run_orbit19_128_grid.py
python render_orbit19_html.py
python render_grid128.py
node verify_orbit19_html.cjs
```

The first command creates the 64-iteration fit and initial whole-slow wavelet result. The second continues that checkpoint to 128 iterations and computes the cropped grid. Some historical renderers and comparison scripts require intermediate arrays or crops from earlier stages; this archive does not claim every old exploratory step is a single-command replay. Consult script inputs and `archive_manifest.json` for omitted arrays. Never run older renderers over a report you intend to preserve without making a copy first.

## Archive contents and limits

All discovered experiment source scripts, notebooks, JSON/CSV summaries, Markdown notes, HTML reports, PNG figures, and GIF animations are included. Raw netCDF/Master files, NumPy caches/checkpoints, ZIP duplicates, logs, and installed runtimes are excluded. `archive_manifest.json` records included and excluded original file paths, sizes, and original hashes for included files. HTML publication copies extract embedded images into adjacent `assets/` folders without changing image bytes, so their HTML hashes differ from the originals.

Earlier reports deliberately keep their original settings and conclusions. The L1F attempt was interrupted, not completed. Separation fits alone do not prove physical layer identity. Local frame fits are not held-out validation. Automated checks cover numerical consistency and report/player logic, not independent scientific validation.

## Attribution

`wave_cluster_limited.py` is derived from the local juwavelet `utils.py`, with an added maximum-cluster stopping limit. It is distributed under GNU AGPL v3; see `JUWAVELET_LICENSE.md`. The original algorithm and coefficient peak ranking are retained. juwavelet implements Morlet continuous wavelet transforms; the upstream README names J. Ungermann as its support contact. Other project files are provided as a research archive; no additional blanket license is asserted here.

## Publishing

The root contains a static report index and `.nojekyll`. It can be served with `python -m http.server` or GitHub Pages from the repository root. The scientific computations run offline in Python; the website displays saved results.
