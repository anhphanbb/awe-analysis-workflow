# Whole-orbit analysis from L1F

Run with Python 3.12 and the workspace `.packages` dependencies:

```powershell
python run_orbit.py --orbit 19 --source l1f --workers 6 --iterations 32
python verify_l1f_adapter.py
```

The runner discovers one v02 Q20 L1F file for the requested orbit, uses the
matching L3C and existing Master grid, prepares the input, inventories every
frame, fits overlapping 50-frame windows, assembles the orbit maps, computes
withheld and all-frame q diagnostics, and verifies the result. Completed
windows are reusable. `--prepare-only` stops after projection and inventory.

Orbit 19 outputs are isolated in `orbit00019_85km_l1f_input`,
`orbit00019_85km_l1f_pilot`, and `orbit00019_85km_l1f_full`.
The `85km` label identifies the reference altitude, not the orbit number.

## What using L1F means here

The native L1F frames are 256 by 256 camera pixels. The existing separation
model assumes translations on the Master grid, so this implementation projects
L1F to 300 by 300 grid windows before fitting. It does not fit directly in
camera coordinates. Neither existing L1R radiance nor its interpolation kernel
is used by the runner.

Projection solves the bilinear latitude/longitude geometry of adjacent camera
pixels, including longitude wrapping. Radiance uses the same four neighbors.
All four must be finite; no extrapolation or interpolation across missing
neighbors is allowed. Source timestamps and temporal gaps are unchanged.
Input and Master file identities are recorded and checked before cache reuse.

This avoids the legacy L1R generator's integer count rounding, but changes its
spatial interpolation rule too. It is not established that this method is
scientifically better. The fitted errors and q scores concern spatially
resampled L1F observations, not untouched detector pixels. Adjacent pixels and
frames remain correlated; q is a selection diagnostic, not independent proof
of cloud identity.

The synthetic adapter tests cover plane reconstruction, curved/rotated
geometry, longitude wrapping, missing neighbors, and footprint boundaries.
Three real orbit 19 frames (0, 960, 1922) have exactly matching L1F/L1R epochs
and correlations of 0.9976–0.9992 on jointly valid pixels. Their radiance RMSE
is 0.040–0.060 kR. These comparisons check alignment and plausibility, not
superiority; see `orbit19_l1f_comparison.json`.
