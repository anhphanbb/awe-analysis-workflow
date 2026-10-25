# Fixed latitude/longitude demonstration

Source: `Y:\soc\l1r\2023\11\awe_l1r_q20_2023326T0547_00004_v01.nc`.

Open `index.html` for playback, manual stepping, and an endpoint blink. `fixed_latlon.gif` is a portable animation. Frame indices are explicitly zero-based: 930–979, with reference index 955. If your frame numbers are one-based, reference index 955 corresponds to frame 956.

Each frame's calibrated radiance is linearly interpolated from its supplied pixel Longitude/Latitude onto the same 0.02-degree rectangular grid. The map bounds come from frame index 955. Source coordinates use the file's 87 km OH layer altitude. No temporal registration, frame-dependent brightness adjustment, or motion correction is applied. Interpolation is within the convex hull of valid source samples; outside it is missing. This is a quick visual demonstration, not a precision retrieval.

The gold outline bounds grid cells available in every frame. The right panel subtracts a Gaussian-smoothed version of each projected frame (sigma 10 grid cells, or 0.2 degrees). The smoothing is normalized by valid-data weight; residuals near missing-data boundaries are suppressed. This panel enhances spatial contrast without shifting features. Both panels use fixed scales computed across the full sequence.

Apparent movement is relative to the 87 km geolocation reference. This demonstration alone does not distinguish height-dependent parallax from actual scene evolution or geolocation error. The jump when playback wraps is not physical motion.

`projected_frames.npz` contains radiance (kRayleigh), spatial_residual (kRayleigh), latitude, longitude, frame_index, epoch (unchanged source values), and common_coverage. `metadata.json` records scales and parameters. No UTC conversion is applied because the file labels Epoch both as GPS seconds and as seconds since a UTC origin; relative times are unambiguous for this short sequence.

Reproduce with `make_demo.py` in the parent directory using Python with numpy, scipy, matplotlib, netCDF4, and Pillow. The script first checks the workspace-local `.packages` installation.
