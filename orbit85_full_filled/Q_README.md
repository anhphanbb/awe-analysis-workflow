# Local fast-component score

{
  "settings": {
    "threshold": 0.6,
    "gaussian_sigma_native_cells": 3.0,
    "min_local_coverage": 0.95,
    "min_weighted_residual_energy_kR2": 0.04,
    "min_training_support": 8,
    "strict_threshold": true,
    "definition": "clip(1 - G[(I-S-F)^2 valid] / G[(I-S)^2 valid],0,1)",
    "evaluation": "Original L1R observations only, exact L3C Pixel_Epoch, odd local indices of dominant accepted two-motion window. No temporal interpolation. This validates templates within windows, not globally independent predictions.",
    "gauge": "Uses offset-aligned dominant window templates. The denominator and fast sign depend on the brightness-split gauge; q is not gauge-invariant. The total reconstruction error is invariant.",
    "interpretation": "Candidate fast-feature mask, not cloud probability. Both signs of fast contrast can qualify. Single-slow and unevaluable pixels remain unavailable.",
    "defaults": "Exploratory, not calibrated. Sigma=3 native cells is nominally 6 km, not identical to three pixels on the old latitude/longitude grid."
  },
  "two_motion_windows": 36,
  "eligible_pixels": 529863,
  "candidate_pixels": 148660,
  "status_counts": {
    "0": 2101626,
    "1": 1101765,
    "2": 540211,
    "3": 24853,
    "4": 33,
    "5": 529863
  }
}

Unavailable scores are NaN and masks are -1, not zero. The source timestamp must be held out from the dominant model. Narrow missing stripes are therefore expected. Energies use identical Gaussian weights on genuine observations with valid predictions and >=8 training observations per component. Values outside the image count as missing, unlike reflected-edge filtering. The >0.6 rule is strict. q compares the final slow component with that same slow component plus fast; this differs from comparing two separately fitted models.

Run `python orbit_full.py score --fill-small-gaps --q-threshold 0.6` to refresh without refitting. Settings: --q-sigma, --q-min-coverage, --q-min-energy. Assembly also refreshes q automatically.
