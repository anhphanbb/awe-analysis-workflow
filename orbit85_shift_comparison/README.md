# Shifted-window comparison

Pass A uses starts 0,25,50,...; B uses 12,37,62,... plus the original endpoint windows. Endpoints shared by both passes do not provide distinct-window confirmation. Same_Dominant_Window flags them, and the final masked-radiance panel excludes them.

Categories distinguish both-above-threshold, only-one-eligible candidates, disagreement when both are eligible, below-threshold evidence, and no evidence. No maximum-q map or undifferentiated union is presented. Ineligible is not cloud-free. A candidate is not a calibrated cloud detection.

Withheld maps only use scores from models that excluded that source timestamp. All-frame maps include training observations and show descriptive fit sensitivity. Overlapping windows share data, so agreement is not independent replication. q is computed in each pass's aligned brightness gauge and is sensitive to the slow/fast brightness split. The paired RMSE comparison uses exactly the same genuine observations held out from both models; its total prediction is gauge-invariant. This comparison is not an untouched test because withheld errors also select models.

Radiance masked in comparison plots is the existing coadded L3C, used as geographic context; q itself uses single-frame L1R at Pixel_Epoch. The two q fields and eligibility provenance remain separate in orbit85_shift_comparison.nc.

Reproduce:
```
python orbit_full.py fit --fill-small-gaps --start-offset 12 --workers 6 --iterations 32
python orbit_full.py assemble --fill-small-gaps --start-offset 12
python combine_orbit_passes.py
```
