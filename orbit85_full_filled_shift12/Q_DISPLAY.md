# Radiance and q display

q_radiance_swaths.png includes L3C radiance, existing MLCloud, slow and fast components, original L1R radiance at scoring timestamps, all-frame q, original L1R radiance masked by q > 0.6, and model availability. Four q_radiance_segment_*.png files provide enlarged views.

Narrow stripes in the withheld-only q_score_whole_orbit.png represent alternate timestamps used for training. The new display evaluates every timestamp directly, including training frames; it is a descriptive diagnostic, not independent validation. No gaps are cosmetically filled or scores interpolated. Q_AllFrames_* variables are separate from the original withheld-only Q_* variables in orbit85_separation.nc. Remaining gaps reflect single-slow models, insufficient support or weak signal.

The mask shows actual original L1R radiance, not a binary color block or the coadded L3C values. It uses the same radiance scale as the unmasked radiance panels. q uses the dominant two-motion window; slow and fast orbit panels are blended across compatible windows.
