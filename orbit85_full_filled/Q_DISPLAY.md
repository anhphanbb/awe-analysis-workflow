# Radiance and q display

q_radiance_swaths.png shows L3C radiance, existing MLCloud, slow and fast components, original L1R radiance at scoring timestamps, all-frame q, original L1R radiance masked by q > 0.6, and model availability. Four q_radiance_segment_*.png files provide enlarged views.

The narrow stripes in the withheld-only score mark training timestamps, not physical structures. The new Q_AllFrames_* variables evaluate all timestamps directly, including training frames, so they are descriptive diagnostics rather than independent validation. The original Q_* withheld-only variables remain separate and unchanged. No temporal interpolation or cosmetic gap filling is applied. Single-slow and unsupported regions remain unavailable.

The masked panel shows original L1R radiance on the same scale as the unmasked radiance panels. q is evaluated using the dominant window while the slow and fast orbit maps blend compatible windows.
