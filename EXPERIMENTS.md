# Experiment progression

| Stage | Scripts / notebooks | Saved results | Status and purpose |
|---|---|---|---|
| Geographic reprojection and parallax | `make_demo.py`, `make_altitude_demo.py` | `parallax_demo`, `parallax_demo_85km_15km` | Geographic-grid visualization and altitude comparison |
| Initial layer separation | `separate_layers.py`, `make_walkthrough.py`, `validate_separation.py` | `layer_separation` | Slow/fast decomposition exploration |
| Synthetic experiments | `synthetic_layer_experiment.py`, `synthetic_layer_experiment.ipynb`, `synthetic_separation.py`, `synthetic_separation_64.py` | `synthetic_demo`, `synthetic_demo_64` | Known-input demonstrations and iteration comparisons |
| Actual-data notebook packages | `package_actual_data.py`, `package_v02.py`, `finish_v02_package.py` | `AWE_actual_data`, `AWE_v02_actual_data` | Actual data and epoch-matched v02 notebook workflows |
| Pilot and whole-orbit processing | `orbit_pilot.py`, `orbit_full.py`, `plot_orbit_pilot.py` | `orbit85_pilot`, `orbit85_full` | Earlier orbit pilot and windowed orbit fit |
| Missing pixels and gap filling | `diagnose_orbit_windows.py`, `report_orbit_diagnostics.py`, `small_gap_fill.py`, animation scripts | `orbit85_diagnostics`, `orbit85_full_filled` | Diagnose coverage and test limited filling |
| Shifted windows | `combine_orbit_passes.py`, `verify_shift_comparison.py` | `orbit85_full_filled_shift12`, `orbit85_shift_comparison` | Compare shifted window support |
| Earlier wave analysis | `whole_frame_wave_clusters.py`, `slow_wave_phase_speed.py`, `wave_cluster_limited.py` | `whole_frame_clusters_0955`, `whole_frame_clusters_1180`, `slow_phase_frame0955` | Whole-frame clustering and phase-speed exploration |
| Orbit 19 L1F trial | `l1f_adapter.py`, `verify_l1f_adapter.py`, `L1F_ORBIT_RUN.md` | `orbit00019_85km_l1f_input`, `orbit00019_85km_l1f_full` | Adapter tested; full processing interrupted when switching to L1R |
| Orbit 19 L1R full orbit | `run_orbit.py`, generalized orbit scripts | `orbit00019_85km_l1r_*` | 1923 frames, 76 windows, 32 iterations; verification summaries saved |
| Orbit 19 near 05:23:40 UTC | `build_orbit19_actual.py`, `render_orbit19_html.py` | `AWE_v02_orbit19_052340` | 50-frame local decomposition; initial 64 iterations |
| Footprint and crop changes | `run_cropped_wavelet.py` | `cropped_wavelet` and historical wavelet summaries under the orbit 19 report directory | Whole finite slow component replaced the extra observed/support gate; then wavelet input restricted to the upper concentric-wave crop |
| Angular step comparison | `compare_crop_angles.py`, `render_angle_comparison.py` | `cropped_wavelet_step10`, `cropped_wavelet_step20` | 5/10/20-degree sampling at threshold 0.05, wavelengths below 80 km |
| Latest continuation and grid | `run_orbit19_128_grid.py`, `render_grid128.py` | `grid128`, `grid128_summary.json`, `grid128_verification.json` | 128 iterations; threshold 0.02/0.05 by 5/10/20 degrees; second spectrum line plot removed |

Paths in the last four rows are relative to `AWE_v02_orbit19_052340` where appropriate. This table describes the experiment progression, not Git commit history: the experiments predate this repository snapshot. Historical 64-iteration HTML and fit summaries remain alongside the latest 128-iteration report.
