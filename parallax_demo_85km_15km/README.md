# AWE 85 km versus 15 km demo

Open index.html in a browser. Images are embedded, so the interactive animation works from this single HTML file. Other download links use adjacent files in this ZIP.

The original demo epochs are matched exactly to v02. Both mappings use identical radiance and timestamps. See metadata.json for the skipped original index and source details. Small enclosed holes are filled before projection, with masks retained in the workspace file source_interpolation_provenance.npz. Projected data and spatial contrasts are retained in workspace file projected_frames.npz. Numerical arrays are omitted from this portable demo ZIP. Larger missing regions remain missing.

Columns compare mapping altitude, not separated physical layers. All panels use the same fixed geographic grid; scales are shared within each row. The common-coverage outline is computed independently for each altitude. Numerical source arrays retain the original disk orientation. Reproduce using make_altitude_demo.py in the parent workspace.
