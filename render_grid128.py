"""Render the six-setting crop-only wavelet comparison after 128 fit iterations."""
import re
from render_orbit19_html import OUT,np,plt,json,figure_uri

def main():
    report=json.loads((OUT/'grid128_summary.json').read_text());runs=report['runs']
    loaded={(r['angle_step'],r['threshold']):(r,np.load(OUT/'grid128'/f"step{r['angle_step']:02d}_thr{int(r['threshold']*1000):03d}"/'wavelet_results.npz')) for r in runs}
    crop=np.load(OUT/'concentric_waves_crop.npz')
    extent=[crop['longitude'][0]-.01,crop['longitude'][-1]+.01,crop['latitude'][0]-.01,crop['latitude'][-1]+.01]
    limit=max(float(np.nanpercentile(abs(z['clusters'][0]),99)) for r,z in loaded.values() if len(z['clusters']))
    fig,axes=plt.subplots(2,3,figsize=(15,11),layout='constrained')
    rows=''
    for iy,threshold in enumerate([.02,.05]):
        for ix,step in enumerate([5,10,20]):
            r,z=loaded[(step,threshold)];ax=axes[iy,ix]
            assert np.all(z['period_km']<80);np.testing.assert_allclose(z['slow'],crop['slow'],equal_nan=True,atol=1e-6)
            if len(z['clusters']):
                c=r['clusters'][0]
                im=ax.imshow(z['clusters'][0],origin='lower',extent=extent,cmap='RdBu_r',vmin=-limit,vmax=limit,aspect=1/np.cos(np.deg2rad(crop['latitude'].mean())))
                ax.set(title=f"{step}° / threshold {threshold:.2f}\n{c['weighted_wavelength_km']:.1f} km · {c['wave_normal_axis_degrees']:.1f}°",xlabel='Longitude (°E)',ylabel='Latitude (°N)')
            else:ax.text(.5,.5,'No cluster above peak threshold',transform=ax.transAxes,ha='center')
            ax.set_facecolor('#e1e6e9')
            description='; '.join(f"{c['weighted_wavelength_km']:.1f} km / {c['wave_normal_axis_degrees']:.1f}°" for c in r['clusters']) or 'None'
            rows+=f"<tr><td>{step}°</td><td>{threshold:.2f}</td><td>{len(r['clusters'])}</td><td>{r['spatial_pixels']:,}</td><td>{description}</td></tr>"
    fig.colorbar(im,ax=axes.ravel().tolist(),label='Reconstructed contrast (kR)',shrink=.65)
    fig.savefig(OUT/'grid128_comparison.png',dpi=155,bbox_inches='tight')
    uri=figure_uri(fig)
    section=f'''<section id="grid-search"><div class="section-title"><span class="number">05</span><h2>Threshold × angular-step grid search</h2></div><p class="intro">Six runs on the same cropped slow component after 128 separation iterations. Rows compare inclusion thresholds 0.02 and 0.05; columns compare angular steps 5°, 10°, and 20°. All runs use 22 wavelength scales from 12 to 74.0 km, peak amplitude threshold 0.25, and the same crop, plane removal, normalization and 40 km taper. The default above is 5° / 0.02.</p><div class="card"><img class="figure" src="{uri}" alt="Six-setting grid search: threshold 0.02 and 0.05 by angular steps 5, 10, and 20 degrees"></div><div class="card details" style="overflow-x:auto"><table style="width:100%;text-align:left;border-collapse:collapse"><thead><tr><th>Angle step</th><th>Inclusion threshold</th><th>Clusters</th><th>Spatial pixels*</th><th>Weighted wavelength / axis</th></tr></thead><tbody>{rows}</tbody></table><p class="caption">*Union of pixels assigned to any cluster, before display masking. The figure shows the strongest cluster per run on one shared contrast scale. A larger footprint can include weaker features or merge patterns; it is not an independent accuracy score.</p></div><p class="caption">This is a sensitivity grid, not a validated optimization. The 0.02 threshold is the requested default. The reconstruction-RMSE plot is retained; the additional power-spectrum line plot has been removed.</p></section>'''
    target=OUT/'orbit19_052340.html';text=target.read_text(encoding='utf-8')
    text=re.sub(r'<section id="grid-search">.*?</section>','',text,flags=re.S)
    text=text.replace('<section id="method"><div class="section-title"><span class="number">05</span>', '<section id="method"><div class="section-title"><span class="number">06</span>')
    text=text.replace('<section id="method">',section+'<section id="method">')
    if '<a href="#grid-search">' not in text:text=text.replace('<a href="#method">','<a href="#grid-search">Grid search</a><a href="#method">')
    target.write_text(text,encoding='utf-8')
    verification=dict(status='PASS',iterations=128,grid_runs=6,thresholds=[.02,.05],angle_steps=[5,10,20],same_crop_all_runs=True,wavelengths_strictly_below_80=True,power_spectrum_line_removed=True)
    (OUT/'grid128_verification.json').write_text(json.dumps(verification,indent=2))
    print(json.dumps(verification),flush=True)

if __name__=='__main__':main()
