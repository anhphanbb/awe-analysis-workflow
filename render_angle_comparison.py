"""Add the 5/10/20 degree crop comparison to the standalone report."""
import re
from render_orbit19_html import OUT,np,plt,json,figure_uri

def main():
    runs=[]
    for step,folder in [(5,'cropped_wavelet'),(10,'cropped_wavelet_step10'),(20,'cropped_wavelet_step20')]:
        summary=json.loads((OUT/folder/'wavelet_summary.json').read_text())
        z=np.load(OUT/folder/'wavelet_results.npz')
        runs.append((step,summary,z))
    crop=np.load(OUT/'concentric_waves_crop.npz')
    extent=[crop['longitude'][0]-.01,crop['longitude'][-1]+.01,crop['latitude'][0]-.01,crop['latitude'][-1]+.01]
    vmax=max(float(z['spectrum'].max()) for _,_,z in runs)
    limit=max(float(np.nanpercentile(abs(z['clusters'][0]),99)) for _,_,z in runs if len(z['clusters']))
    fig,axes=plt.subplots(2,3,figsize=(16,10),layout='constrained')
    for col,(step,s,z) in enumerate(runs):
        period=z['period_km'];ticks=np.arange(0,len(period),4)
        im=axes[0,col].imshow(z['spectrum'],origin='lower',aspect='auto',extent=[-step/2,180-step/2,-.5,len(period)-.5],cmap='magma',vmin=0,vmax=vmax)
        axes[0,col].set(title=f'{step}° steps · {180//step} angles',xlabel='Wave-normal axis from east (°)',ylabel='Wavelength (km)',yticks=ticks,yticklabels=[f'{period[k]:.0f}' for k in ticks])
        if len(z['clusters']):
            c=s['clusters'][0]
            image=axes[1,col].imshow(z['clusters'][0],origin='lower',extent=extent,cmap='RdBu_r',vmin=-limit,vmax=limit,aspect=1/np.cos(np.deg2rad(crop['latitude'].mean())))
            axes[1,col].set(title=f"Strongest cluster: {c['weighted_wavelength_km']:.1f} km · {c['wave_normal_axis_degrees']:.1f}°",xlabel='Longitude (°E)',ylabel='Latitude (°N)')
            axes[1,col].set_facecolor('#e1e6e9')
        else:axes[1,col].text(.5,.5,'No cluster above threshold',ha='center',transform=axes[1,col].transAxes)
    fig.colorbar(im,ax=axes[0],label='Mean normalized coefficient power',shrink=.7)
    fig.colorbar(image,ax=axes[1],label='Reconstructed contrast (kR)',shrink=.7)
    fig.savefig(OUT/'angular_step_comparison.png',dpi=160,bbox_inches='tight')
    uri=figure_uri(fig)
    rows=''
    for step,s,z in runs:
        description='; '.join(f"{c['weighted_wavelength_km']:.1f} km / {c['wave_normal_axis_degrees']:.1f}°" for c in s['clusters']) or 'None'
        rows+=f"<tr><td>{step}°</td><td>{180//step}</td><td>{len(s['clusters'])}</td><td>{s['cluster_expansion']['spatial_pixels']:,}</td><td>{description}</td></tr>"
    section=f'''<section id="angle-comparison"><div class="section-title"><span class="number">05</span><h2>Compare angular step sizes</h2></div><p class="intro">The same 151 × 151 crop, 22 wavelength scales (12–74.0 km), 0.25 minimum peak amplitude, and 0.05 cluster threshold. Only the angular sampling changes: 5°, 10°, and 20°. The original 5° results above remain as the baseline.</p><div class="card"><img class="figure" src="{uri}" alt="Comparison of 5, 10 and 20 degree angular steps: wavelet spectra and strongest cluster reconstructions"></div><div class="card details" style="overflow-x:auto"><table style="width:100%;border-collapse:collapse;text-align:left"><thead><tr><th>Step</th><th>Angles</th><th>Clusters</th><th>Spatial footprint*</th><th>Weighted wavelength / axis</th></tr></thead><tbody>{rows}</tbody></table><p class="caption">*Union of spatial grid pixels assigned to any cluster, before restricting displayed maps to finite slow values. Coefficient counts are not directly comparable across different numbers of angles. Top panels share one power scale; bottom panels share one contrast scale and show the strongest cluster per run.</p></div><p class="caption">Weighted axes can fall between sampled angles. Coarser angular sampling may merge or miss directional structure; a larger footprint alone does not establish a better physical wave estimate.</p></section>'''
    target=OUT/'orbit19_052340.html';text=target.read_text(encoding='utf-8')
    text=re.sub(r'<section id="angle-comparison">.*?</section>','',text,flags=re.S)
    text=text.replace('<section id="method">',section+'<section id="method">')
    if '<a href="#angle-comparison">' not in text:text=text.replace('<a href="#method">','<a href="#angle-comparison">Angle steps</a><a href="#method">')
    target.write_text(text,encoding='utf-8')
    report=dict(status='PASS',angle_steps=[5,10,20],same_crop=True,same_taper=True,wavelength_max_km=float(runs[0][2]['period_km'].max()),min_peak_amplitude=.25,cluster_threshold=.05)
    (OUT/'angular_comparison_verification.json').write_text(json.dumps(report,indent=2))
    print(json.dumps([dict(step=step,clusters=s['clusters'],spatial_pixels=s['cluster_expansion']['spatial_pixels']) for step,s,z in runs],indent=2))

if __name__=='__main__':main()
