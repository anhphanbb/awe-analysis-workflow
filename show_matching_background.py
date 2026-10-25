"""Display the exact sigma=5 background used by motion matching."""
import separate_layers as s
import numpy as np
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt

j=25
a=s.I[j]
ok=np.isfinite(a)
w=gaussian_filter(ok.astype(float),5)
background=gaussian_filter(np.nan_to_num(a),5)/np.maximum(w,1e-9)
background=np.where(ok,background,np.nan)
texture=s.highpass(a)
np.testing.assert_allclose((a-background)[np.isfinite(texture)],texture[np.isfinite(texture)])
np.savez_compressed(s.OUT/'matching_background_0955.npz',original=a,background=background,matching_texture=texture,latitude=s.lat,longitude=s.lon,sigma_grid_pixels=5)
lo,hi=np.nanpercentile(s.I,[2,98])
lim=np.nanpercentile(np.abs(texture),98)
extent=[s.lon[0]-.01,s.lon[-1]+.01,s.lat[0]-.01,s.lat[-1]+.01]
plt.rcParams.update({'font.size':11,'figure.facecolor':'#f5f7fa','axes.facecolor':'#dce0e6'})
fig,axes=plt.subplots(1,3,figsize=(16,6),layout='constrained')
for ax,b,title,cmap,vmin,vmax in zip(axes,[a,background,texture],['Original projected radiance','Smooth spatial background','Texture used for motion matching\n(original − background)'],['gray','gray','RdBu_r'],[lo,lo,-lim],[hi,hi,lim]):
    im=ax.imshow(b,origin='lower',extent=extent,cmap=cmap,vmin=vmin,vmax=vmax,interpolation='nearest',aspect=1/np.cos(np.deg2rad(s.lat.mean())))
    ax.set(title=title,xlabel='Longitude (°E)',ylabel='Latitude (°N)')
    fig.colorbar(im,ax=ax,shrink=.75,label='kRayleigh')
fig.suptitle('Frame 955 (zero-based) · exact motion-matching preprocessing\nGaussian σ = 5 grid pixels = 0.10° (~11 km); original and background share the same brightness scale',fontsize=14)
fig.savefig(s.OUT/'matching_background_0955.png',dpi=130)
plt.close(fig)
p=s.OUT/'index.html'
html=p.read_text(encoding='utf-8')
section='''<section id="matching-background"><h2>Smooth background used for motion matching</h2><img src="matching_background_0955.png" alt="Frame 955 original radiance, Gaussian background, and matching texture"><p>This is the exact preprocessing used to estimate motion: a Gaussian with σ = 5 pixels (0.10°, approximately 11 km), normalized by valid-data weight. The original and background use the same brightness scale. Matching texture is original minus background; pixels with smoothing support ≤98% are excluded. The background is a spatially smoothed observation containing both layers, not the fitted slow airglow component. The earlier projection viewer used a different, broader σ = 10-pixel smooth for its display panel.</p><p><a href="matching_background_0955.png">Full-size comparison</a> · <a href="matching_background_0955.npz">Numerical arrays</a></p></section>'''
if 'id="matching-background"' not in html:
    html=html.replace('<section><h2>Motion measurements</h2>',section+'<section><h2>Motion measurements</h2>')
    p.write_text(html,encoding='utf-8')
print('Saved background comparison; exact subtraction verified.')
