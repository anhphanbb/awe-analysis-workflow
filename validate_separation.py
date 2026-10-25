"""Small numerical checks for motion direction, missing-data handling and outputs."""
import separate_layers as s
import numpy as np
from scipy.ndimage import gaussian_filter
from PIL import Image

rng=np.random.default_rng(2026)
a=gaussian_filter(rng.normal(size=(s.ny,s.nx)),1)
v=np.array([.5,-.25]);dt=12
b=s.warp(a,v,dt)
cy,cx=150,150
t=a[cy-15:cy+16,cx-15:cx+16]
search=b[cy-25:cy+26,cx-25:cx+26]
c=s.ncc(t,search)
py,px=np.unravel_index(np.nanargmax(c),c.shape)
assert (px-10,py-10)==(6,-3),(px,py)
assert np.nanmax(c)>.99999
assert np.isnan(b[:,-1]).sum()==3
assert np.isnan(b[:,:6]).all()
back=s.warp(b,v,-dt)
good=np.isfinite(back)
assert np.max(np.abs(back[good]-a[good]))<1e-10
d=np.load(s.OUT/'separated_layers.npz')
fit=np.isfinite(d['residual'])
assert fit.sum()>500000
np.testing.assert_allclose((d['observed']-d['slow']-d['fast'])[fit],d['residual'][fit],atol=2e-6)
assert np.array_equal(d['fast_mask'],d['fast_coherence']>.6)
assert d['observed'].shape==(50,s.ny,s.nx)
assert Image.open(s.OUT/'separation.gif').n_frames==50
assert len(list(s.OUT.glob('frame_*.png')))==50
print('PASS: known-translation NCC, warp direction/inversion, no edge wrapping, reconstruction identity, mask consistency, 50 output frames.')
