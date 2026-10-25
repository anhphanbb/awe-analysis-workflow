"""Numerical checks for the L1F projection, independent of orbit fitting."""
import l1f_adapter as a
from l1f_adapter import np

def check():
    lat,lon=np.mgrid[0:10,0:12].astype(float)
    radiance=2*lat+3*lon
    ta,tb=np.mgrid[.2:8.5:.7,.3:10.5:.8]
    out,geometry=a.resample(lat,lon,radiance,ta,tb)
    assert geometry.all()
    np.testing.assert_allclose(out,2*ta+3*tb,atol=1e-5)
    radiance[3,3]=np.nan
    out,geometry=a.resample(lat,lon,radiance,np.array([[3.2,15.]]),np.array([[3.3,5.]]))
    assert np.isnan(out).all() and geometry[0,0] and not geometry[0,1]
    wrapped=(lon+175+180)%360-180
    out,_=a.resample(lat,wrapped,2*lat+3*lon,ta,(tb+175+180)%360-180)
    np.testing.assert_allclose(out,2*ta+3*tb,atol=1e-5)
    # Curved, rotated camera coordinates: generate targets from known bilinear cells.
    lat2=20+.02*lat+.004*lon+.0001*lat*lon
    lon2=150-.003*lat+.021*lon+.0002*lat*lon
    ta2=20+.02*ta+.004*tb+.0001*ta*tb
    tb2=150-.003*ta+.021*tb+.0002*ta*tb
    out,geometry=a.resample(lat2,lon2,2*lat+3*lon,ta2,tb2)
    assert geometry.all()
    np.testing.assert_allclose(out,2*ta+3*tb,atol=1e-5)
    return ['Planar radiance reconstruction and orientation','Missing-corner exclusion and no extrapolation','Dateline wrapping','Inverse mapping on curved rotated camera geometry']

if __name__=='__main__':print(check())
