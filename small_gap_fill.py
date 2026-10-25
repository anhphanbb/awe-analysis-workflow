"""Conservative filling of small enclosed missing-radiance components."""
import numpy as np
from scipy.ndimage import label,find_objects

def fill_small(a,geometry):
    out=a.copy();filled=np.zeros(a.shape,bool);labels,n=label(~np.isfinite(a),structure=np.ones((3,3)));objects=find_objects(labels);sizes=np.bincount(labels.ravel());accepted=[]
    for k,sl in enumerate(objects,1):
        if sl is None or sizes[k]>16:continue
        sy,sx=sl
        if sy.start==0 or sx.start==0 or sy.stop==a.shape[0] or sx.stop==a.shape[1]:continue
        if sy.stop-sy.start>5 or sx.stop-sx.start>5:continue
        yy,xx=np.where(labels[sl]==k);yy+=sy.start;xx+=sx.start
        if not geometry[yy,xx].all():continue
        coords=list(zip(yy,xx));idx={pt:i for i,pt in enumerate(coords)};A=np.eye(len(coords))*4;b=np.zeros(len(coords));boundary=[]
        for i,(y,x) in enumerate(coords):
            for pt in [(y-1,x),(y+1,x),(y,x-1),(y,x+1)]:
                if pt in idx:A[i,idx[pt]]-=1
                else:
                    assert np.isfinite(a[pt]);b[i]+=a[pt];boundary.append(float(a[pt]))
        values=np.linalg.solve(A,b);assert values.min()>=min(boundary)-1e-5 and values.max()<=max(boundary)+1e-5
        out[yy,xx]=values;filled[yy,xx]=True;accepted.append(int(sizes[k]))
    return out,filled,accepted
