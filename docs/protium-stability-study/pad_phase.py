"""Extend a reference table with its uniform, stationary ambient exterior."""
import argparse
from pathlib import Path
import numpy as np

p=argparse.ArgumentParser()
p.add_argument('source',type=Path)
p.add_argument('destination',type=Path)
p.add_argument('--factor',type=int,default=2)
a=p.parse_args()
with a.source.open('rb') as f:
    nt,nr=np.fromfile(f,dtype='<i4',count=2)
    times=np.fromfile(f,dtype='<f8',count=nt)
    edges=np.fromfile(f,dtype='<f8',count=nr+1)
    states=np.fromfile(f,dtype='<f8').reshape(nt,nr,5)
assert a.factor>=1 and np.allclose(np.diff(edges),np.diff(edges)[0])
assert np.allclose(states[:,-8:,:],states[:,-1:, :],rtol=1e-10,atol=1e-23)
assert np.max(np.abs(states[:,-1,1]/states[:,-1,0]))<1e-14
assert np.all(states[:,-1,3:]==0)
newnr=nr*a.factor
newstates=np.concatenate([states,np.repeat(states[:,-1:,:],newnr-nr,axis=1)],axis=1)
with a.destination.open('wb') as f:
    np.array([nt,newnr],dtype='<i4').tofile(f)
    times.tofile(f)
    np.linspace(edges[0],edges[-1]*a.factor,newnr+1,dtype='<f8').tofile(f)
    newstates.tofile(f)
print(a.destination,newstates.shape)
