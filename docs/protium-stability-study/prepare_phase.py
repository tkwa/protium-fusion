"""Export conservative radial snapshots for physically evolved timing seeds."""
import argparse,json
from pathlib import Path
import numpy as np
from analyze import snapshot
p=argparse.ArgumentParser();p.add_argument('case',type=Path);p.add_argument('--out',type=Path,required=True);p.add_argument('--backward-case',type=Path);a=p.parse_args()
times=[];states=[];radii=None
for f in sorted(a.case.glob('*.vtk')):
    t,r,theta,d=snapshot([f]);times.append(t)
    if radii is None:radii=r
    assert np.allclose(r,radii)
    states.append(np.stack([d['dens'][0],d['mom'][0,:,0],d['Etot'][0],d['s0'][0],d['s1'][0]],axis=-1))
if a.backward_case:
    backward=[]
    for f in sorted(a.backward_case.glob('*.vtk')):
        t,r,theta,d=snapshot([f]);assert np.allclose(r,radii)
        if t==0:continue
        state=np.stack([d['dens'][0],-d['mom'][0,:,0],d['Etot'][0],d['s0'][0],d['s1'][0]],axis=-1)
        backward.append((-t,state))
    backward.sort(key=lambda x:x[0])
    times=[x[0] for x in backward]+times;states=[x[1] for x in backward]+states
with a.out.open('wb') as f:
    np.array([len(times),len(radii)-1],dtype='<i4').tofile(f)
    np.array(times,dtype='<f8').tofile(f);np.array(radii,dtype='<f8').tofile(f)
    np.array(states,dtype='<f8').tofile(f)
print(json.dumps(dict(snapshots=len(times),radial_zones=len(radii)-1,start_us=times[0]*30,end_us=times[-1]*30)))
