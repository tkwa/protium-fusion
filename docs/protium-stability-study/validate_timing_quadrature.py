"""Audit actual mass-weighted timing errors in Cartesian initialization quadrature."""
import argparse,json,sys
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

base=Path('/var/tmp/tkwa-chatgpt/protium-stability')
sys.path.insert(0,str(base/'athena/vis/python'));import athena_read
p=argparse.ArgumentParser();p.add_argument('case',type=Path);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
cfg=json.loads((a.case/'config.json').read_text())
with Path(cfg['phase_table']).open('rb') as f:
    nt,nr=np.fromfile(f,'<i4',2)
    times=np.fromfile(f,'<f8',nt);radii=np.fromfile(f,'<f8',nr+1)
    states=np.fromfile(f,'<f8').reshape(nt,nr,5)
errors=np.loadtxt(cfg['packet_times']);n=len(errors)
z=(np.arange(n//2)+.5)/(n//2);phi=np.arange(n//2)*np.pi*(3-np.sqrt(5))
v=np.column_stack([np.sqrt(1-z*z)*np.cos(phi),np.sqrt(1-z*z)*np.sin(phi),z])
tree=cKDTree(np.stack([v,-v],axis=1).reshape(-1,3))
sigmas=[0.,3.,4.5,6.]
totals=np.zeros((len(sigmas),4))
nodes=[.5-.5/np.sqrt(3),.5+.5/np.sqrt(3)]
fuel_radius=radii[1:][np.any(states[:,:,4]>0,axis=0)].max()
for file in sorted(a.case.glob('*.out2.00000.vtk')):
    x,y,z,_=athena_read.vtk(file)
    x,y,z=[q.astype(float) for q in [x,y,z]]
    nearest=np.array([max(0,float(q[0]),-float(q[-1])) for q in [x,y,z]])
    if np.linalg.norm(nearest)>fuel_radius:continue
    dv=float((x[1]-x[0])*(y[1]-y[0])*(z[1]-z[0]))/8
    for aa in nodes:
      for bb in nodes:
       for cc in nodes:
        xx,yy,zz=np.meshgrid(x[:-1]+aa*np.diff(x),y[:-1]+bb*np.diff(y),z[:-1]+cc*np.diff(z),indexing='ij')
        xyz=np.column_stack([q.ravel() for q in [xx,yy,zz]]).astype(float)
        radius=np.linalg.norm(xyz,axis=1)
        _,index=tree.query(xyz/radius[:,None],workers=1)
        error=errors[index]
        ir=np.clip(np.searchsorted(radii,radius,side='right')-1,0,nr-1)
        for k,sigma in enumerate(sigmas):
            t=cfg['reference_time_us']/30+sigma*error/30
            it=np.clip(np.searchsorted(times,t,side='right')-1,0,nt-2)
            w=(t-times[it])/(times[it+1]-times[it])
            fuel=(1-w)*states[it,ir,4]+w*states[it+1,ir,4]
            energy=(1-w)*states[it,ir,2]+w*states[it+1,ir,2]
            mass=fuel*dv
            totals[k]+=[mass.sum(),mass@(sigma*error),mass@((sigma*error)**2),energy.sum()*dv]
rows=[]
for sigma,(mass,m1,m2,energy) in zip(sigmas,totals):
    rows.append(dict(nominal_sigma_us=sigma,mean_error_us=m1/mass,
                     actual_sigma_us=float(np.sqrt(max(0,m2/mass-(m1/mass)**2))),
                     initialized_fuel_mass_kg=mass*2.7e12,
                     mass_ratio_to_zero_jitter=mass/totals[0,0]))
a.out.write_text(json.dumps(dict(mesh_case=a.case.name,rows=rows),indent=2)+'\n')
print(json.dumps(rows,indent=2))
