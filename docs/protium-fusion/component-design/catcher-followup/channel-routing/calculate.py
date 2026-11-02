#!/usr/bin/env python3
"""Finite-volume first-passage routing in an ideal anisotropic graphite strip.

Coordinates: z normal inward, y lateral. A half-cell has a reflecting midplane,
a reflecting rear, an absorbing cavity face, and an absorbing side collection
boundary buried beneath a continuous cap. This resolves the cap + channel domain
jointly. It predicts routing of an already-mobile particle only. It contains no
retained fraction, pore backpressure, straggling distribution, damage or seal law.
"""
from pathlib import Path
import hashlib,json,math
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
P=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()

def solve(pitch_um,cap_um,depth_um,anisotropy,nz=160,ny=64):
    # Finite-volume cells, graded in depth. Cap transition is exactly a grid face.
    if cap_um>0:
        edges=np.r_[np.linspace(0,cap_um,max(8,nz//10)+1),np.geomspace(cap_um,depth_um,nz+1)[1:]]
    else:
        edges=np.r_[0.,np.geomspace(1e-4,depth_um,nz+1)]
    z=(edges[1:]+edges[:-1])/2;dz=np.diff(edges);N=len(z);dy=pitch_um/2/ny
    # Units um: choose Dnormal=1um²/s, Dlateral=anisotropy um²/s.
    rows=[];cols=[];values=[];b=np.zeros(N*ny);vol=np.repeat(dz*dy,ny)
    for i in range(N):
      for j in range(ny):
        k=i*ny+j;diag=0.
        if i>0:
            a=dy/(z[i]-z[i-1]);diag+=a;rows.append(k);cols.append(k-ny);values.append(-a)
        else: diag+=dy/z[0] # absorbing front, value zero
        if i<N-1:
            a=dy/(z[i+1]-z[i]);diag+=a;rows.append(k);cols.append(k+ny);values.append(-a)
        # bottom reflects
        if j>0:
            a=anisotropy*dz[i]/dy;diag+=a;rows.append(k);cols.append(k-1);values.append(-a)
        if j<ny-1:
            a=anisotropy*dz[i]/dy;diag+=a;rows.append(k);cols.append(k+1);values.append(-a)
        elif z[i]>=cap_um:
            a=anisotropy*dz[i]/(dy/2);diag+=a;b[k]+=a # side collection value one
        rows.append(k);cols.append(k);values.append(diag)
    A=coo_matrix((values,(rows,cols)),shape=(N*ny,N*ny)).tocsr()
    u=spsolve(A,b).reshape(N,ny)
    # Expected first-exit time, -div(D grad t)=1. Only a scaling reference;
    # trapping and temperature-dependent coefficients require a time-domain model.
    t=spsolve(A,vol).reshape(N,ny)
    assert u.min()>-1e-9 and u.max()<1+1e-9 and t.min()>=0
    return z,u.mean(axis=1),t.mean(axis=1)

def analytic_front_loss(z,pitch,eta):
    # Semi-infinite strip with channels reaching the front; uniform lateral source.
    ns=np.arange(1,4000,2,dtype=float)
    return np.exp(-np.outer(np.asarray(z)*math.sqrt(eta)*math.pi/pitch,ns))@(8/math.pi**2/ns**2)

def main():
    source=P.parent/'stopping/surface-stopping.json';d=json.loads(source.read_text())
    rows=[];checks=[]
    for pitch in [5.,10.,20.,50.]:
      for cap in [.01,.1,1.]:
       for eta in [1.,10.,100.,1000.,10000.,100000.]:
        z,u,t=solve(pitch,cap,150.,eta)
        stages={}
        for stage,sd in d['stages'].items():
          species={}
          for sp,x in sd['species'].items():
            depths=np.array(x['packet_mean_projected_depth_um']);weights=np.array(x['packet_nuclei_m2']);weights/=weights.sum()
            # u(0)=0. Reflecting rear interpolation; depths beyond150um are
            # separately reported instead of silently counted as collected.
            v=np.interp(depths,np.r_[0,z,150.],np.r_[0,u,u[-1]])
            tt=np.interp(depths,np.r_[0,z,150.],np.r_[0,t,t[-1]])
            species[sp]={'mobile_front_exit_fraction':float(weights@(1-v)),
              'mobile_collection_fraction':float(weights@v),
              'expected_exit_s_if_Dnormal_1um2_s':float(weights@tt),
              'number_fraction_depth_beyond_domain':float(weights[depths>150.].sum())}
          stages[stage]=species
        rows.append({'pitch_um':pitch,'buried_cap_um':cap,'active_depth_um':150.,'Dparallel_over_Dnormal':eta,
          'grid_depth_cells':len(z),'grid_half_pitch_cells':64,'stages':stages})
    # Validation against a distinct closed-form infinite-strip solution and mesh
    # refinement at the difficult shallow-source/high-anisotropy operating point.
    z,u,t=solve(20.,0.,150.,1000.,320,128)
    select=(z>.01)&(z<1.)
    analytic=analytic_front_loss(z[select],20.,1000.)
    err=float(np.max(abs((1-u[select])-analytic)))
    checks.append({'name':'uncapped finite-volume vs infinite-strip analytic','max_absolute_error':err,'pass':err<.003})
    refine=[]
    for cap in [.01,.1,1.]:
      for eta in [1000.,100000.]:
        z,u,t=solve(20.,cap,150.,eta,320,128)
        row=next(x for x in rows if x['pitch_um']==20. and x['buried_cap_um']==cap and x['Dparallel_over_Dnormal']==eta)
        for stage,sd in d['stages'].items():
          for sp,x in sd['species'].items():
            depths=np.array(x['packet_mean_projected_depth_um']);w=np.array(x['packet_nuclei_m2']);w/=w.sum()
            fine=float(w@(1-np.interp(depths,np.r_[0,z,150.],np.r_[0,u,u[-1]])))
            coarse=row['stages'][stage][sp]['mobile_front_exit_fraction']
            refine.append({'stage':stage,'species':sp,'cap_um':cap,'eta':eta,'coarse':coarse,'fine':fine,'absolute_difference':abs(fine-coarse)})
    maxerr=max(x['absolute_difference'] for x in refine)
    checks.append({'name':'resolution 160x64 vs320x128','max_absolute_routing_difference':maxerr,'pass':maxerr<.003})
    checks.append({'name':'mobile routing partition','pass':all(abs(v['mobile_front_exit_fraction']+v['mobile_collection_fraction']-1)<1e-12 for row in rows for st in row['stages'].values() for v in st.values())})
    out={'scope':__doc__,'source_sha256':sha(source),'code_sha256':sha(__file__),
      'numerics':'Cell-centered conservative finite-volume elliptic solve, lateral symmetry, finite reflecting rear. Channels have idealzero width; finite roof width/directly irradiated channels require a separate design. Per-energy mean projected depths ignore straggling and reflection. Spatially constant diffusion tensor. No measured reactor retention or engineered capture is inferred.',
      'rows':rows,'refinement':refine,'checks':checks}
    (P/'results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps(checks,indent=2));assert all(x['pass'] for x in checks)
if __name__=='__main__':main()
