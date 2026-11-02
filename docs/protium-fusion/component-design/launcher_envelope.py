"""Cold uniform-beam alternating-gradient envelope existence screen.
Constant speed periodic slices, paraxial nonrelativistic KV-like uniform ellipse.
No halo, granular instabilities, finite emittance, RF phase acceptance or hardware qualification.
"""
from pathlib import Path
import math,json,hashlib
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root
EPS=8.8541878128e-12;qm=.00265625634384;Imean=429.1650213870357

def propagate(y,k=4.,end=1.):
 def f(s,y,sgn):
  x,z,u,w=y
  if x<=0 or z<=0:raise ValueError('crossed envelope')
  return [u,w,-sgn*k*x+2/(x+z),sgn*k*z+2/(x+z)]
 traces=[]
 for lo,hi,sgn in [(0,.5,1),(.5,1,-1)]:
  if end<=lo:break
  hi=min(hi,end);r=solve_ivp(lambda s,y:f(s,y,sgn),(lo,hi),y,rtol=1e-10,atol=1e-12,dense_output=True,max_step=.01);y=r.y[:,-1];traces.append(r.sol(np.linspace(lo,hi,101)))
 return y,np.hstack(traces)

def residual(y):
 try:return propagate(y)[0]-y
 except ValueError:return np.ones(4)*1e6
r=root(residual,[1.8,1.8,1.8,-1.8]);assert r.success
end,trace=propagate(r.x);err=max(abs(end-r.x));checks=[err<1e-8]
h=1e-5;mon=np.column_stack([(propagate(r.x+np.eye(4)[j]*h)[0]-propagate(r.x-np.eye(4)[j]*h)[0])/(2*h) for j in range(4)]);eigs=np.linalg.eigvals(mon);checks.append(max(abs(eigs))<1.00001)
rows=[]
for speed in [1e4,1e5,5e5,2.4e6]:
 for current_factor in [1,8]:
  I=Imean*current_factor;G=1e4;kappa=qm*G/speed**2;period=math.sqrt(4/kappa);K=qm*I/(2*math.pi*EPS*speed**3);scale=period*math.sqrt(K);amax=float(max(trace[:2].ravel())*scale);amin=float(min(trace[:2].ravel())*scale)
  rows.append(dict(speed_m_s=speed,current_A=I,current_factor=current_factor,quadrupole_gradient_V_m2=G,full_FODO_period_m=period,matched_envelope_min_m=amin,matched_envelope_max_m=amax,quadrupole_field_at_envelope_V_m=G*amax,self_edge_field_round_at_max_envelope_V_m=I/(2*math.pi*EPS*speed*amax),minimum_geometric_bore_diameter_m=2*amax,scope='No wall clearance or halo margin. Factor8 illustrates4x temporal masscurrent peak and50% RF bucket filling, not a derived waveform.'))
  checks += [amax>amin>0,math.isclose(kappa*period**2,4,rel_tol=1e-12)]
out={'scope':__doc__,'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'dimensionless_initial_envelope':r.x.tolist(),'dimensionless_peak_radius':float(max(trace[:2].ravel())),'periodicity_residual':float(err),'envelope_monodromy_eigenvalues_real_imag':[[float(z.real),float(z.imag)] for z in eigs],'max_envelope_eigenvalue_modulus':float(max(abs(eigs))),'rows':rows,'checks':{'count':len(checks),'passed':int(sum(checks)),'all_passed':bool(all(checks))}}
Path(__file__).with_name('launcher-envelope.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out['checks']));print(r.x);print(rows[-2:])
