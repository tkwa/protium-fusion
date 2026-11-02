"""Affine low-angular-mode support screen, with exact quadratic energies.
Uniform thin homoeoid under affine deformation F; guided ellipse circulation
conserves action, and the receiver has a declared isotropic tensile membrane.
This is a9DOF global-mode screen, not a finite-sector stability simulation.
"""
from pathlib import Path
import math,json,hashlib
import numpy as np
P=Path(__file__).parent;inp=P/'sector-account.json';j=json.loads(inp.read_text());R=j['reference']['radius_m'];sigma=j['support']['total_areal_mass_kg_m2'];a=j['support']['self_acceleration_m_s2']/R;n2=j['support']['maximum_tidal_acceleration_m_s2']/(3*R)
b=np.array([(a+4*n2)/2,(a-2*n2)/2,(a-4*n2)/2]);t=np.array([3*n2,0.,-n2]);n=math.sqrt(n2)
def energy2(E,E2D,nu=.2,tide=True):
 aa=a;nn=n2 if tide else 0.;bb=np.array([(aa+4*nn)/2,(aa-2*nn)/2,(aa-4*nn)/2]);tt=np.array([3*nn,0.,-nn])
 S=np.trace(E);norm=np.sum(E*E)
 # All energies normalized by M R²/3, the receiver affine modal mass.
 qg=-1.5*aa*(np.trace(E@E)/5-2*norm/15+S*S/5)
 qk=0.
 for i in range(3):
  j,k=[x for x in range(3) if x!=i]
  Q=-.5*(np.sum(E[:,j]**2)+np.sum(E[:,k]**2))+9/8*(E[j,j]**2+E[k,k]**2)+7/4*E[j,j]*E[k,k]+(E[j,k]+E[k,j])**2/8
  qk+=bb[i]*Q
 qt=-.5*np.sum(tt[:,None]*E*E)
 sym=(E+E.T)/2;tr=np.trace(sym);sq=np.trace(sym@sym)
 lam=E2D*nu/(1-nu*nu);mu=E2D/(2*(1+nu))
 qe=3/(sigma*R*R)*((lam/5+mu/15)*tr*tr+(lam/15+7*mu/15)*sq)
 return float(qg+qk+qt+qe)
def hessian(E2D,tide):
 basis=np.eye(9).reshape(9,3,3);H=np.zeros((9,9))
 for i in range(9):
  H[i,i]=2*energy2(basis[i],E2D,tide=tide)
  for k in range(i):H[i,k]=H[k,i]=energy2(basis[i]+basis[k],E2D,tide=tide)-.5*(H[i,i]+H[k,k])
 return H
def case(E2D,tide=True,attitude_days=None):
 H=hessian(E2D,tide);C=np.zeros((9,9));nn=n if tide else 0.
 for col in range(3):C[col,3+col]=-2*nn;C[3+col,col]=2*nn
 Pa=np.zeros((9,9))
 for ii in range(9):
  unit=np.eye(9)[ii].reshape(3,3);Pa[:,ii]=((unit-unit.T)/2).ravel()
 if attitude_days:
  tau=attitude_days*86400;H=H+Pa/tau**2;C=C+1.4*Pa/tau
 # Scale time bysqrt(a) so numerical zero modes remain well conditioned.
 Z=np.zeros((9,9));I=np.eye(9);eig=np.linalg.eigvals(np.block([[Z,I],[-H/a,-C/math.sqrt(a)]]))*math.sqrt(a)
 maxreal=max(eig.real);positive=maxreal>1e-7*math.sqrt(a)
 return {'tides_and_Coriolis':tide,'extensional_modulus_N_m':E2D,'active_attitude_response_days':attitude_days,'H_eigenvalues_in_units_a':np.linalg.eigvalsh(H/a).tolist(),'state_eigenvalues_per_s':[[float(e.real),float(e.imag)] for e in eig],'positive_growth_detected':bool(positive),'largest_growth_per_s':float(maxreal),'fastest_growth_efold_days':float(1/maxreal/86400) if positive else None,'equivalent_carbon_sheet_kg_m2_at100GPa2000density':E2D/1e11*2000,'scope':'Homogeneous affine shape/orientation modes, circular-orbit linear tide, equal opposite streams and adiabatically conserved circulation action. No local guide/mesh modes or control damping.'}
rows=[case(0.,False)]+[case(E) for E in [0.,1e7,3e7,1e8,3e8]]+[case(1e8,True,days) for days in [30.,100.,300.]]
# Exact no-tide benchmarks: diagonal shear7a/40, symmetric off-diagonal−19a/20,
# uniform breathinga, rotationszero. Check the quadratic energy directly.
E=np.diag([1.,-1.,0.]);diagq=energy2(E,0,tide=False)
E=np.zeros((3,3));E[0,1]=E[1,0]=1.;offq=energy2(E,0,tide=False)
E=np.eye(3);common=energy2(E,0,tide=False)
checks={'diagonal_shape_energy_over_a':diagq/a,'offdiagonal_shear_energy_over_a':offq/a,'common_radial_energy_over_a':common/a}
assert abs(diagq/a-7/40)<1e-12
assert abs(offq/a+19/20)<1e-12
assert abs(common/a-1.5)<1e-12
# Verify quadratic energy against a nonlinear homoeoid integral and ellipse
# perimeter integration. This does not share the polynomial implementation.
zz,ww=np.polynomial.legendre.leggauss(96);zz=(zz+1)/2;ww=ww/2
th=np.linspace(0.,2*math.pi,1024,endpoint=False)
def full_energy(F):
 C=F@F.T;integ=2*sum(w/math.sqrt(np.linalg.det(np.eye(3)+z*z*(C-np.eye(3)))) for z,w in zip(zz,ww))
 grav=-1.5*a*integ;kin=0.
 for i in range(3):
  jk=[x for x in range(3) if x!=i];tt=np.zeros((3,len(th)));tt[jk[0]]=np.cos(th);tt[jk[1]]=np.sin(th)
  perimeter_ratio=np.mean(np.linalg.norm(F@tt,axis=0));kin+=b[i]/perimeter_ratio**2
 tide=-.5*np.sum(t[:,None]*F*F)
 return grav+kin+tide
rng=np.random.default_rng(123);fd=[];h=2e-4;f0=full_energy(np.eye(3))
for i in range(12):
 d=rng.normal(size=(3,3));d/=np.linalg.norm(d)
 exact=(full_energy(np.eye(3)+h*d)+full_energy(np.eye(3)-h*d)-2*f0)/(2*h*h)
 poly=energy2(d,0,tide=True)
 fd.append(abs(exact-poly)/a)
checks['max_nonlinear_quadrature_second_variation_absolute_error_in_a']=max(fd)
assert max(fd)<1e-6
out={'scope':'Bounded global affine dynamics. A three-family stream system needs shear coupling even without tides; membrane rows test a physical coupling scale. This is not a complete network proof.','input_sha256':hashlib.sha256(inp.read_bytes()).hexdigest(),'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'a_s2':a,'solar_n2_s2':n2,'cases':rows,'checks':checks}
Path(__file__).with_name('support-modes.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
print(json.dumps({'checks':checks,'rows':[{k:x[k] for k in ['tides_and_Coriolis','extensional_modulus_N_m','positive_growth_detected','fastest_growth_efold_days','active_attitude_response_days','H_eigenvalues_in_units_a']} for x in rows]},indent=2))
