"""Independent Fourier, Lorentz, Maxwellian and spatial-scale checks."""
from pathlib import Path
from math import pi,sqrt,exp
import argparse,json,time
import numpy as np
from scipy.integrate import quad,solve_ivp
from scipy.special import k0,k1,erf
from cold_magnetized_response import pulse_fractions,speed_density,conditional_mu2,rate
from electron_wall_orbits import Cylinder,ME,MP,E_CHARGE as E,EPS0
from lepton_collision_screen import local_state

ROOT=Path(__file__).parent


def main():
 ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=ROOT);ap.add_argument('--outdir',type=Path,default=ROOT);args=ap.parse_args()
 start=time.monotonic();checks=[]
 def check(name,ok,metric=None):checks.append({'name':name,'passed':bool(ok),'metric':metric})
 for a in [.001,.03,.1,1.,3.,10.]:
  even=2*quad(lambda t:(1+t*t)**-1.5,0,np.inf,weight='cos',wvar=a,epsabs=1e-11,limlst=300)[0]
  odd=2*quad(lambda t:t*(1+t*t)**-1.5,0,np.inf,weight='sin',wvar=a,epsabs=1e-11,limlst=300)[0]
  err=max(abs(even-2*a*k1(a)),abs(odd-2*a*k0(a)))
  check(f'Coulomb force Fourier transforms a={a}',err<2e-9,err)
 # Independent direct Lorentz ODE, with a finite force-pulse window.
 for a in [.03,1.,10.]:
  theta=.73;phi=.91
  u=np.array([np.sin(theta),0.,np.cos(theta)])
  b=np.array([np.cos(theta)*np.cos(phi),np.sin(phi),-np.sin(theta)*np.cos(phi)])
  T=80.;epsilon=1e-4;omega=-a
  def rhs(t,y):
   force=epsilon*(b+u*t)/(1+t*t)**1.5
   acceleration=force+np.array([omega*y[1],-omega*y[0],0.])
   return np.r_[acceleration,np.dot(force,y[:3]),force]
  sol=solve_ivp(rhs,[-T,T],np.zeros(7),rtol=2e-10,atol=epsilon*1e-12,method='DOP853',max_step=min(.2,.3/a))
  v=sol.y[:3,-1];work=sol.y[3,-1];impulse=sol.y[4:,-1]
  even=2*quad(lambda t:(1+t*t)**-1.5,0,T,weight='cos',wvar=a,epsabs=1e-12)[0]
  odd=2*quad(lambda t:t*(1+t*t)**-1.5,0,T,weight='sin',wvar=a,epsabs=1e-12)[0]
  amplitude=epsilon*(even*(b[0]+1j*b[1])-1j*odd*(u[0]+1j*u[1]))
  expected_perp=abs(amplitude)**2
  expected_parallel=(2*epsilon*b[2]*T/sqrt(1+T*T))**2
  check(f'Direct Lorentz transverse response a={a}',abs(np.sum(v[:2]**2)-expected_perp)<1e-7*epsilon**2,
      float(abs(np.sum(v[:2]**2)-expected_perp)/epsilon**2))
  check(f'Direct Lorentz parallel response a={a}',abs(v[2]**2-expected_parallel)<1e-8*epsilon**2)
  check(f'Direct Lorentz work/energy closure a={a}',abs(work-.5*np.dot(v,v))<1e-7*epsilon**2)
  check(f'Parallel impulse is magnetic-field independent a={a}',abs(v[2]-impulse[2])<1e-12)
 # Impact-azimuth average of the complex Fourier amplitude.
 phis=(np.arange(64)+.5)*2*pi/64
 for mu in [0.,.3,1.]:
  st=sqrt(1-mu*mu);u=np.array([st,0.,mu])
  b=np.column_stack([mu*np.cos(phis),np.sin(phis),-st*np.cos(phis)])
  for a in [.01,.3,1.,4.]:
   amp=a*k1(a)*(b[:,0]+1j*b[:,1])+1j*a*k0(a)*(u[0]+1j*u[1])
   parallel,perp=pulse_fractions(np.array([a]),mu*mu)
   error=max(abs(np.mean(abs(amp)**2)-perp[0]),abs(np.mean(b[:,2]**2)-parallel[0]))
   check(f'Impact-azimuth average mu={mu:g}, a={a:g}',error<1e-13,float(error))
 for U in [0.,.01,1.,3.,10.]:
  vth=1.
  norm=quad(lambda u:speed_density(u,vth,U),0,U+10)[0]
  check(f'Drifting Maxwell speed normalization U/vth={U:g}',abs(norm-1)<2e-11,norm)
  inverse=quad(lambda u:speed_density(u,vth,U)/u,0,U+10)[0]
  expected=2/sqrt(pi) if U==0 else erf(U)/U
  check(f'Drifting Maxwell inverse-speed moment U/vth={U:g}',abs(inverse/expected-1)<2e-11,float(inverse/expected-1))
  # Angular integral around the drift axis gives one of the two transverse
  # moments; the magnetic-field direction lies in this transverse plane.
  speed=1.7;a=2*speed*U
  den=quad(lambda x:exp(a*(x-1)),-1,1)[0]
  value=quad(lambda x:.5*(1-x*x)*exp(a*(x-1)),-1,1)[0]/den
  theory=float(conditional_mu2(np.array([speed]),1.,U)[0])
  check(f'Conditional pitch moment U/vth={U:g}',abs(value-theory)<1e-11,abs(value-theory))
 # Zero-B, fixed-log limit must equal the independent Landau cold-particle tensor.
 for mass,z in [(MP,1),(49*MP,22)]:
  for U in [0.,1e6,6e6,2e7]:
   row=rate(mass,1e19,z,1.6e6,U,0.,5e-4,fixed_log=15.)
   err=abs(row['total_heating_eV_s']/row['NRL_fixed_log_15_cold_heating_eV_s']-1)
   check(f'Zero-field heating matches Landau Z={z}, U={U:g}',err<1e-10,err)
 # Scale gradients independently from nearby reconstructed equilibria.
 plant=next(p for p in json.loads((args.root/'network-plant.json').read_text())['cases'] if p['assumptions']['B']==140 and p['assumptions']['heat']==.7 and p['assumptions']['end_in']==.999)
 c=Cylinder(plant['selected_source_profiles']['22049'])
 scales=json.loads((args.outdir/'magnetized-lepton-scales.json').read_text())
 for row in scales['cases']:
  r=row['radius_m'];h=2e-6*c.scale
  for species,gradient,key in [('H','proton_density_gradient_length_m','proton_density_m-3'),('Ti49','catalyst_density_gradient_length_m','catalyst_density_m-3')]:
   # Evaluate analytic density using the independent potential spline to avoid
   # a grid-interpolation derivative in the comparison.
   def density(rr):
    x=rr/c.scale
    return np.exp(-float(c.spline(x))-(c.eq['ah']*x*x)) if species=='H' else np.exp(-c.eq['Z']*float(c.spline(x))-c.eq['ac']*x*x)
   derivative=(np.log(density(r+h))-np.log(density(r-h)))/(2*h)
   err=abs(1/abs(derivative)/row[gradient]-1)
   check(f'Independent density gradient {species}, fraction={row["reaction_fraction_inside"]}',err<2e-7,err)
 # Quadrature refinement is a numerical comparison, not physical screening validation.
 mid=next(r for r in scales['cases'] if r['reaction_fraction_inside']==.5);local=mid['local_state']
 for s in mid['species']:
  params=(s['mass_kg'],s['density_m-3'],s['charge'],local['ion_temperature_eV'],local['relative_azimuthal_m_s'],137.14908499885485,local['ion_static_screening_m'])
  a=rate(*params,nv=256,nb=128);b=rate(*params,nv=512,nb=256)
  for metric in ['parallel_heating_eV_s','perpendicular_heating_eV_s']:
   err=abs(a[metric]/b[metric]-1)
   check(f'Cold rate quadrature refinement {s["name"]} {metric}',err<2e-3,err)
 out={'scope':__doc__,'passed':sum(c['passed'] for c in checks),'total':len(checks),'checks':checks,'seconds':time.monotonic()-start}
 (args.outdir/'cold-magnetized-validation.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
 print(out['passed'],'/',out['total'],'checks pass',out['seconds'],'s')
 for c in checks:
  if not c['passed']:print('FAILED',c)
 if out['passed']!=out['total']:raise SystemExit(1)


if __name__=='__main__':main()
