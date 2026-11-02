"""Check linear wall-bias screening against analytic and nonlinear equilibria."""
from pathlib import Path
from math import sqrt
import argparse,json
import numpy as np
from scipy.special import i0
from scipy.interpolate import CubicHermiteSpline
from catalyst_penning import equilibrium
from electron_wall_orbits import Cylinder,ME,E_CHARGE
from axial_bias_screen import response,integrate_response
ROOT=Path(__file__).parent


def main():
 ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=ROOT);ap.add_argument('--outdir',type=Path,default=ROOT);args=ap.parse_args()
 plant=next(p for p in json.loads((args.root/'network-plant.json').read_text())['cases'] if p['assumptions']['B']==140 and p['assumptions']['heat']==.7 and p['assumptions']['end_in']==.999)
 checks=[];derivatives=[]
 def check(name,passed,metric=None):checks.append({'name':name,'passed':bool(passed),'metric':metric})
 x=np.linspace(1e-6,22,5001)
 for chi in [0.,.001,.1,1.5]:
  sol=integrate_response(x,np.full_like(x,chi),22)
  actual=sol.sol(x)[0]/sol.y[0,-1];exact=i0(sqrt(chi)*x)/i0(sqrt(chi)*22)
  error=float(np.max(np.abs(actual/exact-1)))
  check(f'Constant-density modified-Bessel solution chi={chi:g}',sol.success and error<1e-8,error)
 for key,src in plant['selected_source_profiles'].items():
  c=Cylinder(src);ref=response(c);R=src['xwall'];T=src['T_MeV'];z=src['ZA'][0]
  points=np.array([c.eq['x'][0],R*.2,R*.4,R*.7,R]);expect=np.array([ref['field_transfer'](xx*c.scale) for xx in points])
  variants=[]
  for delta in [1e-7,5e-8]:
   potentials=[]
   for sign in [-1,1]:
    h=src['h0']*np.exp(-sign*delta/T);cc=src['c0']*np.exp(-z*sign*delta/T)
    eq=equilibrium(*src['ZA'],h,cc,src['xi'],nout=24001,rtol=2e-12)
    f=CubicHermiteSpline(eq['x'],eq['psi'],eq['dpsi'])
    potentials.append(sign*delta+T*f(points))
   diff=potentials[1]-potentials[0];ratio=diff/diff[-1]
   error=float(np.max(np.abs(ratio/expect-1)))
   variants.append({'axis_perturbation_MV':delta,'maximum_relative_error':error})
   check(f'Finite nonlinear chemical-potential perturbation {key}, delta={delta:g} MV',error<4e-4,error)
  derivatives.append({'processor':key,'linear_axis_response':ref['axis_wall_voltage_response'],'comparisons':variants})
 # Check the purely kinematic constant-acceleration bound without plasma
 # assumptions. Collisions can change residence and are explicitly excluded.
 rng=np.random.default_rng(2891);L=plant['active_length_m'];a=E_CHARGE/ME*.01
 tlimit=2*sqrt(2*L/a);largest=0.
 for _ in range(10000):
  z=rng.uniform(0,L);v=rng.normal()*sqrt(2*a*L)*3
  tplus=(-v+sqrt(v*v+2*a*(L-z)))/a
  t=tplus
  if v<0 and v*v>=2*a*z:
   tminus=(-v-sqrt(v*v-2*a*z))/a
   if tminus>0:t=min(t,tminus)
  largest=max(largest,t/tlimit)
 check('Constant-acceleration absorbing-interval time bound',largest<=1,largest)
 out={'scope':__doc__,'checks':checks,'passed':sum(c['passed'] for c in checks),'total':len(checks),'finite_equilibrium_comparisons':derivatives}
 (args.outdir/'axial-bias-validation.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
 print(out['passed'],'/',out['total'],'checks pass')
 for c in checks:
  if not c['passed']:print('FAILED',c)
 if out['passed']!=out['total']:raise SystemExit(1)


if __name__=='__main__':main()
