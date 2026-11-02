"""Long-wavelength axial bias response of a thermally equilibrated ion column.

Linearized Poisson-Boltzmann response at fixed global chemical potentials and
rotation, with prescribed wall voltage. This tests whether a wall gradient
actually penetrates the reacting core; it is not an electron clearing model.
"""
from pathlib import Path
from math import sqrt,pi
import argparse,json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
from electron_wall_orbits import Cylinder,ME,E_CHARGE
ROOT=Path(__file__).parent

def integrate_response(x,chi,xw):
 susceptibility=PchipInterpolator(x,chi)
 x0=x[0];k0=float(susceptibility(x0))
 def rhs(x,y):return [y[1],susceptibility(x)*y[0]-y[1]/x]
 return solve_ivp(rhs,[x0,xw],[1+k0*x0*x0/4,k0*x0/2],rtol=1e-10,atol=1e-12,dense_output=True,max_step=.02)

def response(c):
 eq=c.eq;z=c.src['ZA'][0];xw=c.src['xwall'];x0=eq['x'][0]
 sol=integrate_response(eq['x'],eq['h']+z*z*eq['c'],xw)
 gain=float(sol.y[0,-1]);x=np.linspace(x0,xw,8001);y=sol.sol(x)[0]/gain
 h=np.interp(x,eq['x'],eq['h']);cat=np.interp(x,eq['x'],eq['c']);w=h*cat*x
 average=float(np.trapezoid(w*y,x)/np.trapezoid(w,x))
 return {'axis_wall_voltage_response':1/gain,'reaction_weighted_voltage_response':average,
  'field_transfer':lambda r:float(sol.sol(r/c.scale)[0]/gain),'success':sol.success}

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=ROOT);ap.add_argument('--outdir',type=Path,default=ROOT);args=ap.parse_args()
 plant=next(p for p in json.loads((args.root/'network-plant.json').read_text())['cases'] if p['assumptions']['B']==140 and p['assumptions']['heat']==.7 and p['assumptions']['end_in']==.999)
 out=[];L=plant['active_length_m']
 for key,src in plant['selected_source_profiles'].items():
  c=Cylinder(src);r=response(c);field=r.pop('field_transfer')
  rates=np.interp(src['xwall'],c.eq['x'],c.eq['ihc']);cdf=c.eq['ihc']/rates
  r['processor']=key;r['at_reaction_quantiles']=[]
  for f in [.1,.5,.9,.99]:
   rr=float(np.interp(f,cdf,c.eq['x'])*c.scale);fraction=field(rr)
   # Maximum nonrelativistic collisionless crossing time for any initial axial
   # velocity in a finite interval with constant nonzero acceleration.
   # A trajectory can turn just inside the upstream boundary, giving 2*sqrt(2L/a).
   E=1e3/L*fraction;bound=2*sqrt(2*L*ME/(E_CHARGE*E))
   r['at_reaction_quantiles'].append({'reaction_fraction_inside':f,'radius_m':rr,'wall_voltage_response':fraction,
     'collisionless_time_bound_s_for_1kV_wall_drop':bound})
  out.append(r);print(key,r['axis_wall_voltage_response'],r['reaction_weighted_voltage_response'],flush=True)
 (args.outdir/'axial-bias-screen.json').write_text(json.dumps({'scope':__doc__,'active_length_m':L,'wall_voltage_drop_V':1000.,'processors':out,
 'limits':['Fixed chemical potentials and common rotation; long axial wavelength compared with the channel radius.',
 'The collisionless time bound omits Coulomb scattering, pitch-angle diffusion, radiation, traps and collective fields.',
 'A wall-voltage change must be solved self-consistently before its full axial field is assigned to electrons.']},indent=2,allow_nan=False)+'\n')
