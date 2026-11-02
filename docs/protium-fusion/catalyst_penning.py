"""Generalized two-species nonneutral equilibrium; conditional field-mass screen.

Same long-cylinder Poisson--Boltzmann model as penning_equilibrium.py,
with arbitrary fully stripped catalyst charge and mass. Closed equilibria
are not driven reactor solutions. All end, collection, coil and shielding
engineering remains separate.
"""
from pathlib import Path
from math import pi,sqrt
from functools import lru_cache
import json,argparse
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar
from catalyst_thermal import case,electricity,CAT,nucmass,MEV,C
from nonneutral_screen import EPS0,MU0,E_CHARGE
ROOT=Path(__file__).parent
MP=nucmass(1,1)*MEV/C**2


def equilibrium(Z,A,h0,c0,xi=.5,xmax=60,rtol=1e-8,nout=2001):
 M=nucmass(Z,A)/nucmass(1,1);kap=Z/M
 ah=(kap*xi-(kap*xi)**2)/2;ac=M*kap**2*xi*(1-xi)/2
 def densities(x,psi):
  return h0*np.exp(np.clip(-psi-ah*x*x,-745,50)),c0*np.exp(np.clip(-Z*psi-ac*x*x,-745,50))
 def rhs(x,y):
  psi,w=y[:2];h,c=densities(x,psi)
  return [w,-w/x-h-Z*c,h*x,c*x,h*c*x,h*c*x**3,h*x**3,c*x**3]
 def blowup(x,y):
  h,c=densities(x,y[0]);return 1000-h-Z*c
 blowup.terminal=True;blowup.direction=-1
 x0=1e-6;s=h0+Z*c0
 y0=[-s*x0*x0/4,-s*x0/2,h0*x0*x0/2,c0*x0*x0/2,h0*c0*x0*x0/2,h0*c0*x0**4/4,h0*x0**4/4,c0*x0**4/4]
 sol=solve_ivp(rhs,[x0,xmax],y0,rtol=rtol,atol=rtol*1e-3,dense_output=True,events=blowup,max_step=.3)
 x=np.linspace(x0,sol.t[-1],nout);y=sol.sol(x);h,c=densities(x,y[0])
 return {'Z':Z,'A':A,'M':M,'kappa':kap,'x':x,'psi':y[0],'dpsi':y[1],'h':h,'c':c,
  'ih':y[2],'ic':y[3],'ihc':y[4],'ihcx2':y[5],'ihx2':y[6],'icx2':y[7],
  'h0':h0,'c0':c0,'xi':xi,'ah':ah,'ac':ac,
  'success':bool(sol.success and not sol.t_events[0].size),'solver_rtol':rtol}


def at_wall(eq,xwall,target,T=1.,B=100,heater=.95,direct=.95,heat=.4,side_cost=50.,
 strength=1e6,mass_multiplier=5,packing=.7,potential_in=.999,potential_out=.999,
 radiation_escape=1.,inelastic_fraction=.5,pn_gamma_fraction=0.):
 def v(k):return float(np.interp(xwall,eq['x'],eq[k]))
 def lifetime(n,f):return 1e300 if f<=n/1e300 else n/f
 Z,A=eq['Z'],eq['A'];MC=MP*eq['M'];tj=T*MEV
 l=sqrt(tj*MP)/(E_CHARGE*B);ns=EPS0*B*B/MP;om=eq['xi']*Z*E_CHARGE*B/MC;R=xwall*l
 nc=case(target,T,inelastic_fraction,pn_gamma_fraction)
 base=electricity(nc,heater,direct,heat,side_cost)
 Rpn=2*pi*l*l*ns*ns*nc['pn_sigma_v_m3_s']*v('ihc')
 nh=2*pi*l*l*ns*v('ih');nC=2*pi*l*l*ns*v('ic')
 r2=l*l*v('ihcx2')/v('ihc')
 rotpair=.5*(MP+MC)*om**2*r2/MEV
 rotneutron=rotpair*nucmass(0,1)/(nucmass(1,1)+nucmass(Z,A))
 g=direct+(1-direct)*heat
 rotnet=g*(rotpair-rotneutron)+heat*rotneutron-rotpair/heater
 fh=2*pi*R*ns*v('h')*sqrt(tj/(2*pi*MP));fc=2*pi*R*ns*v('c')*sqrt(tj/(2*pi*MC))
 wallcirc=(fh*(2*T+.5*MP*om**2*R*R/MEV)+fc*(2*T+.5*MC*om**2*R*R/MEV))/Rpn
 wallpenalty=wallcirc*(1/heater-g)
 inside=eq['x']<xwall;xx=np.r_[eq['x'][inside],xwall]
 hh=np.interp(xx,eq['x'],eq['h']);cc=np.interp(xx,eq['x'],eq['c']);pp=np.interp(xx,eq['x'],eq['psi']);dd=np.interp(xx,eq['x'],eq['dpsi'])
 meanpsi=float(np.trapezoid(hh*cc*pp*xx,xx)/v('ihc'))
 potwork=(Z+1)*T*(meanpsi-v('psi'))
 potpenalty=potwork*(1/potential_in-potential_out-(1-potential_out)*heat)
 er=-tj/(E_CHARGE*l)*dd;rad=0.
 for den,z,m in [(hh,1,MP),(cc,Z,MC)]:
  force2=(er-om*xx*l*B)**2+2*B*B*tj/m
  pone=(z*E_CHARGE)**4/(6*pi*EPS0*m*m*C**3)*force2
  rad+=2*pi*l*l*ns*np.trapezoid(den*pone*xx,xx)/Rpn/MEV
 radpenalty=radiation_escape*rad*(1/heater-heat)
 net=base+rotnet-wallpenalty-potpenalty-radpenalty
 mf=mass_multiplier*B*B/(2*MU0)*pi*R*R/strength/packing
 return {'target':target,'T_MeV':T,'B_T':B,'xi':eq['xi'],'h0':eq['h0'],'c0':eq['c0'],
 'xwall':xwall,'radius_m':R,'reaction_rate_per_m_s':Rpn,'proton_inventory_per_m':nh,'catalyst_inventory_per_m':nC,
 'base_net_MeV_per_pn':base,'pair_rotation_MeV':rotpair,'rotation_net_MeV_per_pn':rotnet,
 'wall_ion_energy_circulated_MeV_per_pn':wallcirc,'wall_penalty_MeV_per_pn':wallpenalty,
 'potential_work_MeV_per_pn':potwork,'potential_penalty_MeV_per_pn':potpenalty,
 'vacuum_ion_Larmor_MeV_per_pn':float(rad),'radiation_penalty_MeV_per_pn':float(radpenalty),
 'net_MeV_per_pn':float(net),'positive_net_case':bool(net>0),'field_mass_per_m_kg':mf,'field_mass_only_net_W_kg':float(Rpn*net*MEV/mf),
 'decay_inventory_per_m_kg':Rpn*nc['decay_inventory_kg_per_pn_s'],
 'axis_to_wall_V':-T*1e6*v('psi'),'wall_E_V_m':-tj/(E_CHARGE*l)*v('dpsi'),
 'proton_wall_lifetime_s':lifetime(nh,fh),'catalyst_wall_lifetime_s':lifetime(nC,fc),
 'rotation_speed_wall_over_c':om*R/C,
 'assumptions':{'heater':heater,'direct':direct,'heat':heat,'net_side_restoration_MeV':side_cost,
  'strength_J_kg':strength,'mass_multiplier':mass_multiplier,'packing':packing,
  'potential_in':potential_in,'potential_out':potential_out,'radiation_escape':radiation_escape,
  'inelastic_fraction':inelastic_fraction,'pn_gamma_fraction':pn_gamma_fraction}}


def scan(target,T=1.):
 Z,A=CAT[target]['Z'],CAT[target]['A'];M=nucmass(Z,A)/nucmass(1,1);S=Z/(2*M)
 allr=[];bad=0
 for xi in [.2,.3,.4,.5,.6,.7,.8]:
  for frac in [.15,.3,.5,.7,.9,1.,1.1,1.3]:
   for share in [.3,.5,.7,.85,.95]:
    eq=equilibrium(Z,A,S*frac*share,S*frac*(1-share)/Z,xi)
    if not eq['success']:bad+=1;continue
    fn=lambda x:-at_wall(eq,x,target,T=T)['field_mass_only_net_W_kg']
    opt=minimize_scalar(fn,bounds=(10,60),method='bounded',options={'xatol':.01})
    row=at_wall(eq,opt.x,target,T=T);row['central_charge_fraction_of_cold_Brillouin']=frac;row['central_proton_charge_share']=share
    allr.append(row)
 allr.sort(key=lambda r:r['field_mass_only_net_W_kg'],reverse=True)
 return {'target':target,'scope':'Main cycle plus specified inclusive side cost, ideal closed-cylinder equilibrium; field mass only. Finite-length confinement, extraction, coil, shielding and D reactor not engineered.','attempted':len(allr)+bad,'blowups':bad,'best':allr[:20]}

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--target',default='49Ti');p.add_argument('--T',type=float,default=1.);a=p.parse_args()
 out=scan(a.target,a.T)
 (ROOT/f'catalyst-penning-{a.target}.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
 print(json.dumps({k:v for k,v in out.items() if k!='best'}));print(json.dumps(out['best'][0],indent=2))
