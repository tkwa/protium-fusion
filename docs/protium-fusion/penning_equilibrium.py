"""Two-ion Poisson--Boltzmann equilibrium screen, not a working Penning reactor.

A long conducting cylinder with uniform imposed B, common temperature and
rigid rotation. Infinite axial length; end losses and driven transport are absent. Ideal global thermal equilibrium
is stable; extraction and boundary asymmetries require separate treatment. Wall flux is a full-Maxwellian-tail
estimate, not a kinetic absorbing-boundary solution. Do not interpret an
axisymmetric equilibrium or a field-only mass allocation as viability.
"""
from pathlib import Path
import json
from math import pi,sqrt
from functools import lru_cache
import numpy as np
from scipy.integrate import solve_ivp
from thermal_neutron_screen import MASS,MEV,case,electricity
from nonneutral_screen import EPS0,MU0,E_CHARGE,C

ROOT=Path(__file__).parent
MP=MASS['p']*MEV/C**2
MC=MASS['c13']*MEV/C**2
MRC=MC/MP
KAPPA=6/MRC
thermal_case=lru_cache(maxsize=64)(case)


def equilibrium(h0,c0,xi=.5,xmax=50,rtol=1e-9,nout=1501):
 """n/nscale at axis h0,c0; omega=xi*omega_c,C; x=r/l."""
 ah=(KAPPA*xi-(KAPPA*xi)**2)/2
 ac=MRC*KAPPA**2*xi*(1-xi)/2
 def densities(x,psi):
  return h0*np.exp(np.clip(-psi-ah*x*x,-745,50)),c0*np.exp(np.clip(-6*psi-ac*x*x,-745,50))
 def rhs(x,y):
  psi,w=y[:2];h,c=densities(x,psi)
  # Integrated number, overlap, reaction radius squared, and total rotation.
  return [w,-w/x-h-6*c,h*x,c*x,h*c*x,h*c*x**3,h*x**3,c*x**3]
 def blowup(x,y):
  h,c=densities(x,y[0]);return 1000-h-6*c
 blowup.terminal=True;blowup.direction=-1
 x0=1e-6;s=h0+6*c0
 y0=[-s*x0*x0/4,-s*x0/2,h0*x0*x0/2,c0*x0*x0/2,h0*c0*x0*x0/2,
     h0*c0*x0**4/4,h0*x0**4/4,c0*x0**4/4]
 sol=solve_ivp(rhs,[x0,xmax],y0,rtol=rtol,atol=rtol*1e-3,dense_output=True,events=blowup,max_step=.3)
 x=np.linspace(x0,sol.t[-1],nout);y=sol.sol(x);h,c=densities(x,y[0])
 return {'x':x,'psi':y[0],'dpsi':y[1],'h':h,'c':c,'ih':y[2],'ic':y[3],'ihc':y[4],
         'ihcx2':y[5],'ihx2':y[6],'icx2':y[7], 'h0':h0,'c0':c0,'xi':xi,'ah':ah,'ac':ac,
         'success':bool(sol.success and not sol.t_events[0].size), 'solver_rtol':rtol}


def at_wall(eq,xwall,T=.6,B=100,heat_eff=.95,direct_eff=.95,thermal_eff=.4,
            side_cost=50,strength=1e6,mass_multiplier=5,packing=.7,
            potential_in_eff=1.,potential_out_eff=1.,radiation_escape_fraction=0.):
 """Evaluate every column per unit length; field volume includes its empty edge."""
 def v(key):return float(np.interp(xwall,eq['x'],eq[key]))
 tj=T*MEV;l=sqrt(tj*MP)/(E_CHARGE*B);ns=EPS0*B*B/MP
 omega=eq['xi']*6*E_CHARGE*B/MC;R=xwall*l
 nc=thermal_case(T);base=electricity(nc,heat_eff,direct_eff,thermal_eff,side_cost)['net_MeV_per_pn']
 Rpn=2*pi*l*l*ns*ns*nc['pn_sigma_v_m3_s']*v('ihc')
 nh=2*pi*l*l*ns*v('ih');nC=2*pi*l*l*ns*v('ic')
 mx2=v('ihcx2')/v('ihc');r2=mx2*l*l
 # Products are extracted at their reaction radius. Common rotation contributes
 # a mass-weighted neutron kinetic share. Injected/extracted charge has the same
 # electrostatic potential at this location in a reversible closed circuit;
 # actual circuits and finite-efficiency potential handling remain unmodeled.
 rotpair=.5*(MP+MC)*omega*omega*r2/MEV
 neutronrot=rotpair*MASS['n']/(MASS['p']+MASS['c13'])
 event_ratio=1+nc['pd_per_pn']+nc['knockout_neutron_per_pn']
 rotprep=event_ratio*rotpair
 rotneutron=neutronrot*(1+nc['knockout_neutron_per_pn'])
 rotcharged=rotprep-rotneutron
 g=direct_eff+(1-direct_eff)*thermal_eff
 rotnet=g*rotcharged+thermal_eff*rotneutron-rotprep/heat_eff
 fluxH=2*pi*R*ns*v('h')*sqrt(tj/(2*pi*MP))
 fluxC=2*pi*R*ns*v('c')*sqrt(tj/(2*pi*MC))
 # Outgoing flux weights the radial kinetic term, giving total random KE 2T.
 wallH=2*T+.5*MP*omega**2*R*R/MEV
 wallC=2*T+.5*MC*omega**2*R*R/MEV
 wallcirculation=(fluxH*wallH+fluxC*wallC)/Rpn
 wallpenalty=wallcirculation*(1/heat_eff-g)
 # Each extracted pn/pd/p,pn event carries total charge +7e. Injection and
 # extraction at the same local potential cancel only for a reversible circuit.
 inside=eq['x']<xwall
 xx=np.r_[eq['x'][inside],xwall]
 hh=np.interp(xx,eq['x'],eq['h']);cc=np.interp(xx,eq['x'],eq['c'])
 pp=np.interp(xx,eq['x'],eq['psi']);dd=np.interp(xx,eq['x'],eq['dpsi'])
 meanpsi=float(np.trapezoid(hh*cc*pp*xx,xx)/v('ihc'))
 potential_work=7*T*(meanpsi-v('psi'))*event_ratio
 potential_recovered=potential_out_eff+(1-potential_out_eff)*thermal_eff
 potentialpenalty=potential_work*(1/potential_in_eff-potential_recovered)
 # Vacuum Larmor estimate. A conducting sub-cutoff tube can strongly alter
 # emission and reabsorption; escape_fraction is a sensitivity, not a spectrum.
 er=-tj/(E_CHARGE*l)*dd;rad=0.
 for dens,z,m in [(hh,1,MP),(cc,6,MC)]:
  force2=(er-omega*xx*l*B)**2+2*B*B*tj/m
  pone=(z*E_CHARGE)**4/(6*pi*EPS0*m*m*C**3)*force2
  rad+=2*pi*l*l*ns*np.trapezoid(dens*pone*xx,xx)/Rpn/MEV
 radpenalty=radiation_escape_fraction*rad*(1/heat_eff-thermal_eff)
 net=base+rotnet-wallpenalty-potentialpenalty-radpenalty
 mf=mass_multiplier*(B*B/(2*MU0))*pi*R*R/strength/packing
 # E=-T/e/l*dpsi, potential measured between axis and wall.
 ewall=-tj/(E_CHARGE*l)*v('dpsi');dV=-T*1e6*v('psi')
 return {'T_MeV':T,'B_T':B,'xi':eq['xi'],'xwall':xwall,'h0':eq['h0'],'c0':eq['c0'],
  'radius_m':R,'reaction_rate_per_m_s':Rpn,'proton_inventory_per_m':nh,'carbon_inventory_per_m':nC,
  'mean_reaction_x2':mx2,'mean_pair_rotation_MeV':rotpair,
  'rotation_net_MeV_per_pn':rotnet,'base_net_MeV_per_pn':base,
  'wall_ion_energy_circulated_MeV_per_pn':wallcirculation,'wall_net_penalty_MeV_per_pn':wallpenalty,
  'net_MeV_per_pn':net,'field_mass_per_m_kg':mf,
  'electrostatic_work_circulated_MeV_per_pn':potential_work,
  'electrostatic_net_penalty_MeV_per_pn':potentialpenalty,
  'potential_in_efficiency':potential_in_eff,'potential_out_efficiency':potential_out_eff,
  'vacuum_ion_Larmor_MeV_per_pn':float(rad),'assumed_radiation_escape_fraction':radiation_escape_fraction,
  'radiation_net_penalty_MeV_per_pn':float(radpenalty),'field_mass_only_net_W_kg':Rpn*net*MEV/mf,
  'proton_wall_lifetime_s':nh/fluxH if fluxH else 1e300,
  'carbon_wall_lifetime_s':nC/fluxC if fluxC else 1e300,
  'wall_E_V_m':ewall,'axis_to_wall_potential_V':dV,
  'proton_tail_ratio':v('h')/eq['h0'],'carbon_tail_ratio':v('c')/eq['c0'],
  'uniform_cold_optimum_overlap_ratio':(2*v('ihc')/xwall**2)/((3/MRC/2)*(3/MRC/12)),
  'rotation_speed_wall_over_c':omega*R/C}


def main():
 eq=equilibrium(3/MRC/2,3/MRC/12,.5)
 print('canonical success',eq['success'],'h0,c0',eq['h0'],eq['c0'])
 rows=[at_wall(eq,x) for x in [10,15,20,25,30,35,40,50]] if eq['success'] else []
 out={'scope':'Two-species cylindrical thermal equilibrium and full-tail wall estimate; field mass only. Ideal closed thermal equilibrium is stable; driven extraction, boundary asymmetries and end plugs are unmodeled.',
      'baseline':rows,'profile':{k:eq[k].tolist() for k in ['x','psi','h','c']}}
 (ROOT/'penning-equilibrium.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
 for row in rows:
  print('x',row['xwall'],'W/kg',round(row['field_mass_only_net_W_kg'],8),'net',round(row['net_MeV_per_pn'],4),
        'rot',round(row['mean_pair_rotation_MeV'],3),'wall',round(row['wall_net_penalty_MeV_per_pn'],3),
        'overlap',round(row['uniform_cold_optimum_overlap_ratio'],4),'potentialMV',round(row['axis_to_wall_potential_V']/1e6,2))

if __name__=='__main__':main()
