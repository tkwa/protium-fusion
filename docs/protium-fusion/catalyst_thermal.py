"""Conditional main-cycle thermal screen for arbitrary (p,n)/EC catalysts.

Only the pn -> cold decay -> catalyst and neutron -> D -> closed D-burn
cycle is closed. All destructive channels have a specified inclusive net
restoration cost, and are diagnostics rather than exclusive reaction rates.
No useful nuclear-energy credit is taken for those side channels.
"""
from pathlib import Path
from functools import lru_cache
from math import pi,sqrt,log
import json,gzip
import numpy as np
from catalyst_catalog import thermal,nucmass,atomicmass,nuclear
from nonneutral_screen import MEV,C,Q,DD_CHARGED_PER_NEW_N,DD_THERMAL_PER_NEW_N
ROOT=Path(__file__).parent
CAT={r['target']:r for r in json.loads((ROOT/'catalyst-catalog.json').read_text())['candidates'] if 'thermal_cases' in r}
CURVES=json.loads(gzip.decompress((ROOT/'catalyst-pn-curves.json.gz').read_bytes()))['evaluations']

@lru_cache(maxsize=512)
def case(target,T,inelastic_energy_fraction=.5,pn_gamma_fraction=0.):
 r=CAT[target];Z,A=r['Z'],r['A'];cs=CURVES[r['source_file']]['channels']
 pnkey=f'{1000*(Z+1)+A}_0';same=f'{1000*Z+A}_0'
 sv=thermal(cs[pnkey],T,Z,A);Ecm=thermal(cs[pnkey],T,Z,A,1)/sv
 kp=Ecm+1.5*T
 mi=nucmass(1,1)+nucmass(Z,A);mo=nucmass(0,1)+nucmass(Z+1,A)
 kcom=1.5*T*mi/mo
 # The gamma fraction brackets the unknown excitation cascade. Small photon
 # recoil is omitted; all cascade energy thermalizes outside the source.
 avail=max(0,kp+r['Q_pn_MeV']-kcom)
 gamma=pn_gamma_fraction*avail;rel=avail-gamma
 kn=nucmass(0,1)/mo*kcom+nucmass(Z+1,A)/mo*rel
 kdaughter=kp+r['Q_pn_MeV']-kn-gamma
 insv=thermal(cs[same],T,Z,A) if same in cs else 0.
 inel=insv/sv
 inel_E=inelastic_energy_fraction*thermal(cs[same],T,Z,A,1)/insv if insv else 0.
 # For EC-only return, taking the full atomic decay Q as escaping neutrino
 # energy underestimates recoverable heat by the atomic cascade energy.
 nu=.706 if target=='13C' else r['decay_Q_atomic_MeV']
 qheat=.5*Q['net_H']-nu
 prep=kp+inel*inel_E
 charged=DD_CHARGED_PER_NEW_N+kdaughter
 heat=DD_THERMAL_PER_NEW_N+Q['np_capture']+kn+gamma+r['decay_Q_atomic_MeV']-nu+inel*inel_E
 sides={}
 for key,tab in cs.items():
  za,iso=map(int,key.split('_'));zz,aa=divmod(za,1000)
  if aa<=4 or key in [pnkey,same]:continue
  val=thermal(tab,T,Z,A)/sv
  if val>1e-10:sides[key]=val
 incln=thermal(cs['1_0'],T,Z,A)/sv if '1_0' in cs else None
 incld=thermal(cs['1002_0'],T,Z,A)/sv if '1002_0' in cs else None
 tau=nuclear(Z+1,A).tau/log(2)
 return {'target':target,'Z':Z,'A':A,'daughter':r['daughter'],'T_MeV':T,
  'pn_sigma_v_m3_s':sv,'mean_reacting_Ecm_MeV':Ecm,'pair_kinetic_removed_MeV':kp,
  'neutron_kinetic_MeV':kn,'daughter_kinetic_MeV':kdaughter,
  'assumed_neutrino_MeV':nu,'nuclear_heat_MeV':qheat,
  'pn_cascade_gamma_MeV':gamma,'pn_gamma_fraction':pn_gamma_fraction,
  'inelastic_events_per_pn':inel,'inelastic_excitation_MeV':inel_E,
  'inelastic_energy_fraction':inelastic_energy_fraction,
  'input_ion_energy_MeV':prep,'recoverable_ion_energy_MeV':charged,'heat_MeV':heat,
  'energy_closure_MeV':charged+heat-prep-qheat,
  'other_heavy_residuals_per_pn':sides,'side_residual_sum_per_pn':sum(sides.values()),
  'inclusive_neutrons_per_pn':incln,'inclusive_D_per_pn':incld,
  'neutral_daughter_mean_life_s':tau,
  'decay_inventory_kg_per_pn_s':tau*atomicmass(Z+1,A)*MEV/C**2,
  'limits':['main pn weak-conversion plus D-burn cycle only; side restoration unspecified',
    'PADF independent residual-production data are not a closed reaction network',
    'full stripping, prompt extraction, neutral cold decay and ion reinjection assumed',
    'nonrelativistic thermal distribution; pn excitation-gamma fraction is a bracket',
    'D-reactor auxiliaries and actual energy-recovery hardware excluded']}

def electricity(d,heater=.95,direct=.95,heat=.4,side_cost=50.):
 g=direct+(1-direct)*heat
 return g*d['recoverable_ion_energy_MeV']+heat*d['heat_MeV']-d['input_ion_energy_MeV']/heater-side_cost*d['side_residual_sum_per_pn']

if __name__=='__main__':
 out=[]
 for name in ['13C','49Ti','51V','53Cr','55Mn','57Fe','65Cu','93Nb']:
  for T in [.4,.6,.8,1.,1.2,1.5,2.,2.5]:
   d=case(name,T);d['net_MeV_per_pn']=electricity(d);out.append(d)
   print(name,T,'sv',f"{d['pn_sigma_v_m3_s']:.3g}",'side',round(d['side_residual_sum_per_pn'],5),'inel',round(d['inelastic_events_per_pn'],4),'net',round(d['net_MeV_per_pn'],3))
 (ROOT/'catalyst-thermal.json').write_text(json.dumps({'scope':'Conditional cycle, no trap or whole plant;95% heater and ion recovery,40% heat,50MeV net cost per other heavy residual. C13 pd counted as generic side cost, unlike dedicated C13 model.','cases':out},indent=2,allow_nan=False)+'\n')
 print('closure',max(abs(d['energy_closure_MeV']) for d in out))
