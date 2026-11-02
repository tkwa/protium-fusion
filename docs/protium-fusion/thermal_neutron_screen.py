"""Thermal nonneutral C13 neutron source linked to a closed D-fusion cycle.

Common Maxwellian ion temperature; all pn/pd products are promptly extracted.
Coulomb thermalization is not charged as beam stopping. Reaction-selected
kinetic energies, external neutron capture and minimum inelastic radiation
are counted. This is an optimistic process screen, not a confinement model.
"""
from pathlib import Path
from math import pi,sqrt
import json
import numpy as np
from numpy.polynomial.legendre import leggauss
from nonneutral_screen import (DATA,Q,mass,MEV,MP,MC,MU,EPS0,MU0,
 DD_CHARGED_PER_NEW_N,DD_THERMAL_PER_NEW_N)

ROOT=Path(__file__).parent
ME_MEV=.51099895
MASS={'p':mass('p')-ME_MEV,'c13':mass('c13')-6*ME_MEV,
      'n':mass('n'),'n13':mass('n13')-7*ME_MEV}
MU_PC=(MASS['p']*MASS['c13']/(MASS['p']+MASS['c13']))*MEV/299792458.**2
ALPHA=MASS['c13']/(MASS['p']+MASS['c13'])
NODES,WEIGHTS=leggauss(24)


def thermal_integrals(T,za,extra_power=0):
 """Integrate sigma(E_lab) E_cm^(1+extra_power) exp(-E_cm/T)."""
 d=DATA['channels'][f'{za}_0']
 e=ALPHA*np.array(d['energy_MeV']);s=np.array(d['sigma_barn'])
 lo=e[:-1,None];hi=e[1:,None];half=(hi-lo)/2
 x=(hi+lo)/2+half*NODES;sig=(s[:-1,None]+s[1:,None])/2+(s[1:,None]-s[:-1,None])/2*NODES
 val=np.sum(half*WEIGHTS*sig*x**(1+extra_power)*np.exp(-x/T))
 return float(val)


def moments(T,za):
 i0=thermal_integrals(T,za)
 return {'sigma_v_m3_s':sqrt(8*MEV/(pi*MU_PC))*T**-1.5*i0*1e-28,
         'mean_reacting_Ecm_MeV':thermal_integrals(T,za,1)/i0 if i0 else 0.}


def neutron_kinetic(T,Ecm,Qrx):
 # The Maxwellian COM velocity is independent of relative velocity. Its
 # cross term with outgoing relative velocity averages to zero even for
 # anisotropic differential scattering. Tiny mass-defect recoil is retained.
 mi=MASS['p']+MASS['c13'];mo=MASS['n']+MASS['n13']
 kcm=1.5*T*mi/mo;erel=Ecm+Qrx+1.5*T-kcm
 return MASS['n']/mo*kcm+MASS['n13']/mo*erel


def case(T,excitation=3.5,knockout_partition=.5):
 m={za:moments(T,za) for za in [1,1002,5010,6012,6013,7013,7014]}
 pn=m[7013]['sigma_v_m3_s'];ratio={za:x['sigma_v_m3_s']/pn for za,x in m.items()}
 kpn=m[7013]['mean_reacting_Ecm_MeV']+1.5*T
 kpd=m[1002]['mean_reacting_Ecm_MeV']+1.5*T
 kn=neutron_kinetic(T,m[7013]['mean_reacting_Ecm_MeV'],Q['pn'])
 kN=kpn+Q['pn']-kn
 y=ratio[1002];inel=ratio[6013]*excitation
 # Extra inclusive neutrons are provisionally assigned to p,pn below other
 # major thresholds. Both charged products are extracted; C12 captures that
 # neutron and both stable ions are returned. This adds no new D or net Q.
 knock=max(0,ratio[1]-1)
 if knock>1e-10:
  e_knock=(m[1]['sigma_v_m3_s']*m[1]['mean_reacting_Ecm_MeV']-
          pn*m[7013]['mean_reacting_Ecm_MeV'])/(pn*knock)
  kk=e_knock+1.5*T
  # A bracketed share of the available kinetic energy, not measured spectra.
  knout=knockout_partition*max(0,kk-Q['c12_capture'])
 else: e_knock=kk=knout=0.
 prep=kpn+y*kpd+inel+knock*kk
 charged=(DD_CHARGED_PER_NEW_N+kN+y*(kpd+Q['pd'])+
          knock*(kk-Q['c12_capture']-knout))
 heat=(DD_THERMAL_PER_NEW_N+kn+Q['n13_decay_atomic']-Q['n13_neutrino']+
       (1-y)*Q['np_capture']+y*Q['c12_capture']+inel+
       knock*(Q['c12_capture']+knout))
 return {'ion_temperature_MeV':T,'pn_sigma_v_m3_s':pn,
  'mean_reacting_Ecm_MeV':m[7013]['mean_reacting_Ecm_MeV'],
  'mean_pair_kinetic_removed_MeV_per_pn':kpn,
  'pn_neutron_kinetic_MeV':kn,'N13_extracted_kinetic_MeV':kN,
  'pd_per_pn':y,'raw_C12_residual_per_pn':ratio[6012],
  'B10_per_pn':ratio[5010],'N14_per_pn':ratio[7014],
  'inelastic_per_pn':ratio[6013],'assumed_excitation_MeV':excitation,
  'knockout_neutron_per_pn':knock,'assumed_knockout_neutron_energy_share':knockout_partition,
  'thermal_preparation_MeV_per_pn':prep,
  'recoverable_ion_MeV_per_pn':charged,'heat_MeV_per_pn':heat,
  'neutron_capture_substitution_valid':y<=1,
  'energy_closure_MeV':charged+heat-prep-Q['heat_per_new_neutron'],
  'limits':['prompt product extraction and cold catalyst return',
   'no additional Coulomb loss for a single common thermal ion population',
   'specified ion-heating and recovery efficiencies; no actual converter',
   'rotation, electrostatic-potential circulation, confinement and electron loads absent',
   'deuterium reactor confinement, startup and auxiliaries absent',
   'independent PADF production channels reconciled only in a reduced network',
   'all neutrons needed for capture are retained; no parasitic captures',
   'B10/N14 restoration remains a net-cost parameter',
   'higher nuclear channels and neutron spectra require a consistent evaluation',
   'nonrelativistic Maxwellian kinematics; relativistic corrections unmodeled' ]}


def electricity(d,eta_heater=.8,eta_direct=.9,eta_thermal=.4,side_cost=50.,other=0.):
 g=eta_direct+(1-eta_direct)*eta_thermal
 w=(g*d['recoverable_ion_MeV_per_pn']+eta_thermal*d['heat_MeV_per_pn']-
    d['thermal_preparation_MeV_per_pn']/eta_heater-
    side_cost*(d['B10_per_pn']+d['N14_per_pn'])-other)
 return {'ion_heating_efficiency':eta_heater,'ion_direct_recovery_fraction':eta_direct,
  'heat_conversion_efficiency':eta_thermal,'net_side_restoration_cost_MeV':side_cost,
  'other_net_load_MeV_per_pn':other,'net_MeV_per_pn':w}


def field_only(d,w,B=100,fill=.2,strength=1e6,mult=5):
 S=6*EPS0*B*B/(2*MC);h=S/2;c=S/12
 r=h*c*d['pn_sigma_v_m3_s'];u=B*B/(2*MU0)
 return {'B_T':B,'plasma_fill_fraction':fill,'pn_rate_per_plasma_m3_s':r,
         'optimistic_net_W_kg_field_mass_only':fill*r*w*MEV/(mult*u/strength)}


def main():
 rows=[]
 for T in [.1,.2,.3,.4,.5,.6,.8,1.,1.2,1.5,2.,3.]:
  d=case(T)
  d['energy_cases']=[electricity(d,h,dc,.4,s) for h,dc,s in [(.5,.9,0),(.8,.9,0),(.8,.9,50),(.95,.95,50),(.8,.9,230)]]
  w=d['energy_cases'][2]['net_MeV_per_pn']
  d['field_mass_only']=field_only(d,w)
  rows.append(d)
 out={'scope':'Common-temperature nonneutral source, prompt product extraction, separate closed D reactor; no complete plant.',
      'source_evaluation':DATA['evaluation'],'source_sha256':DATA['source_sha256'],
      'cases':rows}
 (ROOT/'thermal-neutron-screen.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
 for d in rows:
  print('T',d['ion_temperature_MeV'],'sv',f"{d['pn_sigma_v_m3_s']:.3e}",'pd',round(d['pd_per_pn'],3),'B10',round(d['B10_per_pn'],3),'prep',round(d['thermal_preparation_MeV_per_pn'],2),'net',[round(v['net_MeV_per_pn'],2) for v in d['energy_cases']],'fieldW/kg',f"{d['field_mass_only']['optimistic_net_W_kg_field_mass_only']:.3g}")
 print('closure',max(abs(d['energy_closure_MeV']) for d in rows))

if __name__=='__main__':main()
