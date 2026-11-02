"""Conditional C13(p,n) neutron-catalysis screen, not a reactor design.

The useful cycle is 4 H -> He4, with two N13 beta decays. Ion scattering,
beam consumption, competing (p,d), minimum inelastic loss, and component
energy recovery are included. Other side-channel restoration is an explicit
unproved penalty. PADF production channels are independently evaluated and
not a consistent exclusive reaction library; the inconsistency is exposed.
"""
from pathlib import Path
import json
from math import pi,sqrt
import numpy as np
import pynucastro as pyn

ROOT=Path(__file__).parent
MEV=1.602176634e-13
E_CHARGE=1.602176634e-19
EPS0=8.8541878128e-12
MU0=1.25663706212e-6
C=299792458.
MP=1.67262192369e-27
ME=9.1093837015e-31
MU=1.66053906660e-27
MC=13*MU
HBARC=197.3269804
COULOMB=1.4399645478 # MeV fm
N13_TAU=9.951*60/np.log(2)
DATA=json.loads((ROOT/'padf2-c13.json').read_text())


def mass(n): return pyn.Nucleus(n).mass

def qvalue(left,right):
 return sum(mass(n)*v for n,v in left.items())-sum(mass(n)*v for n,v in right.items())

Q={
 'pn':qvalue({'c13':1,'p':1},{'n13':1,'n':1}),
 'pd':qvalue({'c13':1,'p':1},{'c12':1,'d':1}),
 'palpha':qvalue({'c13':1,'p':1},{'b10':1,'he4':1}),
 'n13_decay_atomic':mass('n13')-mass('c13'),
 'np_capture':qvalue({'n':1,'p':1},{'d':1}),
 'c12_capture':qvalue({'c12':1,'n':1},{'c13':1}),
 'DD_Tp':qvalue({'d':2},{'t':1,'p':1}),
 'DD_He3n':qvalue({'d':2},{'he3':1,'n':1}),
 'DT':qvalue({'d':1,'t':1},{'he4':1,'n':1}),
 'DHe3':qvalue({'d':1,'he3':1},{'he4':1,'p':1}),
 'net_D':qvalue({'d':2},{'he4':1}),
 'net_H':qvalue({'p':4},{'he4':1})}
Q['n13_neutrino']=.706
Q['heat_per_new_neutron']=(Q['net_H']-2*Q['n13_neutrino'])/2
# 6D -> 2He + 2p + 2n, followed by 2(n+p -> D), gives net 4D -> 2He.
# Equal DD branching, all T and He3 subsequently burned; zero incident KE.
DD_NEUTRON_PER_HE=.5*(Q['DD_He3n']*mass('he3')/(mass('he3')+mass('n'))+
                            Q['DT']*mass('he4')/(mass('he4')+mass('n')))
DD_THERMAL_PER_NEW_N=.5*(DD_NEUTRON_PER_HE+Q['np_capture'])
DD_CHARGED_PER_NEW_N=.5*Q['net_D']-DD_THERMAL_PER_NEW_N


def xs(za,energy):
 d=DATA['channels'][f'{za}_0']
 assert all(law==2 for _,law in d['interpolation'])
 return np.interp(energy,d['energy_MeV'],d['sigma_barn'],left=0,right=0)


def spectrum(mean,std):
 if std==0: return np.array([mean]),np.array([1.])
 e=np.linspace(max(.05,mean-6*std),mean+6*std,1201)
 w=np.exp(-.5*((e-mean)/std)**2); w[[0,-1]]*=.5; w/=sum(w)
 return e,w


def nuclear_and_drag(mean=4.5,std=.1,ln_lambda=25.):
 e,w=spectrum(mean,std)
 gamma=1+e*MEV/(MP*C*C); v=C*np.sqrt(1-gamma**-2)
 pn=xs(7013,e); sv=float(sum(w*pn*v)) # barns m/s
 if sv<=0: raise ValueError('No useful reaction')
 def avg(f): return float(sum(w*v*f)/sv)
 rpd=avg(xs(1002,e)); rC12=avg(xs(6012,e)); ralpha=avg(xs(5010,e))
 rpg=avg(xs(7014,e)); inel=avg(xs(6013,e))
 # Minimum inelastic excitation is 3.089 MeV. Higher states cost more.
 inelastic_loss=inel*3.089443
 # Nonrelativistic Rutherford recoil on stationary fully stripped C13.
 # E*sigma_energy = 4*pi*k^2/(M_C*v_nonrel^2)*ln Lambda.
 ratio=MC/MP; ecm=e*ratio/(1+ratio); b90=6*COULOMB/(2*ecm)
 recoil_fraction=4*ratio/(1+ratio)**2
 sigma_energy=2*pi*recoil_fraction*b90*b90*ln_lambda/100 # barn
 drag=avg(e*sigma_energy)
 # Small-angle laboratory deflection variance, integrated over impact parameter.
 # Its ratio to fractional recoil-energy transfer is M_C/m_p.
 angular_variance=avg(8*pi*(ratio/(1+ratio))**2*b90*b90*ln_lambda/100)
 # Input energy of every proton consumed in pn or pd, plus replenished drag.
 e_pn=avg(e*pn); e_pd=avg(e*xs(1002,e))
 # Two-body pn kinematics, isotropic in CM and ground-state final nucleus.
 # The carbon-13 catalyst recoil is directly recoverable; neutron KE is heat.
 # Mass-number kinematics is consistent to better than 1% at these energies.
 pn_neutron=avg(pn*((13/14)*(13/14*e+Q['pn'])+e/196))
 pn_recoil=e_pn+Q['pn']-pn_neutron
 charged=DD_CHARGED_PER_NEW_N+pn_recoil+e_pd+rpd*Q['pd']+drag
 thermal=(DD_THERMAL_PER_NEW_N+pn_neutron+Q['n13_decay_atomic']-Q['n13_neutrino']+
          (1-rpd)*Q['np_capture']+rpd*Q['c12_capture']+inelastic_loss)
 prep=e_pn+e_pd+drag+inelastic_loss
 closure=charged+thermal-prep-Q['heat_per_new_neutron']
 return {'mean_lab_MeV':mean,'rms_energy_width_MeV':std,'coulomb_log':ln_lambda,
  'pn_sigma_v_m3_s':sv*1e-28,'effective_pn_sigma_barn':sv/float(sum(w*v)),
  'pd_events_per_pn':rpd,'C12_residual_per_pn':rC12,'B10_residual_per_pn':ralpha,
  'N14_capture_per_pn':rpg,'inelastic_C13_events_per_pn':inel,
  'pd_C12_difference_per_pn':rpd-rC12,
  'beam_consumption_MeV_per_pn':e_pn+e_pd,'ion_recoil_drag_MeV_per_pn':drag,
  'minimum_inelastic_loss_MeV_per_pn':inelastic_loss,
  'small_angle_deflection_variance_rad2_per_pn':angular_variance,
  'pn_neutron_kinetic_MeV':pn_neutron,'pn_N13_recoil_MeV':pn_recoil,
  'beam_preparation_MeV_per_pn':prep,'directly_recoverable_ion_energy_MeV_per_pn':charged,
  'thermal_energy_MeV_per_pn':thermal,'energy_closure_MeV':closure,
  'neutron_regeneration_feasible_in_reduced_cycle':rpd<=1,
  'neutrons_left_for_hydrogen_capture_per_pn':1-rpd,
  'electron_neutralization_fraction_for_equal_cold_electron_and_carbon_drag':6*ME/MC,
  'limitations':['independent inclusive production data do not close exclusive branches',
   'p,d represented by inclusive deuteron production; invalid when other D channels open',
   'higher inelastic states, nuclear elastic scattering and p,p collisions omitted',
   'B10/N14 restoration handled only as a separate specified penalty',
   'all useful and secondary neutrons captured without parasitic loss',
   'all generated deuterium burned with zero additional driver cost',
   'ion recoil is classified as potentially recoverable, not demonstrated recovery',
   'nonrelativistic scattering, stationary carbon and ground-state pn kinematics']}


def energy_case(d,eta_a=.9,eta_th=.4,eta_direct=.9,side_cost=0.,other=0.):
 # eta_direct is the fraction sent directly to electricity. Its unrecovered
 # share becomes heat and may receive the thermal efficiency once only.
 g=eta_direct+(1-eta_direct)*eta_th
 side_rate=d['B10_residual_per_pn']+d['N14_capture_per_pn']
 w=(g*d['directly_recoverable_ion_energy_MeV_per_pn']+
    eta_th*d['thermal_energy_MeV_per_pn']-d['beam_preparation_MeV_per_pn']/eta_a-
    side_cost*side_rate-other)
 # side_cost is net electricity charged per restored catalyst, NOT a modeled
 # process. No hidden thermal refund or fuel credit is assigned to it.
 return {'accelerator_efficiency':eta_a,'thermal_efficiency':eta_th,
  'ion_direct_electric_recovery_fraction':eta_direct,'side_restoration_net_cost_MeV':side_cost,
  'other_net_input_MeV_per_pn':other,'net_MeV_per_pn_before_unmodeled_plant_loads':w,
  'reduced_cycle_valid':d['neutron_regeneration_feasible_in_reduced_cycle']}


def confinement(d,w,B=20.,fill=.2,specific_strength=1e6,mass_factor=5.,pipe_radius=.02):
 # Necessary cold Brillouin force-balance bound, not an equilibrium solution.
 # Carbon q/m limits total charge density: nH+6*nC <= 6 eps0 B^2/(2 M_C).
 # nH=6*nC maximizes reaction rate at fixed charge when N13 is removed promptly.
 charge_bound=6*EPS0*B*B/(2*MC); nh=charge_bound/2; nc=charge_bound/12
 rate=nh*nc*d['pn_sigma_v_m3_s']; field_u=B*B/(2*MU0)
 mass=mass_factor*field_u/specific_strength
 pnet=rate*w*MEV
 # Ideal centered uniform column potential/field, for a column of half pipe radius.
 r=pipe_radius/2; edge_field=charge_bound*E_CHARGE*r/(2*EPS0)
 edge_to_axis_potential=edge_field*r/2
 # Proton cyclotron power uses mean energy and an isotropic pitch distribution.
 E=d['mean_lab_MeV']*MEV; gam=1+E/(MP*C*C); v=C*sqrt(1-gam**-2)
 cyclotron_one=E_CHARGE**4*B*B*(2/3*v*v)*gam**2/(6*pi*EPS0*MP*MP*C**3)
 cyclotron_per_pn=nh*cyclotron_one/rate/MEV
 # If N13 remains in the same field volume through its decay, it occupies
 # charge inventory. Optimize H/C ratio while solving S=2sqrt(6R/k)+7R*tau.
 k=d['pn_sigma_v_m3_s']; a=7*N13_TAU; b=2*sqrt(6/k)
 y=2*charge_bound/(b+sqrt(b*b+4*a*charge_bound)); retained_rate=y*y
 # Homogeneous cold longitudinal two-fluid dispersion with an ideal conducting
 # cylinder's transverse fundamental k_perp=2.4048/a. Sufficient stability
 # condition within this model: drift exceeds the maximum coupling velocity.
 wp=sqrt(nh*E_CHARGE**2/(EPS0*MP)); wc=sqrt(nc*(6*E_CHARGE)**2/(EPS0*MC))
 kperp=2.4048255577/pipe_radius
 drift_threshold=(wp**(2/3)+wc**(2/3))**1.5/kperp
 return {'B_T':B,'plasma_field_volume_fraction':fill,
  'magnetic_structure_specific_strength_J_kg':specific_strength,'field_mass_multiplier':mass_factor,
  'optimistic_proton_density_m3':nh,'optimistic_carbon_density_m3':nc,
  'pn_rate_per_plasma_m3_s':rate,'net_power_per_plasma_m3_W':pnet,
  'field_energy_per_field_m3_J':field_u,'allocated_field_mass_per_field_m3_kg':mass,
  'optimistic_net_W_kg_field_mass_only':fill*pnet/mass,
  'N13_cold_decay_inventory_kg_per_plasma_m3':rate*N13_TAU*MC,
  'if_N13_retained_reaction_rate_fraction':retained_rate/rate,
  'proton_cyclotron_energy_MeV_per_pn_isotropic_pitch':cyclotron_per_pn,
  'pipe_radius_m':pipe_radius,'assumed_column_radius_m':r,
  'column_edge_electric_field_V_m':edge_field,
  'axis_to_edge_potential_V':edge_to_axis_potential,
  'proton_90deg_pitch_Larmor_radius_m':gam*MP*v/(E_CHARGE*B),
  'cold_pipe_two_stream_drift_threshold_m_s':drift_threshold,
  'proton_axial_speed_m_s':v,'cold_pipe_model_above_instability_band':v>drift_threshold,
  'limitations':['field-mass-only upper performance screen; no complete plant mass',
   'cold Brillouin bound does not establish a finite-temperature mixed-ion equilibrium',
   'radioactive N13 must be extracted promptly for stated reactant density',
   'pipe dispersion omits realistic transverse profiles, rotation, finite temperature, wall resistance and electromagnetic modes',
   'the stated column/pipe geometry must separately satisfy emittance and orbit constraints']}


def angular_recycle(d,theta,eta_a=.95,eta_direct=.95,eta_th=.4):
 # Rotational diffusion first passage from pole to a circular absorbing cap.
 # Resetting all escaped ions to the original beam is a specified control
 # scenario, not a necessary cost for a trap accepting all pitch angles.
 # t_escape = -8 ln(cos(theta/2)) / (d<theta^2>/dt).
 factor=-8*np.log(np.cos(theta/2))
 events=d['small_angle_deflection_variance_rad2_per_pn']/factor
 energy=events*d['mean_lab_MeV']
 g=eta_direct+(1-eta_direct)*eta_th
 return {'accepted_half_angle_rad':theta,'recycling_events_per_pn':events,
  'additional_recirculated_beam_energy_MeV_per_pn':energy,
  'electrical_penalty_MeV_per_pn':energy*(1/eta_a-g),
  'scope':'Small-angle rotational diffusion with absorbing acceptance cone; not all-angle confinement.'}


def main():
 cases=[]
 for E in [3.3,3.775,3.99,4.5,5.,6.,7.]:
  for sd in [0.,.1]:
   d=nuclear_and_drag(E,sd)
   d['energy_cases']=[energy_case(d,a,t,c,s) for a,t,c,s in [(.9,.4,0,0),(.9,.6,0,0),(.9,.4,.8,0),(.9,.4,.9,0),(.95,.4,.95,0),(.95,.4,.95,50),(.95,.4,.95,230)]]
   d['angular_recycling_scenarios']=[angular_recycle(d,t) for t in [.3,1.,pi/2]]
   cases.append(d)
 selected=nuclear_and_drag(4.5,.1)
 selected_w=energy_case(selected,eta_a=.95,eta_direct=.95,side_cost=50)['net_MeV_per_pn_before_unmodeled_plant_loads']
 density=[confinement(selected,selected_w,B,fill,pipe_radius=.4/B) for B in [20,50,100] for fill in [.05,.2,.5]]
 out={'scope':'Conditional C13(p,n) process and field-mass screens. No validated complete reactor.',
  'data_evaluation':DATA['evaluation'],'data_sha256':DATA['source_sha256'],'Q_MeV':Q,
  'DD_charged_MeV_per_new_neutron':DD_CHARGED_PER_NEW_N,
  'DD_thermal_MeV_per_new_neutron':DD_THERMAL_PER_NEW_N,'cases':cases,
  'selected_field_screen':{'energy_MeV':4.5,'width_MeV':.1,'net_MeV':selected_w,
   'accelerator_efficiency':.95,'direct_recovery_fraction':.95,
   'side_restoration_cost_MeV':50,'cases':density}}
 (ROOT/'nonneutral-screen.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
 for d in cases:
  if d['rms_energy_width_MeV']:
   print('E',d['mean_lab_MeV'],'pn sigma',round(d['effective_pn_sigma_barn'],4),'pd/pn',round(d['pd_events_per_pn'],3),'drag',round(d['ion_recoil_drag_MeV_per_pn'],2),'net cases',[round(x['net_MeV_per_pn_before_unmodeled_plant_loads'],3) for x in d['energy_cases']])
 print('field-only selected',[(x['B_T'],x['plasma_field_volume_fraction'],round(x['optimistic_net_W_kg_field_mass_only'],6)) for x in density])
 print('energy closure',max(abs(d['energy_closure_MeV']) for d in cases))

if __name__=='__main__': main()
