#!/usr/bin/env python3
"""Pressure/barrier requirements; no assumed connected-hole production yield."""
import json,math
from pathlib import Path
KB=1.380649e-23; MU=1.66053906660e-27; YEAR=365.25*86400
M={'H':1.00782503223*MU,'H2':2*1.00782503223*MU,'He':4.00260325413*MU}
ROOT=Path(__file__).resolve().parent

def speed(T,m):return math.sqrt(8*KB*T/(math.pi*m))
def wall_flux(P,T,m):return P/math.sqrt(2*math.pi*m*KB*T)
def aperture_speed(area,T,m):return area*speed(T,m)/4

def fraction_leaked(phi,beta):
 assert phi>=0 and beta>=0 and phi+beta>0
 return phi/(phi+beta)

def pressure_decay(P0,h,phi,beta,T,m,t):
 """Isothermal compartment, two molecular outlets at zero pressure."""
 k=aperture_speed(phi+beta,T,m)/h
 return P0*math.exp(-k*t)

def hold_phi_max(f,h,t,T,m):return -math.log1p(-f)*4*h/(speed(T,m)*t)
def serial_effective(C1,C2,pump=0):return C1*C2/(C1+C2+pump)

def sqrt_pressure_leak_fraction(z):
 # Pump throughput aP, membrane b*sqrt(P), z=a*sqrt(P0)/b.
 assert z>0
 if z<1e-4:return 1-2*z/3+z*z/2-2*z**3/5
 return 2*(z-math.log1p(z))/z**2

def evaluate(snap):
 cand=snap['candidate'];T=cand['temperature_K']; f=cand['conditional_front_fraction'];mouth=cand['slit_mouth_fraction'];channel_beta=cand['channel_beta_reference'];pump=cand['downstream_pump_speed_m_s_reference'];S_channel={k:aperture_speed(channel_beta,T,m) for k,m in M.items()};S={k:1/(1/v+1/pump) for k,v in S_channel.items()};beta_species={k:S[k]/(speed(T,m)/4) for k,m in M.items()};beta=min(beta_species[k] for k in ['H2','He'])
 flowH=sum(x['H_atoms_m2']/x['period_s'] for x in snap['stages'].values());flowHe=sum(x['He_atoms_m2']/x['period_s'] for x in snap['stages'].values());rates={'H':flowH,'H2':flowH/2,'He':flowHe}
 fixed=[]
 for name,m in M.items():
  for P in [.001,1.,1000.,100000.]:
   J=wall_flux(P,T,m);phi=f*rates[name]/J
   fixed.append({'species':name,'partial_P_Pa':P,'thermal_particles_m2_s':J,'equivalent_phi_max_at_mean_source':phi,'per_cap_area_phi_max':phi/mouth,'one_nm_radius_zero_thickness_holes_per_receiver_m2':phi/(math.pi*1e-18),'scope':'maintained partial pressure; leakage budget f times mean source for this species'})
 phi=f/(1-f)*beta
 continuous={'beta_min_H2_He':beta,'channel_beta':channel_beta,'downstream_pump_speed_m_s':pump,'beta_species':beta_species,'phi_max_by_species':{k:f/(1-f)*v for k,v in beta_species.items()},'S_per_wall_m2_m_s':S,'phi_max':phi,'per_cap_area_transmission_max':phi/mouth,'mean_partial_P_Pa':{k:rates[k]*KB*T/S[k] for k in M},'H2_He_mean_total_P_Pa':rates['H2']*KB*T/S['H2']+rates['He']*KB*T/S['He'],'assumption':'geometry channel+via conductance in series with stated downstream pump speed; finite ultimate pressure and additional plenum resistance remain to specify; no source-retention fraction credited'}
 stores=[];stages={};carbon_per_cap_area=2000*cand['cap_m']/(12*MU)
 for name,x in snap['stages'].items():
  g=x['H_atoms_m2']/2+x['He_atoms_m2'];PV=g*KB*T
  stages[name]={**x,'full_H2_He_pulse_PV_per_wall_Pa_m':PV,'P_full_pulse_local_post_corrected_plenum_Pa':PV/cand['local_plenum_void_depth_m'],'P_full_pulse_all_void_isothermal_sensitivity_Pa':PV/cand['connected_void_depth_if_all_isothermal_m'],'NRT25_cap_dpa_proxy_per_pulse':x['cap_0_to_100nm_NRT25_events_m2']/carbon_per_cap_area,'NRT40_cap_dpa_proxy_per_pulse':x['cap_0_to_100nm_NRT40_events_m2']/carbon_per_cap_area,'upper_peak_P_if_immediate_H2_He_release_Pa':KB*T*(x['H_peak_flux_m2_s']/2/S['H2']+x['He_peak_flux_m2_s']/S['He'])}
  for h in cand['plenum_depth_m']:
   stores.append({'stage':name,'plenum_h_m':h,'full_pulse_pressure_Pa':PV/h,'H2_emptying_time_s':h/S['H2'],'He_emptying_time_s':h/S['He']})
 holds=[{'species':k,'h_m':h,'hold_s':t,'phi_max':hold_phi_max(f,h,t,T,m)}for k,m in M.items() for h in cand['plenum_depth_m'] for t in [.01,1.,100.,7698.198541942018]]
 ionrate=sum(x['charged_ion_atoms_m2']/x['period_s'] for x in snap['stages'].values());nrtrate=sum(x['cap_0_to_100nm_NRT25_events_m2']/x['period_s'] for x in snap['stages'].values())
 damage=[]
 for radius in [.5e-9,1e-9,5e-9]:
  hole=math.pi*radius**2
  for name,x in snap['stages'].items():
   damage.append({'stage':name,'hole_radius_m':radius,'hole_open_area_m2':hole,'max_y_per_ion_for_whole_pulse_persistent_holes':phi/(mouth*x['charged_ion_atoms_m2']*hole),'mean_y_times_heal_time_s_max':phi/(mouth*ionrate*hole),'peak_y_times_heal_time_s_max':phi/(mouth*x['charged_peak_flux_m2_s']*hole),'long_cylindrical_pore_transmission_scenario':min(1,8*radius/(3*cand['cap_m'])),'scope':'unit-transmission equivalent holes; y unknown, not inferred from NRT; mean condition is insufficient during bursts'})
 cadence=sum(1/x['period_s'] for x in snap['stages'].values());masscap=2000*cand['cap_m']
 return {'scope':'analytical conditional requirements; f=.001 is a comparison point, not a universal plant acceptance criterion','T_K':T,'f':f,'fixed_partial_pressure':fixed,'continuous_extraction':continuous,'full_pulse_storage':stores,'sealed_hold':holds,'stages':stages,'damage_conditions':damage,'mean_charged_ion_flux_m2_s':ionrate,'cap_NRT25_dpa_proxy_per_year':nrtrate*YEAR/carbon_per_cap_area,'refurbishment_sensitivities':{'mouth_caps_kg_m2_per_event':masscap*mouth,'whole_100nm_face_kg_m2_per_event':masscap,'mouth_caps_W_m2_per_MJ_kg_input_every_event':masscap*mouth*cadence*1e6,'whole_100nm_face_W_m2_per_MJ_kg_input_every_event':masscap*cadence*1e6},'nonlinear_H_pressure_law':{'scope':'independent mathematical pressure-law sensitivity; no hot-laminate permeability assigned','single_cohort_fraction_lost_for_pump_to_sqrt_leak_initial_ratio':{str(z):sqrt_pressure_leak_fraction(z) for z in [10.,100.,1000.,2000.,10000.]}},'two_barriers':{'formula':'Ceff=C1*C2/(C1+C2+S_interspace); coefficients per receiver area, equal T','equal_unpumped_resistance_factor':.5},'checks':{'continuous_fraction_identity':abs(fraction_leaked(phi,beta)-f)<1e-15,'mass_convention':'neutral atomic masses for thermal velocities; nuclei counts for inventory'}}

if __name__=='__main__':
 d=evaluate(json.loads((ROOT/'input-snapshot.json').read_text()));(ROOT/'results.json').write_text(json.dumps(d,indent=2)+'\n')
 print('continuous pressure',d['continuous_extraction']['H2_He_mean_total_P_Pa'],'phi',d['continuous_extraction']['phi_max'],'annual capNRT',d['cap_NRT25_dpa_proxy_per_year'])
