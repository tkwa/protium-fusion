"""Shared-solenoid geometry and conditional complete-mass allocation.

The winding is a distributed hoop-supported continuous current layer. Local
radial stress is set to zero; conductor and circumferential structure fill the
cross section. Axial compression, global Euler and local shell buckling are
separate checks. This force-balanced idealization is not an engineered coil.
Cryogenic/shield, extraction, regeneration and D-plant allowances remain explicit.
"""
from pathlib import Path
from functools import lru_cache
from math import pi,sqrt
import json,argparse
import numpy as np
from scipy.integrate import solve_ivp
from catalyst_penning import equilibrium,at_wall,MP
from catalyst_thermal import case,CAT,MEV,C,nucmass
from nonneutral_screen import MU0,EPS0,Q,DD_CHARGED_PER_NEW_N,DD_THERMAL_PER_NEW_N
from ion_radiation_and_ends import ion_brems
ROOT=Path(__file__).parent
SB=5.670374419e-8
YEAR=365.25*86400


def winding(B=100,a=105,J=1e8,sigma_hoop=8e9,rho_hoop=1600,rho_conductor=8000,rtol=1e-10):
 # t decreases from B to 0. dr/dB is negative; all accumulated masses/energy
 # increase as the integration advances toward zero field.
 def rhs(b,y):
  r=y[0];jav=1/(1/J+b*r/sigma_hoop);fc=jav/J;fs=jav*b*r/sigma_hoop
  dr=-1/(MU0*jav)
  return [dr,2*pi*r*rho_conductor*fc*dr,2*pi*r*rho_hoop*fs*dr,
          2*pi*r*b*b/(2*MU0)*dr]
 ucore=B*B/(2*MU0)*pi*a*a
 sol=solve_ivp(rhs,[B,0],[a,0,0,ucore],rtol=rtol,atol=rtol*1e-3,dense_output=True,max_step=B/200)
 bs=np.linspace(B,0,1001);ys=sol.sol(bs)
 r,mc,mh,u=map(float,sol.y[:,-1])
 assert sol.success
 return {'B_T':B,'inner_radius_m':a,'outer_radius_m':r,'operating_engineering_J_A_m2':J,
  'hoop_working_stress_Pa':sigma_hoop,'hoop_density_kg_m3':rho_hoop,
  'conductor_density_kg_m3':rho_conductor,'conductor_mass_kg_per_m':mc,
  'hoop_mass_kg_per_m':mh,'field_energy_J_per_m':u,'bore_field_energy_J_per_m':ucore,
  'hoop_virial_relative_error':mh*sigma_hoop/rho_hoop/(2*u)-1,
  '_B':bs,'_r':ys[0]}


def finite_axis_factor(w,L,z):
 # Integrate finite-length on-axis contributions of every cylindrical current
 # sheet. This evaluates the uniform-long-cylinder winding, without claiming
 # that finite end forces or off-axis fields remain identical.
 r=w['_r'];b=w['_B'];zp=z+L/2;zm=z-L/2
 f=.5*(zp/np.sqrt(zp*zp+r*r)-zm/np.sqrt(zm*zm+r*r))
 return float(np.trapezoid(f,x=-b)/w['B_T'])


@lru_cache(maxsize=8)
def source_equilibrium(profile='refined'):
 name='titanium-source-refined.json' if profile=='refined' else 'catalyst-penning-49Ti.json'
 b=json.loads((ROOT/name).read_text())['best'][0]
 z,a=CAT['49Ti']['Z'],CAT['49Ti']['A']
 eq=equilibrium(z,a,b['h0'],b['c0'],b['xi'],nout=6001,rtol=1e-10)
 assert eq['success']
 return eq


@lru_cache(maxsize=2048)
def source(T=1.7,B=100,heat=.6,profile='refined',side_cost=50.,heater=.95,direct=.95,potential_in=.999,potential_out=.999):
 from scipy.optimize import minimize_scalar
 eq=source_equilibrium(profile)
 fn=lambda x:-at_wall(eq,x,'49Ti',T=T,B=B,heat=heat,side_cost=side_cost,heater=heater,direct=direct,potential_in=potential_in,potential_out=potential_out)['field_mass_only_net_W_kg']
 opt=minimize_scalar(fn,bounds=(15,40),method='bounded',options={'xatol':1e-4})
 row=at_wall(eq,opt.x,'49Ti',T=T,B=B,heat=heat,side_cost=side_cost,heater=heater,direct=direct,potential_in=potential_in,potential_out=potential_out)
 row['profile_family']=profile
 return row


def plant(core_radius=100,aspect=30,B=100,T=1.7,heat=.6,shield=5,packing=.7,
 J=1e8,sigma_hoop=8e9,rho_hoop=1600,sigma_comp=4e9,rho_comp=3500,Ecomp=1e12,
 comp_mass_factor=1.3,coil_cold_K=4.2,radiator_K=30,hot_K=300,
 radiator_areal_mass=1.,radiator_emissivity=.9,ambient_K=3.,
 cryo_Carnot_fraction=.3,neutron_leak=1e-6,gamma_leak=1e-5,
 warm_to_cold_W_m2=.05,outer_to_cold_W_m2=.001,neutron_capture=.995,
 J_fraction_of_critical=.1,power_law_n=30,Ec_V_m=1e-4,
 extraction_MeV_per_pn=.2,other_fraction=.1,
 conversion_kg_W=.1,Dplant_kg_W=.2,processing_kg_W=.1,
 warm_water_density=1000,wall_E=1e12,wall_density=3500,water_pressure=1e4,
 wall_buckling_safety=3,end_guard_radii=3,end_system_fraction=.15,
 euler_safety=2,stress_poisson=.2,annual_replacement_fraction=1e-4,embodied_J_kg=1e8,
 source_profile='refined',side_cost=50.,finite_field_derating=True,
 heater=.95,direct=.95,potential_in=.999,potential_out=.999,
 manufacturing_kg_W=1.,replacement_inventory_years=10.,
 end_barrier_kT=0.,end_transfer_in=.9999,end_transfer_out=.9999):
 L=aspect*core_radius;a=core_radius+shield
 w=winding(B,a,J,sigma_hoop,rho_hoop)
 activeL=max(0,L-2*end_guard_radii*w['outer_radius_m'])
 finite_min=finite_axis_factor(w,L,activeL/2) if activeL else finite_axis_factor(w,L,0)
 Bsource=B*finite_min if finite_field_derating else B
 src=source(T,Bsource,heat,source_profile,side_cost,heater,direct,potential_in,potential_out);nc=case('49Ti',T)
 rr=src['radius_m'];channels=packing*(core_radius/rr)**2
 rate=channels*activeL*src['reaction_rate_per_m_s']
 base=rate*src['net_MeV_per_pn']*MEV
 ion_brem_MeV=ion_brems(T=T)['radiated_MeV_per_pn']
 ion_brem_W=rate*ion_brem_MeV*MEV
 ion_brem_penalty=ion_brem_W*(1/heater-heat)
 g=direct+(1-direct)*heat
 lost_D_yield=(g*DD_CHARGED_PER_NEW_N+heat*(DD_THERMAL_PER_NEW_N+Q['np_capture']))
 capture_load=rate*(1-neutron_capture)*lost_D_yield*MEV
 # Heating/kinetic terms already charged in the trap model are not charged
 # again here. Product selection/end losses are a new inclusive parameter.
 extraction=rate*extraction_MeV_per_pn*MEV
 # Additional coherent energy circulation if incoming charge23 and V23
 # traverse the same high end barrier. Default0 retains the unspecified
 # selective-extraction hypothesis; the explicit25kT branch is separate.
 end_circulation_MeV=23*T*end_barrier_kT
 end_loss_MeV=end_circulation_MeV*(1/end_transfer_in-end_transfer_out-(1-end_transfer_out)*heat)
 end_transfer_load=rate*end_loss_MeV*MEV
 neutron_rotation=src['pair_rotation_MeV']*nucmass(0,1)/(nucmass(1,1)+nucmass(22,49))
 qn=nc['neutron_kinetic_MeV']+neutron_rotation
 qg=Q['np_capture']*neutron_capture+nc['pn_cascade_gamma_MeV']+nc['inelastic_events_per_pn']*nc['inelastic_excitation_MeV']
 radiation_cold=rate*(neutron_leak*qn+gamma_leak*qg)*MEV
 inner_area=2*pi*a*L+2*pi*a*a
 outer_area=2*pi*w['outer_radius_m']*L+2*pi*w['outer_radius_m']**2
 inner_leak=warm_to_cold_W_m2*inner_area;outer_leak=outer_to_cold_W_m2*outer_area
 vc=w['conductor_mass_kg_per_m']/8000*L
 Ec=Ec_V_m*J_fraction_of_critical**power_law_n
 creep=Ec*J*vc
 cold_load=radiation_cold+inner_leak+outer_leak+creep
 fridge=cold_load*(radiator_K/coil_cold_K-1)/cryo_Carnot_fraction
 heat_credit_lost=heat*(radiation_cold+inner_leak)
 preaux=base-capture_load-extraction-heat_credit_lost-fridge-creep-ion_brem_penalty-end_transfer_load
 other=other_fraction*max(0,preaux)
 # Distributed hoop support carries all circumferential Lorentz force.
 # Long-solenoid midplane axial Maxwell force equals field energy per length.
 F=w['field_energy_J_per_m'];area_comp=F/sigma_comp
 rcomp=sqrt(w['outer_radius_m']**2+area_comp/pi)
 IoverA=(rcomp*rcomp+w['outer_radius_m']**2)/4
 sigma_euler=pi*pi*Ecomp*IoverA/(L*L)
 tcomp=rcomp-w['outer_radius_m']
 sigma_shell=Ecomp*tcomp/(sqrt(3*(1-stress_poisson**2))*rcomp)
 comp=area_comp*rho_comp*L*comp_mass_factor
 hoop=w['hoop_mass_kg_per_m']*L;conductor=w['conductor_mass_kg_per_m']*L
 # Long cylindrical channel wall, external-pressure ovalization estimate.
 electric_wall_pressure=EPS0*src['wall_E_V_m']**2/2
 total_wall_pressure=water_pressure+electric_wall_pressure
 wall_t=rr*(4*(1-stress_poisson**2)*total_wall_pressure*wall_buckling_safety/wall_E)**(1/3)
 wall_volume=channels*pi*((rr+wall_t)**2-rr**2)*activeL
 core_vol=pi*core_radius**2*activeL
 water_vol=max(0,core_vol*(1-packing)-wall_volume)
 water_mass=warm_water_density*(water_vol+pi*(a*a-core_radius**2)*L+2*pi*a*a*shield)
 wall_mass=wall_volume*wall_density
 decay=rate*nc['decay_inventory_kg_per_pn_s']
 active_ions=channels*activeL*(src['proton_inventory_per_m']*nucmass(1,1)+src['catalyst_inventory_per_m']*nucmass(22,49))*MEV/C**2
 heat_source=rate*nc['nuclear_heat_MeV']*MEV
 nuclear_after_capture=heat_source-rate*(1-neutron_capture)*(DD_CHARGED_PER_NEW_N+DD_THERMAL_PER_NEW_N+Q['np_capture'])*MEV
 Dgross=rate*neutron_capture*(DD_CHARGED_PER_NEW_N+DD_THERMAL_PER_NEW_N)*MEV
 gross_energy_processed=rate*(nc['input_ion_energy_MeV']+src['pair_rotation_MeV']+src['potential_work_MeV_per_pn']+nc['nuclear_heat_MeV']+end_circulation_MeV)*MEV
 # Converter mass is rated on circulated energy, including potential handling.
 conversion=conversion_kg_W*gross_energy_processed
 Dmass=Dplant_kg_W*Dgross;process=processing_kg_W*max(0,heat_source)
 end_mass=end_system_fraction*(hoop+conductor+comp)
 fixed_mass=hoop+conductor+comp+water_mass+wall_mass+decay+active_ions+conversion+Dmass+process+end_mass
 # Close the maintenance electricity and radiator-mass feedback analytically.
 rad_flux=radiator_emissivity*SB*(radiator_K**4-ambient_K**4)
 rad_kg_W=radiator_areal_mass/rad_flux
 k=annual_replacement_fraction*embodied_J_kg/YEAR
 no_maintenance=preaux-other
 bus_mass=fixed_mass+(nuclear_after_capture-no_maintenance)*rad_kg_W
 bus_specific=no_maintenance/bus_mass
 # M = fixed + radiator(Q-Pnet) + factory(kM) + spares; Pnet=P0-kM.
 spare_fraction=annual_replacement_fraction*replacement_inventory_years
 denom=1-k*(rad_kg_W+manufacturing_kg_W)-spare_fraction
 if denom<=0: raise ValueError('No finite positive maintenance/radiator closure')
 total=(fixed_mass+(nuclear_after_capture-no_maintenance)*rad_kg_W)/denom
 maintenance=k*total;net=no_maintenance-maintenance
 manufacturing_mass=manufacturing_kg_W*maintenance;spares_mass=spare_fraction*total
 waste=nuclear_after_capture-net
 radiator_mass=waste*rad_kg_W;area_rad=waste/rad_flux
 finite_min=finite_axis_factor(w,L,activeL/2) if activeL else 0
 fields={k:v for k,v in w.items() if not k.startswith('_')}
 return {'scope':'Conditional shared-solenoid plant allocation; side factory, extraction, conductor at100T, 3D stresses, transport and lifetime unvalidated.',
 'core_radius_m':core_radius,'length_m':L,'active_length_m':activeL,'active_fraction':activeL/L,
 'source':src,'winding':fields,'reaction_rate_s':rate,'channels':channels,
 'net_W':net,'total_mass_kg':total,'net_W_kg':net/total,
 'plant_bus_net_W':no_maintenance,'plant_bus_mass_kg':bus_mass,'plant_bus_net_W_kg':bus_specific,'jupiter_masses':total/1.89813e27,
 'base_trap_surplus_W':base,'capture_penalty_W':capture_load,'extraction_load_W':extraction,
 'additional_end_circulation_MeV_per_pn':end_circulation_MeV,
 'additional_end_transfer_loss_MeV_per_pn':end_loss_MeV,'additional_end_transfer_load_W':end_transfer_load,
 'cold_radiation_W':radiation_cold,'inner_heat_leak_W':inner_leak,'outer_heat_leak_W':outer_leak,
 'classical_ion_bremsstrahlung_W':ion_brem_W,'ion_bremsstrahlung_net_penalty_W':ion_brem_penalty,
 'flux_creep_W':creep,'total_cold_load_W':cold_load,'refrigeration_W':fridge,
 'heat_conversion_credit_removed_W':heat_credit_lost,'other_auxiliary_W':other,
 'replacement_manufacturing_W':maintenance,'waste_heat_W':waste,'radiator_area_m2':area_rad,
 'mass_breakdown_kg':{'circumferential_structure':hoop,'conductor':conductor,'axial_compression_support':comp,
 'water_moderator_and_shield':water_mass,'channel_walls':wall_mass,'cold_decay_inventory':decay,'active_plasma_ions':active_ions,
 'conversion_and_potential_handling':conversion,'closed_D_plant':Dmass,'isotope_processing':process,
 'end_system_allowance':end_mass,'radiators':radiator_mass,
 'replacement_manufacturing_equipment':manufacturing_mass,'replacement_inventory':spares_mass},
 'axial_force_N':F,'axial_support_area_m2':area_comp,'axial_outer_radius_m':rcomp,
 'compression_working_stress_Pa':sigma_comp,'euler_critical_stress_Pa':sigma_euler,
 'local_shell_critical_stress_Pa':sigma_shell,
 'global_buckling_margin':sigma_euler/(euler_safety*sigma_comp),
 'local_buckling_margin':sigma_shell/(euler_safety*sigma_comp),
 'buckling_checks_pass':bool(min(sigma_euler,sigma_shell)>=euler_safety*sigma_comp),
 'channel_wall_m':wall_t,'electric_wall_pressure_Pa':electric_wall_pressure,'total_wall_pressure_Pa':total_wall_pressure,'water_volume_fraction_in_active_bundle':water_vol/core_vol if core_vol else 0,
 'on_axis_field_fraction_at_active_end':finite_min,'on_axis_center_field_fraction':finite_axis_factor(w,L,0),
 'stored_field_energy_J':F*L,'field_energy_per_net_power_s':F*L/net if net else None,
 'embodied_energy_per_bus_power_years':bus_mass*embodied_J_kg/no_maintenance/YEAR if no_maintenance>0 else None,
 'heat_Carnot_fraction_used':heat/(1-radiator_K/hot_K),
 'mass_closure_kg':sum([fixed_mass,radiator_mass,manufacturing_mass,spares_mass])-total,
 'energy_closure_W':nuclear_after_capture-net-waste,
 'assumptions':{'shield_m':shield,'packing':packing,'nominal_B_T':B,'effective_source_B_T':Bsource,
 'T_MeV':T,'heat_efficiency':heat,'finite_field_derating':finite_field_derating,'source_profile':source_profile,
 'heater_efficiency':heater,'direct_conversion_efficiency':direct,
 'potential_in':potential_in,'potential_out':potential_out,'side_restoration_MeV_per_event':side_cost,'source_equilibrium_reoptimized_for_plant':False,'coil_cold_K':coil_cold_K,'radiator_K':radiator_K,
 'hot_K':hot_K,'radiator_areal_mass_kg_m2':radiator_areal_mass,'neutron_leak':neutron_leak,'gamma_leak':gamma_leak,
 'neutron_capture':neutron_capture,'J_fraction_of_critical':J_fraction_of_critical,'power_law_n':power_law_n,
 'sigma_hoop_Pa':sigma_hoop,'sigma_comp_Pa':sigma_comp,'rho_hoop_kg_m3':rho_hoop,'rho_comp_kg_m3':rho_comp,
 'compression_E_Pa':Ecomp,'compression_mass_factor':comp_mass_factor,'euler_safety':euler_safety,
 'end_system_fraction':end_system_fraction,'extraction_MeV_per_pn':extraction_MeV_per_pn,
 'other_fraction':other_fraction,'annual_replacement_fraction':annual_replacement_fraction,'embodied_J_kg':embodied_J_kg,
 'end_barrier_kT':end_barrier_kT,'end_transfer_in':end_transfer_in,'end_transfer_out':end_transfer_out,
 'manufacturing_kg_W':manufacturing_kg_W,'replacement_inventory_years':replacement_inventory_years,
 'conversion_kg_W':conversion_kg_W,'Dplant_kg_W':Dplant_kg_W,'processing_kg_W':processing_kg_W,
 'water_density_kg_m3':warm_water_density,'wall_E_Pa':wall_E,'wall_density_kg_m3':wall_density,
 'water_pressure_Pa':water_pressure,'wall_buckling_safety':wall_buckling_safety,'stress_poisson':stress_poisson,
 'end_guard_radii':end_guard_radii,'radiator_emissivity':radiator_emissivity,'ambient_K':ambient_K,
 'Ec_V_m':Ec_V_m,'warm_to_cold_W_m2':warm_to_cold_W_m2,'outer_to_cold_W_m2':outer_to_cold_W_m2,'cryo_Carnot_fraction':cryo_Carnot_fraction}}

if __name__=='__main__':
 rows=[]
 for sh in [2e9,4e9,8e9,16e9]:
  for sc in [1e9,2e9,4e9,8e9]:
   r=plant(sigma_hoop=sh,sigma_comp=sc);rows.append(r)
   print('hoop',sh/1e9,'comp',sc/1e9,'mW/kg',round(r['net_W_kg']*1e3,4),'active',round(r['active_fraction'],3),'buckling',round(r['global_buckling_margin'],3),'GW',round(r['net_W']/1e9,3),flush=True)
 (ROOT/'bundle-solenoid.json').write_text(json.dumps({'cases':rows},indent=2,allow_nan=False)+'\n')
