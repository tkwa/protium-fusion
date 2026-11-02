"""Finite heat-engine conductance and an explicitly conditional cold-loop screen.

A central 300 K engine feeding a distant 30 K spherical radiator is one layout,
not a universal reactor requirement. The loop uses ideal-gas He4, a prescribed
Darcy friction factor, thin pressure-tube walls and a single or parallel pair
of straight supply/return pipes. Branches, manifolds, insulation, pumps and
panel heat exchangers require additional equipment. Existing radiator and
auxiliary allowances may cover part of these requirements, so costs are not
silently charged twice. Distributed engines or other energy transport need a
separate layout and can avoid this central-cold-loop penalty.
"""
from pathlib import Path
from math import pi,sqrt
import json,hashlib
from scipy.optimize import minimize_scalar,brentq
from capture_competition import read_reference,ROOT,HERE,YEAR,MEV,SB
from nonneutral_screen import Q,DD_CHARGED_PER_NEW_N,DD_THERMAL_PER_NEW_N
RHE=8.31446261815324/.004002602
CPHE=2.5*RHE

def thermal_power(r):
 """Heat input actually credited by the source/circulation equations."""
 e=json.loads((ROOT/'network-ion-budget.json').read_text());n=e['network'];a=r['assumptions']
 terms={'network_thermal':e['thermal_MeV_per_primary'],
        'primary_charged_converter_loss':(1-a['direct'])*e['charged_MeV_per_primary']}
 for key in ['rotation_converter_loss','neutron_rotation','wall_converter_loss','potential_return_loss','end_return_loss','ion_radiation']:terms[key]=0.
 for key,p in r['selected_source_profiles'].items():
  f=n['processor_pn_rates'][key]
  terms['rotation_converter_loss']+=f*(1-a['direct'])*(p['rotation_input_MeV_per_pn']-p['neutron_rotation_MeV_per_pn'])
  terms['neutron_rotation']+=f*p['neutron_rotation_MeV_per_pn']
  terms['wall_converter_loss']+=f*(1-a['direct'])*p['wall_circulation_MeV_per_pn']
  terms['potential_return_loss']+=f*(1-a['potential_out'])*p['potential_circulation_MeV_per_pn']
  terms['end_return_loss']+=f*(1-a['end_out'])*p['end_circulation_MeV_per_pn']
  terms['ion_radiation']+=f*(p['Larmor_MeV_per_pn']+p['ion_brems_MeV_per_pn'])
 light=n['net_neutrons_available_for_light_fuels']-n['light_net_before_D_completion'].get('2003',0.)
 terms['failed_H_capture']=-light*(1-a['neutron_capture'])*((1-a['direct'])*DD_CHARGED_PER_NEW_N+DD_THERMAL_PER_NEW_N+Q['np_capture'])
 area=2*pi*(a['core_radius']+a['shield'])*r['length_m']+2*pi*(a['core_radius']+a['shield'])**2
 terms['shield_and_thermal_cold_leak']=-(r['cold_radiation_W']+a['warm_to_cold_W_m2']*area)/(r['primary_pn_rate_s']*MEV)
 H=sum(terms.values())*r['primary_pn_rate_s']*MEV
 return {'terms_MeV_per_primary':terms,'hot_input_W':H,'electric_generation_W':a['heat']*H,'engine_cold_rejection_W':(1-a['heat'])*H}

def exchanger(Qh,eta=.7,relative_Carnot=1.,Th=300.,Tr=30.,engine_hot_K=None):
 # eta = relative_Carnot*(1-Tc/Tx). Phi is a scenario, not measured.
 ratio=1-eta/relative_Carnot
 if ratio<=Tr/Th:return {'feasible':False,'reason':'No finite temperature approaches'}
 b=sqrt((1-eta)*ratio)
 Tx=(Tr+b*Th)/(ratio+b) if engine_hot_K is None else engine_hot_K
 if not Tr/ratio<Tx<Th:raise ValueError('engine temperature outside finite-conductance interval')
 Tc=ratio*Tx
 Kh=Qh/(Th-Tx);Kc=(1-eta)*Qh/(Tc-Tr)
 # This choice minimizes Kh+Kc at fixed eta and relative-Carnot efficiency.
 return {'feasible':True,'relative_Carnot':relative_Carnot,'eta':eta,'hot_reservoir_K':Th,'radiator_surface_K':Tr,
  'selection':'minimum total conductance' if engine_hot_K is None else 'specified engine hot temperature','engine_hot_K':Tx,'engine_cold_K':Tc,'hot_approach_K':Th-Tx,'cold_approach_K':Tc-Tr,
  'hot_conductance_W_K':Kh,'cold_conductance_W_K':Kc,'total_conductance_W_K':Kh+Kc,
  'engine_entropy_generation_W_K':(1-eta)*Qh/Tc-Qh/Tx}

def loop(Qcold,temperature_span_K,mean_K,path_length_m,velocity_m_s,pressure_Pa=1e5,
         darcy_f=.005,parallel_pipes=1,pump_efficiency=.8,wall_density=1600.,wall_stress_Pa=1e9,gas_R=RHE,gas_cp=CPHE):
 # A complete out-and-back length enters both fluid inventory and friction.
 rho=pressure_Pa/(gas_R*mean_K);h=gas_cp*temperature_span_K;v=velocity_m_s;L=path_length_m
 # For each of N identical pipes D = sqrt(4 mdot/(N pi rho v)).
 # Darcy pressure work / kg = f L v^2/(2D) = K/sqrt(mdot).
 K=darcy_f*L*v*v/2*sqrt(parallel_pipes*pi*rho*v/4)
 # Fluid heat balance: mdot*h = Qcold + mdot*pressure_work_per_kg.
 y=(K+sqrt(K*K+4*h*Qcold))/(2*h);mdot=y*y
 A=mdot/(rho*v);D=sqrt(4*A/(parallel_pipes*pi))
 specific_friction=darcy_f*L/D*v*v/2
 fluid_dissipation=mdot*specific_friction;pump=fluid_dissipation/pump_efficiency
 fluid_mass=mdot*L/v
 # Thin circular tube: M_wall/L = 2 rho_wall P A/sigma.
 wall_mass=2*wall_density*pressure_Pa/wall_stress_Pa*A*L
 pressure_drop=rho*specific_friction
 # A low fractional pressure drop supports the constant-density approximation.
 return {'Qcold_W':Qcold,'gas_R_J_kg_K':gas_R,'gas_cp_J_kg_K':gas_cp,'temperature_span_K':temperature_span_K,'mean_K':mean_K,'path_length_m':L,
  'velocity_m_s':v,'pressure_Pa':pressure_Pa,'darcy_f':darcy_f,'parallel_pipes':parallel_pipes,
  'pump_efficiency':pump_efficiency,'wall_density_kg_m3':wall_density,'wall_stress_Pa':wall_stress_Pa,
  'density_kg_m3':rho,'mass_flow_kg_s':mdot,'total_flow_area_m2':A,'diameter_per_pipe_m':D,
  'wall_thickness_m':pressure_Pa*D/(2*wall_stress_Pa),'fluid_mass_kg':fluid_mass,'pressure_wall_mass_kg':wall_mass,
  'pipe_and_fluid_mass_kg':fluid_mass+wall_mass,'fluid_dissipation_W':fluid_dissipation,'pump_electric_W':pump,
  'pressure_drop_Pa':pressure_drop,'pressure_drop_over_absolute_pressure':pressure_drop/pressure_Pa,
  'Mach_number':v/sqrt(gas_cp/(gas_cp-gas_R)*gas_R*mean_K),'constant_density_screen_pass':bool(pressure_drop/pressure_Pa<.1),
  'fluid_heat_relative_residual':(mdot*h-Qcold-fluid_dissipation)/Qcold}

def optimize_loop(r,H,x,span_fraction=.5,pressure_Pa=1e5,darcy_f=.005,parallel_pipes=1,distance_factor=1.):
 radius=sqrt(r['radiator_area_m2']/(4*pi));L=2*radius*distance_factor
 span=x['cold_approach_K']*span_fraction;mean=r['assumptions']['radiator_K']+x['cold_approach_K']/2
 Qcold=H['engine_cold_rejection_W'];k=r['assumptions']['annual_replacement_fraction']*r['assumptions']['embodied_J_kg']/YEAR
 # Minimize the marginal demand against the user's 1mW/kg target. This picks
 # the velocity balancing pump electricity and added mass/replacement burden.
 def candidate(logv):
  return loop(Qcold,span,mean,L,10**logv,pressure_Pa,darcy_f,parallel_pipes)
 def score(logv):
  z=candidate(logv)
  return z['pump_electric_W']+(.001+k)*z['pipe_and_fluid_mass_kg']
 opt=minimize_scalar(score,bounds=(-2,2),method='bounded',options={'xatol':1e-11})
 z=candidate(opt.x)
 z.update({'span_fraction':span_fraction,'distance_factor':distance_factor,'radiator_sphere_radius_m':radius,
  'relative_Carnot':x['relative_Carnot'],'objective_W_equivalent_at_target':score(opt.x),
  'mass_fraction_of_allocated_radiators':z['pipe_and_fluid_mass_kg']/r['mass_breakdown_kg']['radiators'],
  'pump_fraction_of_allocated_auxiliaries':z['pump_electric_W']/r['other_auxiliary_W'],
  'mass_fraction_of_total_plant':z['pipe_and_fluid_mass_kg']/r['total_mass_kg']})
 return z

def additional_cost_case(r,z,radiator_mass_fraction_available=0.,auxiliary_power_fraction_available=0.):
 """Optional budget scenario: only uncovered portions are charged as additions."""
 a=r['assumptions'];dm=max(0,z['pipe_and_fluid_mass_kg']-radiator_mass_fraction_available*r['mass_breakdown_kg']['radiators'])
 dp=max(0,z['pump_electric_W']-auxiliary_power_fraction_available*r['other_auxiliary_W'])
 oldmass=r['mass_breakdown_kg'];fixed=sum(v for k,v in oldmass.items() if k not in ['radiators','replacement_manufacturing_equipment','replacement_inventory'])+dm
 bus=r['plant_bus_net_W']-dp;Qn=r['nuclear_after_capture_W']
 radkg=a['radiator_areal_mass']/(a['radiator_emissivity']*SB*(a['radiator_K']**4-a['ambient_K']**4))
 k=a['annual_replacement_fraction']*a['embodied_J_kg']/YEAR;spare=a['annual_replacement_fraction']*a['replacement_inventory_years']
 M=(fixed+(Qn-bus)*radkg)/(1-k*(radkg+a['manufacturing_kg_W'])-spare);P=bus-k*M
 return {'radiator_mass_fraction_available':radiator_mass_fraction_available,'auxiliary_power_fraction_available':auxiliary_power_fraction_available,
  'extra_mass_kg':dm,'extra_electric_load_W':dp,'total_mass_kg':M,'net_W':P,'net_W_kg':P/M,
  'scope':'Conditional accounting only: pump and pipe deficits are charged once beyond stated available shares. Added pump/exchanger machinery, geometry updates, and altered heat-engine operation are not designed.'}

def main():
 r=read_reference();H=thermal_power(r);xs=[exchanger(H['hot_input_W'],relative_Carnot=p) for p in [.777,.78,.8,.85,.9,1.]]
 loops=[]
 for x in xs:
  if not x['feasible']:continue
  for span in [1.,.5]:
   z=optimize_loop(r,H,x,span);z['additional_cost_scenarios']=[additional_cost_case(r,z,*f) for f in [(0.,0.),(.5,.5),(1.,1.)]];loops.append(z)
 alternatives=[];envelopes=[]
 for phi in [.8,.85]:
  for hot_gap in [1.,5.,10.]:
   ex=exchanger(H['hot_input_W'],relative_Carnot=phi,engine_hot_K=300-hot_gap)
   z=optimize_loop(r,H,ex,.5);z['exchanger']=ex;alternatives.append(z)
  xmax={'relative_Carnot':phi,'cold_approach_K':(1-.7/phi)*300-30}
  for span in [1.,.5]:
   z=optimize_loop(r,H,xmax,span);z['scope']='Limiting internal hot temperature approaches300K, so hot exchanger conductance diverges. Optimistic cold-loop envelope within this central ideal-He/Darcy model, not an achievable complete system or universal bound.';envelopes.append(z)
 sensitivity=[];x=next(x for x in xs if x.get('relative_Carnot')==.8)
 for f in [.0025,.005,.01,.02]:sensitivity.append(optimize_loop(r,H,x,darcy_f=f))
 for n in [16,256]:sensitivity.append(optimize_loop(r,H,x,parallel_pipes=n))
 for d in [.1,.01]:sensitivity.append(optimize_loop(r,H,x,distance_factor=d))
 for p in [1e4,1e6]:sensitivity.append(optimize_loop(r,H,x,pressure_Pa=p))
 out={'scope':__doc__,'thermal_engine':H,'gas_properties':{'specific_gas_constant_J_kg_K':RHE,'ideal_monatomic_cp_J_kg_K':CPHE},
  'reference':{'net_W':r['net_W'],'mass_kg':r['total_mass_kg'],'radiator_area_m2':r['radiator_area_m2'],
   'radiator_mass_kg':r['mass_breakdown_kg']['radiators'],'other_auxiliary_W':r['other_auxiliary_W']},
  'heat_exchangers':xs,'central_loops':loops,'alternative_exchanger_loops':alternatives,'optimistic_cold_loop_envelopes':envelopes,'sensitivities':sensitivity,
  'sources':[{'url':'https://pmc.ncbi.nlm.nih.gov/articles/PMC13649922/','use':'2026 primary He viscosity/EOS study, Table9: 0.1MPa,25K density1.9246kg/m3; sanity comparison for ideal-gas assumption only. No viscosity is imported; Darcy friction remains a scenario.'},
             {'url':'https://doi.org/10.1680/ijoti.1939.13150','use':'Historical pipe-friction reference. The current screen prescribes Darcy f rather than claiming a measured value for kilometre-diameter pipes.'}],
  'input_sha256':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in ['network-plant.json','network-ion-budget.json','network_plant.py','nonneutral_screen.py']}}
 (HERE/'heat-transport.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
 print('heat engine',H)
 for z in loops:
  print('phi',z['relative_Carnot'],'span',z['span_fraction'],'v',z['velocity_m_s'],'diameter',z['diameter_per_pipe_m'],'mass',z['pipe_and_fluid_mass_kg'],'pump',z['pump_electric_W'],'radfraction',z['mass_fraction_of_allocated_radiators'],'auxfraction',z['pump_fraction_of_allocated_auxiliaries'],'dp/p',z['pressure_drop_over_absolute_pressure'])
if __name__=='__main__':main()
