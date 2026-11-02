"""Rough explicit alternative to a central engine and very long cold pipes.

290->150K warm-gas spokes feed N radiator-sector engines. Each sector has its
own35K engine sink,30K radiator and a local cold-helium circuit. A finite-hot-
stream thermodynamic screen includes a2% hot-side temperature approach and a
specified85% fraction of local Carnot efficiency. Helium is ideal monatomic;
warm hydrogen and helium use CoolProp8.0.0 equations of state. These are requirements
allocations, not turbine/heat-exchanger designs. Warm pipes use a prescribed
MLI-effective emissivity; surface radiation is charged conservatively. Direct
spokes are used rather than claiming an optimized branched network.
"""
from math import pi,sqrt,log
from functools import lru_cache
from CoolProp.CoolProp import PropsSI
import CoolProp
import json,hashlib
from scipy.optimize import minimize_scalar,brentq
from capture_competition import read_reference,ROOT,HERE,YEAR,SB
from heat_transport import thermal_power,loop,additional_cost_case,RHE,CPHE
GASES={'He4':{'fluid':'Helium'},'H2':{'fluid':'Hydrogen'}}

@lru_cache(None)
def stream_properties(gas,Thi,Tlo,pressure):
 fluid=GASES.get(gas,{'fluid':gas})['fluid'];mean=(Thi+Tlo)/2
 def at(key,T):return PropsSI(key,'T',T,'P',pressure,fluid)
 dh=at('Hmass',Thi)-at('Hmass',Tlo);ds=at('Smass',Thi)-at('Smass',Tlo)
 return {'fluid':fluid,'pressure_Pa':pressure,'inlet_K':Thi,'outlet_K':Tlo,
  'enthalpy_span_J_kg':dh,'entropy_span_J_kg_K':ds,'cp_effective_J_kg_K':dh/(Thi-Tlo),
  'mean_density_kg_m3':at('Dmass',mean),'hot_density_kg_m3':at('Dmass',Thi),'cold_density_kg_m3':at('Dmass',Tlo),
  'hot_sound_speed_m_s':at('A',Thi),'cold_sound_speed_m_s':at('A',Tlo),
  'mean_viscosity_Pa_s':PropsSI('V','T',mean,'P',pressure,'Hydrogen' if fluid=='OrthoHydrogen' else fluid),
  'viscosity_scope':'normal-hydrogen approximation for ortho isomer; used only for diagnostic Reynolds number, not friction' if fluid=='OrthoHydrogen' else 'CoolProp fluid transport model'}

def finite_stream_efficiency(Thi=290.,Tlo=150.,Tc=35.,hot_fraction=.98,relative_Carnot=.85,gas='He4',pressure=1e5):
 p=stream_properties(gas,Thi,Tlo,pressure)
 return relative_Carnot*(1-Tc*p['entropy_span_J_kg_K']/(hot_fraction*p['enthalpy_span_J_kg']))

def warm_pair_loop(Q,Thi,Tlo,L,v,gas,pressure=1e5,darcy_f=.005,N=1,stress=1e9):
 # Same-bore outgoing and returning pipes each have length L/2. Each leg
 # stays near its specified bulk temperature until the remote heat exchanger.
 # A mean-density reference velocity merely parameterizes common area.
 p=stream_properties(gas,Thi,Tlo,pressure);rh=p['hot_density_kg_m3'];rl=p['cold_density_kg_m3'];rm=p['mean_density_kg_m3']
 vh=v*rm/rh;vl=v*rm/rl;Fv=(vh*vh+vl*vl)/(2*v*v);h=p['enthalpy_span_J_kg']
 K=darcy_f*L*v*v*Fv/2*sqrt(N*pi*rm*v/4)
 y=(K+sqrt(K*K+4*h*Q))/(2*h);m=y*y;A=m/(rm*v);D=sqrt(4*A/(N*pi))
 Re=rm*v*D/p['mean_viscosity_Pa_s']
 fsmooth=brentq(lambda ff:1/sqrt(ff)+2*log(2.51/(Re*sqrt(ff)),10),.0001,.1)
 w=darcy_f*L*(vh*vh+vl*vl)/(4*D);diss=m*w
 inventory=(rh+rl)/2*A*L;wall=2*1600*pressure/stress*A*L
 dp=darcy_f*L*(rh*vh*vh+rl*vl*vl)/(4*D)
 return {'Qcold_W':Q,'gas_properties':p,'temperature_span_K':Thi-Tlo,'mean_K':(Thi+Tlo)/2,
  'path_length_m':L,'velocity_m_s':v,'hot_velocity_m_s':vh,'cold_velocity_m_s':vl,'pressure_Pa':pressure,
  'darcy_f':darcy_f,'parallel_pipes':N,'pump_efficiency':.8,'wall_density_kg_m3':1600.,'wall_stress_Pa':stress,
  'density_kg_m3':rm,'mass_flow_kg_s':m,'total_flow_area_m2':A,'diameter_per_pipe_m':D,
  'wall_thickness_m':pressure*D/(2*stress),'fluid_mass_kg':inventory,'pressure_wall_mass_kg':wall,
  'pipe_and_fluid_mass_kg':inventory+wall,'fluid_dissipation_W':diss,'pump_electric_W':diss/.8,
  'pressure_drop_Pa':dp,'pressure_drop_over_absolute_pressure':dp/pressure,
  'Mach_number':max(vh/p['hot_sound_speed_m_s'],vl/p['cold_sound_speed_m_s']),
  'constant_density_screen_pass':bool(dp/pressure<.1),'fluid_heat_relative_residual':(m*h-Q-diss)/Q,
  'mean_Reynolds_number':Re,'smooth_Colebrook_Darcy_f':fsmooth,'prescribed_f_exceeds_smooth_estimate':bool(darcy_f>=fsmooth)}

def optimize_network_leg(r,Q,span,mean,L,N,gas,heat_loss_temperature_range=None,pressure=1e5,emissivity=.001,insulation_areal_mass=1.,stress=1e9,darcy_f=.005):
 a=r['assumptions'];k=a['annual_replacement_fraction']*a['embodied_J_kg']/YEAR
 def point(logv):
  z=warm_pair_loop(Q,heat_loss_temperature_range[1],heat_loss_temperature_range[0],L,10**logv,gas,pressure,darcy_f,N,stress) if heat_loss_temperature_range else loop(Q,span,mean,L,10**logv,pressure_Pa=pressure,darcy_f=darcy_f,parallel_pipes=N,wall_stress_Pa=stress)
  surface=pi*(z['diameter_per_pipe_m']+2*z['wall_thickness_m'])*N*L
  if heat_loss_temperature_range:
   low,high=heat_loss_temperature_range;meanT4=(high**4+low**4)/2
   radiation=surface*emissivity*SB*(meanT4-3**4);insmass=surface*insulation_areal_mass
  else:radiation=0.;insmass=0.
  z.update({'gas':gas,'external_pipe_surface_m2':surface,'assigned_effective_emissivity':emissivity if heat_loss_temperature_range else None,
    'radiative_heat_loss_W':radiation,'insulation_mass_kg':insmass,
    'total_allocated_transport_mass_kg':z['pipe_and_fluid_mass_kg']+insmass})
  z['new_electric_load_without_recovery_W']=z['pump_electric_W']+a['heat']*radiation
  z['target_penalty_W']=z['new_electric_load_without_recovery_W']+(.001+k)*z['total_allocated_transport_mass_kg']
  return z
 opt=minimize_scalar(lambda v:point(v)['target_penalty_W'],bounds=(-2,2),method='bounded',options={'xatol':1e-10})
 return point(opt.x)

def layout(r,N=1024,gas='He4',pressure=1e5,emissivity=.001,stress=1e9,insulation_areal_mass=1.,darcy_f=.01,cold_darcy_f=.005,hot_inlet=290.,hot_outlet=150.,engine_cold=35.,relative_Carnot=.85,hot_fraction=.98):
 H=thermal_power(r);Q=H['hot_input_W'];Arad=r['radiator_area_m2'];sphereR=sqrt(Arad/(4*pi))
 warm=optimize_network_leg(r,Q,hot_inlet-hot_outlet,(hot_inlet+hot_outlet)/2,2*sphereR,N,gas,(hot_outlet,hot_inlet),pressure,emissivity,insulation_areal_mass,stress=stress,darcy_f=darcy_f)
 # Each radiator sector is represented by a circular emitting patch. Outward
 # and return legs each span one sector radius: conservative for radial trees,
 # optimistic about manifolds/fittings/coverage. All sectors operate in parallel.
 sectorR=sqrt(Arad/(pi*N));cold=optimize_network_leg(r,.3*Q,(engine_cold-30)/2,(engine_cold+30)/2,2*sectorR,N,'He4',None,pressure,stress=stress,darcy_f=cold_darcy_f)
 eta=finite_stream_efficiency(hot_inlet,hot_outlet,engine_cold,hot_fraction,relative_Carnot,gas,pressure);rad=warm['radiative_heat_loss_W'];pump=warm['pump_electric_W']+cold['pump_electric_W']
 z={'pipe_and_fluid_mass_kg':warm['total_allocated_transport_mass_kg']+cold['total_allocated_transport_mass_kg'],
    'pump_electric_W':pump+.7*rad}
 accounting_reference=dict(r);accounting_reference['nuclear_after_capture_W']-=rad
 # Escaped warm-pipe radiation has already left the system and does not also
 # load the 30K panels. Auxiliary capacity can be explicitly assigned to the
 # combined pump and lost-generation burden.
 return {'hot_inlet_K':hot_inlet,'hot_outlet_K':hot_outlet,'engine_cold_K':engine_cold,'relative_Carnot':relative_Carnot,'hot_fraction':hot_fraction,'sector_count':N,'gas':gas,'pressure_Pa':pressure,'effective_warm_emissivity':emissivity,'pipe_working_stress_Pa':stress,'warm_insulation_kg_m2':insulation_areal_mass,'darcy_f':darcy_f,'cold_darcy_f':cold_darcy_f,
  'radiator_sphere_radius_m':sphereR,'sector_equivalent_radius_m':sectorR,'sector_heat_input_W':Q/N,
  'warm_spokes':warm,'local_cold_loops':cold,'finite_stream_efficiency_before_transport':eta,
  'hot_exchanger_conductance_W_K':Q*warm['gas_properties']['entropy_span_J_kg_K']/((1-hot_fraction)*warm['gas_properties']['enthalpy_span_J_kg']),
  'cold_exchanger_conductance_W_K':.3*Q/(engine_cold-30),
  'warm_radiation_fraction_of_source_heat':rad/Q,
  'source_to_generation_fraction_after_radiative_loss':eta*(1-rad/Q),
  'source_to_generation_after_pumps_no_pump_heat_credit':(eta*(Q-rad)-pump)/Q,
  'total_transport_mass_kg':z['pipe_and_fluid_mass_kg'],'total_pump_electric_W':pump,
  'lost_reference_generation_from_warm_radiation_W':.7*rad,
  'transport_mass_over_allocated_radiator_mass':z['pipe_and_fluid_mass_kg']/r['mass_breakdown_kg']['radiators'],
  'pump_and_radiative_load_over_auxiliary_allowance':z['pump_electric_W']/r['other_auxiliary_W'],
  'additional_cost_scenarios':[additional_cost_case(accounting_reference,z,*f) for f in [(0.,0.),(.5,.5),(1.,1.)]],
  'temperature_loss_screen_K':rad/(warm['mass_flow_kg_s']*warm['gas_properties']['cp_effective_J_kg_K']),
  'scope':'Fixed nominal stream temperatures; radiation/friction losses treated as perturbations, no recovered pump-heat credit. No turbine, junction, sector-view-factor, electric transmission, or pipe deployment design.'}

def main():
 r=read_reference();rows=[layout(r,4**n,g) for g in GASES for n in range(2,9)]
 best={g:min((d for d in rows if d['gas']==g),key=lambda d:d['warm_spokes']['target_penalty_W']+d['local_cold_loops']['target_penalty_W']) for g in GASES}
 better_insulation=[layout(r,4**n,g,emissivity=.0001,insulation_areal_mass=10.) for g in GASES for n in range(2,9)]
 better_best={g:min((d for d in better_insulation if d['gas']==g),key=lambda d:d['warm_spokes']['target_penalty_W']+d['local_cold_loops']['target_penalty_W']) for g in GASES}
 sensitivity=[]
 for g,d in best.items():
  N=d['sector_count']
  for eps in [.0001,.01]:sensitivity.append(layout(r,N,g,emissivity=eps))
  for p in [1e4,1e6]:sensitivity.append(layout(r,N,g,pressure=p))
  for sig in [4e9,8e9]:sensitivity.append(layout(r,N,g,stress=sig))
  for friction in [.005,.02]:sensitivity.append(layout(r,better_best[g]['sector_count'],g,emissivity=.0001,insulation_areal_mass=10.,darcy_f=friction))
 same_quality=[];same_best={}
 for phi in [.8,.85,.9]:
  for gas in GASES:
   candidates=[]
   for tout in [90.,120.,150.,180.,210.,230.,250.,270.]:
    props=stream_properties(gas,290.,tout,1e5)
    tc=(1-.7/phi)*.98*props['enthalpy_span_J_kg']/props['entropy_span_J_kg_K']
    if tc<=30:continue
    for sectors in [64,256,1024,4096,16384]:
     z=layout(r,sectors,gas,emissivity=.0001,insulation_areal_mass=10.,hot_outlet=tout,engine_cold=tc,relative_Carnot=phi)
     candidates.append(z)
     same_quality.append({k:z[k] for k in ['relative_Carnot','gas','hot_outlet_K','engine_cold_K','sector_count','total_transport_mass_kg','total_pump_electric_W','lost_reference_generation_from_warm_radiation_W']})
   same_best[str(phi)+'_'+gas]=min(candidates,key=lambda d:d['warm_spokes']['target_penalty_W']+d['local_cold_loops']['target_penalty_W'])
 spin=[]
 for fluid in ['ParaHydrogen','OrthoHydrogen']:
  z=same_best['0.85_H2'];spin.append(layout(r,z['sector_count'],fluid,emissivity=.0001,insulation_areal_mass=10.,hot_outlet=z['hot_outlet_K'],engine_cold=z['engine_cold_K'],relative_Carnot=.85))
 out={'scope':__doc__,'CoolProp_version':CoolProp.__version__,'spin_isomer_sensitivity':spin,'gas_scenarios':GASES,'fluid_properties_source':'https://doi.org/10.1063/1.3160306','stream_assumptions':{'hot_inlet_K':290.,'hot_outlet_K':150.,'warm_source_K':300.,'engine_cold_K':35.,'radiator_K':30.,'hot_working_temperature_over_stream_temperature':.98,'fraction_of_local_Carnot':.85,'credited_heat_conversion':.7},
  'equal_quality_family_summary':same_quality,'equal_quality_best':same_best,'cases':rows,'best_in_discrete_family':best,'better_insulation_cases':better_insulation,'best_with_better_insulation':better_best,'sensitivities':sensitivity,
  'input_sha256':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in ['network-plant.json','network-ion-budget.json']}}
 (HERE/'warm-distribution.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
 print('eta stream',finite_stream_efficiency())
 for key,z in same_best.items():print('same quality',key,'N',z['sector_count'],'Tout',z['hot_outlet_K'],'Tc',z['engine_cold_K'],'mass',z['total_transport_mass_kg'],'pump plus radgen',z['total_pump_electric_W']+z['lost_reference_generation_from_warm_radiation_W'])
 for g,z in {**{g+' baseline':z for g,z in best.items()},**{g+' better MLI':z for g,z in better_best.items()}}.items():
  print(g,'N',z['sector_count'],'sector km',z['sector_equivalent_radius_m']/1000,'mass',z['total_transport_mass_kg'],'pump',z['total_pump_electric_W'],'radiation',z['warm_spokes']['radiative_heat_loss_W'],'radmassfraction',z['transport_mass_over_allocated_radiator_mass'],'auxfraction',z['pump_and_radiative_load_over_auxiliary_allowance'],'postpumpeta',z['source_to_generation_after_pumps_no_pump_heat_credit'],'mWkgfulladditional',z['additional_cost_scenarios'][0]['net_W_kg']*1000)
if __name__=='__main__':main()
