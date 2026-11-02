"""Independent quadrature, geometry, hydraulic and budget checks for warm spokes."""
import json,math
from scipy.integrate import quad
from scipy.optimize import brentq
from CoolProp.CoolProp import PropsSI
from warm_distribution import finite_stream_efficiency,warm_pair_loop
from capture_competition import read_reference,HERE,YEAR,SB
from heat_transport import thermal_power,loop
checks=[]
def check(name,ok,detail=None):
 checks.append({'name':name,'pass':bool(ok),'detail':detail})
 if not ok:print('FAIL',name,detail)
def close(name,x,y,rtol=1e-8,atol=1e-9):check(name,abs(x-y)<=atol+rtol*abs(y),{'actual':x,'expected':y})
r=read_reference();Q=thermal_power(r)['hot_input_W'];data=json.loads((HERE/'warm-distribution.json').read_text())
all_cases=data['cases']+data['better_insulation_cases']+data['sensitivities']+list(data['equal_quality_best'].values())+data['spin_isomer_sensitivity']
for i,c in enumerate(all_cases):
 w=c['warm_spokes'];props=w['gas_properties'];fluid=props['fluid'];p=c['pressure_Pa'];Thi=c['hot_inlet_K'];Tlo=c['hot_outlet_K'];Tc=c['engine_cold_K'];phi=c['relative_Carnot'];hf=c['hot_fraction'];N=c['sector_count'];D=w['diameter_per_pipe_m'];L=w['path_length_m'];m=w['mass_flow_kg_s'];A=N*math.pi*D*D/4
 # Integrate heat capacity independently from the enthalpy/entropy difference.
 cp=lambda T:PropsSI('Cpmass','T',T,'P',p,fluid)
 dh=quad(cp,Tlo,Thi,epsabs=1e-4,epsrel=1e-11)[0];ds=quad(lambda T:cp(T)/T,Tlo,Thi,epsabs=1e-6,epsrel=1e-11)[0]
 close('enthalpy quadrature '+str(i),props['enthalpy_span_J_kg'],dh)
 close('entropy quadrature '+str(i),props['entropy_span_J_kg_K'],ds)
 eta=quad(lambda T:phi*(1-Tc/(hf*T))*cp(T),Tlo,Thi)[0]/dh
 close('finite-stream engine quadrature '+str(i),c['finite_stream_efficiency_before_transport'],eta)
 hotK=Q/dh*quad(lambda T:cp(T)/((1-hf)*T),Tlo,Thi)[0]
 close('hot conductance quadrature '+str(i),c['hot_exchanger_conductance_W_K'],hotK)
 close('cold conductance power '+str(i),c['cold_exchanger_conductance_W_K']*(Tc-30),.3*Q)
 check('positive local engine efficiency '+str(i),hf*Tlo>Tc)
 rh=props['hot_density_kg_m3'];rl=props['cold_density_kg_m3'];rm=props['mean_density_kg_m3'];v=w['velocity_m_s'];f=c['darcy_f']
 close('hot mass flux '+str(i),A*rh*w['hot_velocity_m_s'],m)
 close('cold mass flux '+str(i),A*rl*w['cold_velocity_m_s'],m)
 close('two-leg inventory '+str(i),w['fluid_mass_kg'],A*L/2*(rh+rl))
 work=f*(L/2)/D*((m/(A*rh))**2+(m/(A*rl))**2)/2
 close('two-leg friction '+str(i),w['fluid_dissipation_W'],m*work)
 close('transport heat balance '+str(i),m*dh,Q+m*work)
 def residual(mm):
  diam=math.sqrt(4*mm/(N*math.pi*rm*v));wh=f*(L/2)/diam*((rm*v/rh)**2+(rm*v/rl)**2)/2
  return mm*dh-Q-mm*wh
 close('independent mass-flow root '+str(i),m,brentq(residual,Q/dh,Q/dh*100))
 area_ext=N*math.pi*(D+2*w['wall_thickness_m'])*L
 radiate=area_ext*c['effective_warm_emissivity']*SB*((Thi**4+Tlo**4)/2-3**4)
 close('two-leg radiation '+str(i),w['radiative_heat_loss_W'],radiate)
 exactwall=math.pi*((D/2+w['wall_thickness_m'])**2-(D/2)**2)*N*L*1600
 close('pressure wall thin limit '+str(i),w['pressure_wall_mass_kg'],exactwall,rtol=.0011)
 close('sector area partition '+str(i),N*math.pi*c['sector_equivalent_radius_m']**2,r['radiator_area_m2'])
 cold=c['local_cold_loops'];close('total mass '+str(i),c['total_transport_mass_kg'],w['pipe_and_fluid_mass_kg']+w['insulation_mass_kg']+cold['pipe_and_fluid_mass_kg'])
 close('pumps sum '+str(i),c['total_pump_electric_W'],w['pump_electric_W']+cold['pump_electric_W'])
 # Reconstruct mass ledger including that warm radiation already escaped.
 for j,z in enumerate(c['additional_cost_scenarios']):
  a=r['assumptions'];k=a['annual_replacement_fraction']*a['embodied_J_kg']/YEAR;sp=a['annual_replacement_fraction']*a['replacement_inventory_years'];M=z['total_mass_kg']
  fixed=sum(x for n,x in r['mass_breakdown_kg'].items() if n not in ['radiators','replacement_manufacturing_equipment','replacement_inventory'])+z['extra_mass_kg']
  bus=r['plant_bus_net_W']-z['extra_electric_load_W'];net=bus-k*M
  radiator=(r['nuclear_after_capture_W']-radiate-net)*a['radiator_areal_mass']/(a['radiator_emissivity']*SB*(a['radiator_K']**4-a['ambient_K']**4))
  close('mass feedback ledger '+str(i)+'/'+str(j),M,fixed+radiator+a['manufacturing_kg_W']*k*M+sp*M)
  close('electric ledger '+str(i)+'/'+str(j),z['net_W'],net)
 check('subsonic '+str(i),max(w['Mach_number'],cold['Mach_number'])<.1)
 check('low pressure drop '+str(i),max(w['pressure_drop_over_absolute_pressure'],cold['pressure_drop_over_absolute_pressure'])<.1)
for k,c in data['equal_quality_best'].items():
 close('matched conversion '+k,c['finite_stream_efficiency_before_transport'],.7)
 check('nominal friction exceeds smooth pipe estimate '+k,c['warm_spokes']['prescribed_f_exceeds_smooth_estimate'])
 fs=c['warm_spokes']['smooth_Colebrook_Darcy_f'];Re=c['warm_spokes']['mean_Reynolds_number']
 close('smooth Colebrook residual '+k,1/math.sqrt(fs)+2*math.log10(2.51/(Re*math.sqrt(fs))),0.,atol=1e-8)
# Local velocity perturbation at a representative best point catches optimizer failures.
for key,c in data['equal_quality_best'].items():
 w=c['warm_spokes'];k=r['assumptions']['annual_replacement_fraction']*r['assumptions']['embodied_J_kg']/YEAR
 for factor in [.99,1.01]:
  z=warm_pair_loop(Q,c['hot_inlet_K'],c['hot_outlet_K'],w['path_length_m'],w['velocity_m_s']*factor,c['gas'],c['pressure_Pa'],c['darcy_f'],c['sector_count'],c['pipe_working_stress_Pa'])
  surface=math.pi*(z['diameter_per_pipe_m']+2*z['wall_thickness_m'])*c['sector_count']*w['path_length_m']
  rad=surface*c['effective_warm_emissivity']*SB*((c['hot_inlet_K']**4+c['hot_outlet_K']**4)/2-3**4)
  cost=z['pump_electric_W']+.7*rad+(.001+k)*(z['pipe_and_fluid_mass_kg']+surface*c['warm_insulation_kg_m2'])
  check('local velocity minimum '+key+'/'+str(factor),cost>=w['target_penalty_W']*(1-1e-10))
out={'count':len(checks),'passed':sum(z['pass'] for z in checks),'all_pass':all(z['pass'] for z in checks),'checks':checks}
(HERE/'warm-distribution-validation.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
print({k:v for k,v in out.items() if k!='checks'})
if not out['all_pass']:raise SystemExit(1)
