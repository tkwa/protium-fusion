"""Check finite-temperature, transport and accounting equations independently."""
import json,gzip,math
import numpy as np
from scipy.optimize import brentq,minimize_scalar
from heat_transport import thermal_power,exchanger,loop,optimize_loop,additional_cost_case,RHE,CPHE
from capture_competition import read_reference,ROOT,HERE,YEAR
checks=[]
def check(name,ok,detail=None):
 checks.append({'name':name,'pass':bool(ok),'detail':detail})
 if not ok:print('FAIL',name,detail)
def close(name,a,b,rtol=1e-9,atol=1e-9):check(name,abs(a-b)<=atol+rtol*abs(b),{'actual':float(a),'expected':float(b),'relative_error':float(a/b-1) if b else None})
r=read_reference();H=thermal_power(r)
# In the unchanged selected-source family, the derivative of preaux electricity
# with respect to thermal conversion is exactly the credited hot heat stream.
from network_plant import plant
energy=json.loads((ROOT/'network-ion-budget.json').read_text());family=json.loads(gzip.decompress((ROOT/'network-channel-family.json.gz').read_bytes()))
eps=1e-6;lo=plant(energy,family,heat=.7-eps);hi=plant(energy,family,heat=.7+eps)
derivative=(hi['plant_bus_net_W']-lo['plant_bus_net_W'])/(2*eps*(1-r['assumptions']['other_fraction']))
close('credited heat from independent plant derivative',H['hot_input_W'],derivative,rtol=1e-7)
check('selected source family unchanged for derivative',lo['primary_pn_rate_s']==hi['primary_pn_rate_s'])
close('ideal helium density sanity 25K 0.1MPa',1e5/(RHE*25),1.9246,rtol=.002)
close('heat engine energy',H['electric_generation_W']+H['engine_cold_rejection_W'],H['hot_input_W'])
report=json.loads((HERE/'heat-transport.json').read_text())
for x in report['heat_exchangers']:
 if not x['feasible']:
  check('finite exchanger rejects no headroom',x['reason']=='No finite temperature approaches');continue
 phi=x['relative_Carnot'];ratio=1-.7/phi;Th=300.;Tr=30.;Qh=H['hot_input_W'];Qc=.3*Qh
 def f(T):return Qh/(Th-T)+Qc/(ratio*T-Tr)
 opt=minimize_scalar(f,bounds=(Tr/ratio+1e-9,Th-1e-9),method='bounded',options={'xatol':1e-11})
 close('analytic exchanger minimum '+str(phi),x['total_conductance_W_K'],opt.fun,rtol=1e-10)
 close('exchanger eta '+str(phi),phi*(1-x['engine_cold_K']/x['engine_hot_K']),.7)
 close('exchanger hot power '+str(phi),x['hot_conductance_W_K']*(300-x['engine_hot_K']),Qh)
 close('exchanger cold power '+str(phi),x['cold_conductance_W_K']*(x['engine_cold_K']-30),Qc)
 check('second law '+str(phi),x['engine_entropy_generation_W_K']>=-1e-4)
for i,z in enumerate(report['central_loops']+report['sensitivities']+report['alternative_exchanger_loops']+report['optimistic_cold_loop_envelopes']):
 p=z['pressure_Pa'];rho=p/(RHE*z['mean_K']);v=z['velocity_m_s'];L=z['path_length_m'];n=z['parallel_pipes'];f=z['darcy_f'];Qc=z['Qcold_W'];dT=z['temperature_span_K']
 def residual(m):
  D=math.sqrt(4*m/(math.pi*n*rho*v));specific=f*(L/D)*v*v/2
  return m*CPHE*dT-Qc-m*specific
 m0=Qc/(CPHE*dT);m=brentq(residual,m0,m0*1000)
 close('fluid balance numeric root '+str(i),z['mass_flow_kg_s'],m)
 close('flow geometry '+str(i),math.pi*z['diameter_per_pipe_m']**2/4*n*rho*v,m)
 close('fluid inventory '+str(i),z['fluid_mass_kg'],rho*math.pi*z['diameter_per_pipe_m']**2/4*n*L)
 # Direct annular tube mass, compared with the stated thin-wall approximation.
 rr=z['diameter_per_pipe_m']/2;tw=z['wall_thickness_m'];exact=math.pi*((rr+tw)**2-rr**2)*L*n*z['wall_density_kg_m3']
 close('thin-wall approximation '+str(i),z['pressure_wall_mass_kg'],exact,rtol=.0011)
 close('Darcy fluid power '+str(i),z['fluid_dissipation_W'],f*L/z['diameter_per_pipe_m']*rho*v*v/2*(m/rho))
 close('pump efficiency '+str(i),z['pump_electric_W']*z['pump_efficiency'],z['fluid_dissipation_W'])
 close('fluid heat closure '+str(i),z['fluid_heat_relative_residual'],0.,atol=1e-12)
 # Perturb velocity on either side; optimizer must have a local cost minimum.
 k=r['assumptions']['annual_replacement_fraction']*r['assumptions']['embodied_J_kg']/YEAR
 cost=z['objective_W_equivalent_at_target']
 for factor in [.99,1.01]:
  other=loop(Qc,dT,z['mean_K'],L,v*factor,p,f,n)
  oc=other['pump_electric_W']+(.001+k)*other['pipe_and_fluid_mass_kg']
  check('velocity local optimum '+str(i)+'/'+str(factor),oc>=cost*(1-1e-10))
zero={'pipe_and_fluid_mass_kg':0.,'pump_electric_W':0.}
p=additional_cost_case(r,zero)
close('zero additional costs mass',p['total_mass_kg'],r['total_mass_kg'])
close('zero additional costs power',p['net_W'],r['net_W'])
# Reducing path length at fixed velocity, zero friction removes pump loss while
# fluid inventory remains a direct residence-time calculation.
z=loop(1e9,10.,35.,1e4,10.,darcy_f=0.)
close('zero friction mass flow',z['mass_flow_kg_s'],1e9/(CPHE*10))
close('zero friction pump',z['pump_electric_W'],0.)
check('all saved optimum loops subsonic',all(z['Mach_number']<.1 for z in report['central_loops']+report['sensitivities']+report['alternative_exchanger_loops']+report['optimistic_cold_loop_envelopes']))
out={'checks':checks,'count':len(checks),'passed':sum(c['pass'] for c in checks),'all_pass':all(c['pass'] for c in checks)}
(HERE/'heat-transport-validation.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
print({k:v for k,v in out.items() if k!='checks'})
if not out['all_pass']:raise SystemExit(1)
