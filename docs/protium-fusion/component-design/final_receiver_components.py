"""Final delivered-B1/CH2 receiver support, guide and separated-prompt extension.
Exact final source128 family; hot capture and cold routing come from their separate plant ledger.
Existing guide geometry/attenuation and frame are retained as requirements.
"""
from pathlib import Path
import ast,json,math,hashlib
import numpy as np
P=Path(__file__).parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
# Reuse the frozen prompt/slow combiner, replacing only input/destination names.
original=P/'source128_thermal_combine.py';source=original.read_text()
for before,after in {"P/'source128-thin-thermal.json'":"P/'final-receiver-thermal.json'","P/'source128-thermal-combined.json'":"P/'final-receiver-combined.json'","P/'source128-prompt-proxy.json'":"P/'final-receiver-prompt-proxy.json'","P/'source128-prompt-radiance.json'":"P/'final-receiver-prompt-radiance.json'"}.items():
    assert source.count(before)==1
    source=source.replace(before,after)
exec(compile(source,str(original),'exec'),{'__file__':str(original),'__name__':'receiver250_combine'})
thermal=read(P/'final-receiver-thermal.json');combined=read(P/'final-receiver-combined.json')
trade_path=Path('/home/tkwa/src/protium-fusion/output/final-coupled-tradeoff-extended.json');trade=read(trade_path)
G=6.67430e-11;R=1.1e9;area=4*math.pi*R**2;n2=G*1.98847e30/(5*149597870700.)**3;n=math.sqrt(n2);E2D=1e8;nu=.2
# Pure audited affine energy functions only; no archived main program is run.
ns=dict(np=np,math=math,R=R,n2=n2,n=n,E2D=E2D,nu=nu)
tree=ast.parse((P/'receiver_prestress.py').read_text());subset=ast.Module(body=[z for z in tree.body if isinstance(z,ast.FunctionDef) and z.name in ['energy2','modes']],type_ignores=[]);exec(compile(subset,'prestress-pure-functions','exec'),ns)
gamma=read(P/'guide-gamma.json');guide=read(P/'guide-component.json');buffers=read(P/'guide-buffer-materials.json')
line=gamma['guide_line_density_m_m2'];coil=guide['coil_current_and_mass']['rows'][1]['finite_cos2theta_winding_outer_radius_m'];oldrad=gamma['original_cryostat_radius_m'];cop=gamma['COP65to200'];old_ir=guide['cryogenic_radiation'][1]['electric_W_m2_receiver'];kappa=gamma['carbon_minimum_mass_energy_absorption_m2_kg_through20MeV']
# Pure existing geometric buffer solver, with all inputs explicitly above.
nsb=dict(math=math,line=line,coil=coil,oldrad=oldrad,cop=cop,old_ir=old_ir,kappa=kappa,area=area)
tree=ast.parse((P/'guide_buffer_materials.py').read_text());exec(compile(ast.Module(body=[z for z in tree.body if isinstance(z,ast.FunctionDef) and z.name=='solve'],type_ignores=[]),'buffer-pure-function','exec'),nsb)
service_packet=read(P/'sector-allocations.json')['cases']['mixture_scan_003_top']
rows=[]
for sigma in [200.,250.,350.,500.]:
    cases=[r for r in trade['rows'] if r['receiver_sigma_kg_m2']==sigma and r['additional_missed_hot_fraction_HplusUprecursor'] in [0.,.005] and r['missed_hot_routing']=='cold_equipment']
    for row in cases:
        p=row['plant'];M=p['mass_kg'];body=M/area;a=G*M/(2*R**3);T=max(1e4,(4.4*n2-a)*body*R**2/2);ap=2*T/(body*R**2);b=np.array([(a+ap+4*n2)/2,(a+ap-2*n2)/2,(a+ap-4*n2)/2]);ss=body*b/(1e4/R)**2;assert min(ss)>0
        massstream=float(sum(ss)*area);energy=massstream*1e8/3;loss=.001*energy/(2*math.pi*R/1e4)
        support={'body_mass_kg':M,'body_kg_m2':body,'a_over_n2':a/n2,'effective_a_over_n2':(a+ap)/n2,'pretension_N_m':T,'pressure_Pa':2*T/R,'prestrain':.8*T/E2D,'stress_MPa_in8kg_m2_at2000density_load_layer':T/.004/1e6,'stream_axis_kg_m2':ss.tolist(),'stream_kg_m2':float(sum(ss)),'stream_vs_scalar_budget':massstream/p['mass_breakdown_kg']['support_streams'],'extra_stream_mass_kg_m2':max(0,massstream-p['mass_breakdown_kg']['support_streams'])/area,'stream_dissipation_W_m2':loss/area,'affine_cases':[ns['modes'](a,body,T,days) for days in [30.,100.,300.]]}
        # Conservatively allocate all prompt directly to cold guide surroundings,
        # even for EM20. Delayed whole heat receives D1.1 and the existing
        # localized-cascade attenuation screen. Extra.005 applies to whole pulse.
        pulse={key:1.1*sum(s['decay_fluence'].values())*math.exp(-kappa*sigma)+s['prompt_fluence']+.005*s['fluence'] for key,s in thermal['pulse_fluences'].items()}
        prep=read(P/'final-family/delivered-B1-C13H2/integration.json')['preparations']['n15']['preparation_escaped_radiation_J']/area
        pulse['B1']+=prep
        gr={'receiver_kg_m2':sigma,'mean_missed_hot_to_cold_W_m2':row['heat_in_cold_equipment_W']/area,'pulse_photon_upper_fluence_J_m2_after_main_wall':pulse}
        gb=[nsb['solve'](gr,mat,ang) for mat in buffers['materials'] if mat['name'] in ['linear_polyethylene_effective_density900','linear_polyethylene_half_accessible_capacity'] for ang in [math.pi/2,3.]]
        named=dict(service_packet['named_allocations_within_services_kg_m2'])
        named['local_guide_photon_shield_and_polyethylene_buffer']=gb[0]['shield_plus_buffer_kg_m2']
        aux=dict(service_packet['named_electric_loads_within_ordinary_aux_W_m2'])
        aux['additional_guide_photon_and_buffer_IR_refrigeration']=gb[0]['total_extra_electric_W_m2']
        allocations={'named_services_kg_m2':named,'subtotal_kg_m2':sum(named.values()),'unassigned_at_service_scenarios_kg_m2':{str(sv):sv-sum(named.values()) for sv in [400.,700.,1000.,1500.]},'services_fraction_of_current_plant_mass':1000*area/M,'ordinary_auxiliary_W_m2':p['ordinary_auxiliary_W']/area,'named_auxiliary_W_m2':aux,'remaining_auxiliary_W_m2':p['ordinary_auxiliary_W']/area-sum(aux.values()),'hot_loop_reference_W_m2':service_packet['hot_loop']['heat_W_m2'],'current_hot_loop_W_m2':p['engine_hot_W']/area,'cold_loop_reference_W_m2':service_packet['cold_loop']['heat_W_m2'],'current_cold_with_environment_W_m2':(p['cold_rejection_W']+p['environment_absorbed_W'])/area,'interpretation':'Retain the older higher-load helium/loop-HX allocations as conservative named subdivisions. Genericservices include these functions once. Separate source/converter/storage/radiator/insulation/stationHX/factory budgets remain outside this subtotal. The unassigned remainder is an allowance, not a measured bill of materials or a probabilistic bound.'}
        assert allocations['current_hot_loop_W_m2']<=allocations['hot_loop_reference_W_m2']
        assert allocations['current_cold_with_environment_W_m2']<=allocations['cold_loop_reference_W_m2']*(1+1e-9)
        rows.append({'additional_missed_hot_fraction_HplusUprecursor':row['additional_missed_hot_fraction_HplusUprecursor'],'service_allocations':allocations,'sigma_kg_m2':sigma,'prompt_policy':row['prompt_policy'],'mW_kg':row['mW_kg'],'plant_mass_kg':M,'plant_power_W':p['net_W'],'He_s':p['helium_per_s'],'occupancy':p['sequential_occupancy_mass999'],'B1_interval_s':p['pulse_interval_s']['n15'],'faster_thermal_cadence_ratio':thermal['reference']['cadence_He_s']/p['helium_per_s'],'support':support,'guide_input':gr,'guide_buffer_cases':gb})
out={'scope':__doc__,'code_sha256':sha(__file__),'inputs_sha256':{str(x):sha(x) for x in [trade_path,P/'final-receiver-thermal.json',original,P/'final-receiver-prompt-proxy.json',P/'final-receiver-prompt-radiance.json',P/'receiver_prestress.py',P/'guide-gamma.json',P/'guide-component.json',P/'guide_buffer_materials.py',P/'guide-buffer-materials.json',P/'sector-allocations.json']},'rows':rows,'limits':['Thermal load retains all source heat at faster archivedsource1282000/.99 cadence; pulse transport and guide histories remain explicitly distinct from plant capture.','Prestress carries over the20kg/m2 frame,E2D1e8N/m,nu.2,8kg/m2 load layer and remaining joints. These are material/geometry acceptance conditions, not proof of local guide-network control.','Guide shield/buffer inventories must be named within1000kg/m2 services. Their refrigeration must fit the retained ordinaryaux allowance or receive a separate debit; no free added cooling.','Guide attenuation assumes closed local graphite shield and localized secondary energy deposition. Coupledguide-specific transport, heat access/quench and radiation lifetime remain unqualified.','The entireprompt pulse is conservatively put into cold surroundings for guidebuffer sizing even inEM20 cases. Uniform angular interception is explicit; changed sourceplacement needs re-evaluation.','Exact deliveredB1/CH2 pulse family used; pulse guidebuffer sizing conservatively retains an extra.005 wholepulse allowance even for primary extra0 plant rows.']}
(P/'final-receiver-components.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
print(json.dumps({'summary':[{'sigma':r['sigma_kg_m2'],'prompt':r['prompt_policy'],'mWkg':r['mW_kg'],'mass':r['plant_mass_kg'],'B1_period':r['B1_interval_s'],'support':r['support'],'guide':r['guide_buffer_cases'][0]} for r in rows]},indent=2))
