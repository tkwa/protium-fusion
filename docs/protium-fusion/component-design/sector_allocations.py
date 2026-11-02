"""Named subdivisions of existing services and auxiliary budget; no new plant solve.
The mixture snapshot remains at its original receiver mass and pulse cadence.
"""
from pathlib import Path
import json, math, hashlib
P=Path(__file__).parent
paths={'thermal':P.parent/'cold-sink-redesign/cases.json','mixture':Path('/home/tkwa/src/protium-fusion/output/mixture-runs/plant-scan-003.json'),'guide':P/'guide-component.json'}
raw={k:json.loads(v.read_text()) for k,v in paths.items()}
points={'frozen_thermal_reference':raw['thermal']['cases']['working_200K_60pct_IR11_bypass1'],'mixture_scan_003_top':raw['mixture']['top_detailed'][0]['plant']}
cp=2.5*8.314462618/.004002602

def loop(q,T,u,L=1000.,pressure=2e6,stress=1e8):
    rho=pressure*.004002602/(8.314462618*T);mdot=q*L*L/(cp*20)
    area=mdot/(rho*u);D=math.sqrt(4*area/math.pi);length=2*L;t=pressure*D/(2*stress)
    wall=math.pi*((D/2+t)**2-(D/2)**2)*length*8000/(L*L)
    gas=area*length*rho/(L*L)
    dp=(.015*length/D+6)*rho*u*u/2
    return {'heat_W_m2':q,'mean_K':T,'deltaT_K':20.,'gas_speed_m_s':u,'pressure_Pa':pressure,'density_kg_m3':rho,'pipe_inner_diameter_m':D,'pipe_wall_thickness_m':t,'pressure_wall_allowable_Pa':stress,'wall_kg_m2':wall,'inventory_kg_m2':gas,'pump_W_m2':mdot/rho*dp/.8/(L*L)}

def lmtd(a,b): return (a-b)/math.log(a/b)

cases={}
for name,p in points.items():
    A=p['area_m2'];qh=p['hot_engine_input_W_m2'];qc=p['cold_rejection_W_m2']+p['environment_absorbed_W']/A
    hot=loop(qh,980.,10.);cold=loop(qc,215.,5.)
    ua=qh/lmtd(30,10)+qh/lmtd(40,20)+.4*qh/lmtd(25,5)+qc/lmtd(25,5)
    oldexternal=p['external_exchanger_minimum_UA_W_K']/A
    extraHX=max(0.,ua-oldexternal)/10
    flowmass=sum(x['wall_kg_m2']+x['inventory_kg_m2'] for x in [hot,cold])
    pump=sum(x['pump_W_m2'] for x in [hot,cold])
    cryo=raw['guide']['cryogenic_radiation'][1]['electric_W_m2_receiver']
    ordinary=p['ordinary_auxiliary_W']/A
    allocations={'hot_and_cold_helium_headers_and_inventory':flowmass,'additional_loop_interface_HX_over_retained_external_HX':extraHX,'cold_cross_braced_mesh_package_requirement':20.,'guide_coils_cryostats_trim_control_package_requirement':50.}
    service=sum(allocations.values());remaining=1000-service
    ports=[]
    for N in [100000,1000000]:
      for r in [100.,250.]:
        frac=N*math.pi*r*r/A
        ports.append({'spatial_outlet_count':N,'radius_m':r,'area_fraction':frac,'black_1000_to_200K_unbaffled_W_m2':frac*5.670374419e-8*(1000**4-200**4)})
    cases[name]={'original_mass_kg':p['mass_kg'],'original_net_W_kg':p['net_W_kg'],'radius_m':p['radius_m'],'area_m2':A,'selected_cases':p.get('selected_cases'),
       'hot_loop':hot,'cold_loop':cold,'loop_interface_UA_W_K_m2':ua,'retained_external_UA_W_K_m2':oldexternal,
       'named_allocations_within_services_kg_m2':allocations,'named_subtotal_kg_m2':service,'unallocated_remainder_of_1000_services_kg_m2':remaining,
       'lower_current_density_guide_package_200kg_m2_remaining_services':remaining-150,
       'named_electric_loads_within_ordinary_aux_W_m2':{'external_helium_loop_pumps':pump,'guide_radiative_cryocooling':cryo},
       'named_auxiliary_subtotal_W_m2':pump+cryo,'retained_ordinary_auxiliary_W_m2':ordinary,'remaining_ordinary_auxiliary_W_m2':ordinary-pump-cryo,
       'if_instead_additional_fixed_rate_net_W_kg':p['net_W_kg']-(pump+cryo)*A/p['mass_kg'],
       'additional_power_convention_warning':'Extra subtraction is a separate conservative sensitivity; it also raises radiator load. A radiator-saturated point must reduce cadence or raise radiator temperature. Do not both itemize within the retained margin and subtract a second time.',
       'ports':ports,
       'unchanged_separate_allocations':{'receiver_kg_m2':5000.,'insulation_kg_m2':100.,'radiator_kg_m2_emitting_area':20.,'dedicated_converter_and_internal_HX':'Retained from the station; not paid again from services.'}}
    assert remaining>0 and ordinary>pump+cryo
out={'scope':'Conditional named allocations and explicit loads inside previously budgeted services/ordinary auxiliaries. They do not establish a qualified heat engine, guide controller, catcher recycling process, or pressure-pipe material.',
 'input_sha256':{k:hashlib.sha256(p.read_bytes()).hexdigest() for k,p in paths.items()},'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'assumptions':{'sector_side_m':1000,'helium_cp_J_kg_K':cp,'hot_pipe_material_density_kg_m3':8000,'hot_and_cold_pressure_allowable_Pa':1e8,'Darcy_f':.015,'minor_loss_coefficient':6,'pump_efficiency':.8,'HX_rating_W_K_kg':10.,'mesh_ideal_Eh_N_m':1e8,'mesh_ideal_kg_m2_at100GPa_density2000':2.,'mesh_package_multiplier':10.,'guide_selected_engineering_current_A_mm2':30.,'guide_complete_conductor_kg_m2':raw['guide']['coil_current_and_mass']['rows'][1]['finite_winding_conductor_kg_m2_receiver'],'guide_effective_radiative_emissivity':.01,'manufacturing_boundary':'No blanket atomization cost. Existing factory throughput capacity and 2 MJ/kg routine maintenance sensitivity retained. Surface recovery process throughput is separate and not derived here.'},'cases':cases}
Path(__file__).with_name('sector-allocations.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
print(json.dumps({k:{f:v[f] for f in ['named_subtotal_kg_m2','unallocated_remainder_of_1000_services_kg_m2','named_auxiliary_subtotal_W_m2','remaining_ordinary_auxiliary_W_m2']} for k,v in cases.items()},indent=2))
