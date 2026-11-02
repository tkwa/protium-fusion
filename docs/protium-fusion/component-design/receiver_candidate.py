"""Package a separate2000kg/m² conditional receiver; do not optimize the plant.
Every number links to a saved spectrum/thermal/support calculation.
"""
from pathlib import Path
import json,hashlib
P=Path(__file__).parent
files={'thermal':P/'receiver-thickness.json','thermal_validation':P/'receiver-thickness-validation.json','support':P/'receiver-support-family.json','photons':P/'photon-containment/containment.json','spectrum':P/'photon-containment/spectrum.json','sector_allocations':P/'sector-allocations.json'}
x={k:json.loads(v.read_text()) for k,v in files.items()}
thermal=next(r for r in x['thermal']['thermal_rows'] if r['sigma_kg_m2']==2000 and r['k_W_m_K']==50 and r['cp_J_kg_K']==1800)
support=next(r for r in x['support']['rows'] if r['receiver_sigma_kg_m2']==2000)
photons=next(r for r in x['photons']['aggregate_localized_cascade_bounds'] if r['sigma_kg_m2']==2000)
photon=next(r for r in photons['Doppler_boost_scenarios'] if r['source_Doppler_factor']==1.1)
selected=x['thermal']['reference']['selected_cases']
assert selected==x['support']['reference_cases']
assert x['spectrum']['mixture']=={'A':3,'B1':2,'B2':2.5,'C13':4}
assert support['positive_three_family_stream_solution'] and not any(r['positive_growth_detected'] for r in support['global_affine_control_cases'])
assert x['thermal_validation']['passed']
assert photon['remaining_nonphoton_after_full_hot_reserve_fraction_for_99pct_capture']>0
out={'status':'Separate conditional component candidate;5000kg/m² remains the controlled mixture-comparison receiver.',
 'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'input_sha256':{k:hashlib.sha256(v.read_bytes()).hexdigest() for k,v in files.items()},
 'source':{'selected_cases':selected,'radius_m':x['thermal']['reference']['radius_m'],'held_comparison_helium_per_s':x['thermal']['reference']['helium_per_s'],'hydro_zones':64,'hydro_opacity_m2_kg':.001},
 'receiver':{'material':'Graphite body; renewable catcher skin remains a separate interface.','areal_mass_kg_m2':2000.,'density_kg_m3':2000.,'thickness_m':1.,'whole_body_graphite_kg':2000*4*3.141592653589793*x['thermal']['reference']['radius_m']**2,'controlled_thermal_property_scenario':{'k_W_m_K':50.,'cp_J_kg_K':1800.},'nominal_volume_average_minimum_K':1000.,'all_stage_bulk_temperature_swing_bound_K':thermal['all_stage_bulk_swing_upper_K_at_constant_mean_extraction'],'conditional_front_temperature_screen_K':thermal['finite_buffer_front_temperature_bound_K_if_each_bulk_mode_minimum_1000K'],'large_pulse_guard_s':x['thermal']['reference']['large_pulse_guard_s'],'mean_front_to_back_gradient_K':thermal['mean_front_back_gradient_K'],'thermal_scope':'Full-retention surface-deposition/ballistic/decay screen. Bulk storage is finite; local pickup temperature and variable-temperature engine operation remain to be matched.'},
 'capture':{'performance_target':.99,'conditional_local_cascade_loss_bound_fraction_of_HplusU':photon['plus_entire_D1p10_hot_direct_bypass_loss_fraction_thermal'],'remaining_loss_fraction_of_HplusU_for_target':photon['remaining_nonphoton_after_full_hot_reserve_fraction_for_99pct_capture'],'qualified_total_capture':False,'conditions':x['photons']['bound_assumptions'],'source_geometry_requirement':'Radioactive skin and processing hold-up must remain inside the credited shield or have separately enclosed shielding.'},
 'support':support,
 'retained_budgets':{'services_kg_m2':1000.,'insulation_kg_m2':100.,'source':'Archived source hardware unchanged.','electrical_storage':'Archived store unchanged. No free receiver-energy buffer in that store.','converter_and_radiator':'The parent station ledger may resize according to its declared ratings and capture/rate assumptions; none is optimized in this packet.','guide_package_kg_m2_within_services':50.,'cold_mesh_package_kg_m2_within_services':20.},
 'next_requirements':['Coupled photon/charged-secondary containment and source placement','Actual hot graphite property/joint/thermal-shock design and renewable catcher process','Finite-buffer heat pickup and converter performance over the resulting temperature cycle','Local tile/subframe/guide controller implementation within the retained services and auxiliary allowance'],
 'lighter_1000_branch':'Remains a transport-optimization sensitivity. Its current three-family support construction also requires prestress or another force topology before use.',
 'checks':{'same_selected_stage_family':True,'thermal_validation_passed':True,'positive_stream_family_densities':True,'tested_affine_controls_non_growing':True,'positive_conditional_capture_headroom':True}}
Path(__file__).with_name('receiver-candidate-2000.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
print(json.dumps({'status':out['status'],'thermal':out['receiver'],'capture_loss_upper':out['capture']['conditional_local_cascade_loss_bound_fraction_of_HplusU'],'remaining_loss_headroom':out['capture']['remaining_loss_fraction_of_HplusU_for_target'],'checks':out['checks']},indent=2))
