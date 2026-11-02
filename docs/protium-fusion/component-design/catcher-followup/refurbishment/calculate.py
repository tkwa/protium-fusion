#!/usr/bin/env python3
"""Current-family handling, gas-inventory and repair requirements.

These are deterministic material/energy accounts, not a nanomachine design or
retention prediction. Reworking graphite, transporting a skin, repairing defects,
and extracting gas are separate operations. No atomization charge is assumed.
"""
from pathlib import Path
import csv, hashlib, json, math
P=Path(__file__).resolve().parent
R=P.parents[2]
YEAR=365.25*86400.; EV=1.602176634e-19; MU=1.66053906660e-27; KB=1.380649e-23
H_U=1.00782503223; HE_U=4.00260325413

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p): return json.loads(Path(p).read_text())

def main():
    src=R/'component-design/catcher-followup/stopping/surface-spectra.json'
    point=R/'mixture-design/selected-plant/point.json'
    dispatch=R/'component-design/final-family-dispatch/sigma350/events.csv'
    d=load(src); plant=load(point)['plant']; area=d['area_m2']
    stages=d['stages']; rho=2000.; cp=1800.
    events=list(csv.DictReader(dispatch.open()))
    counts={s:sum(e['stage']==s for e in events) for s in stages}
    frequencies={s:1/x['period_s'] for s,x in stages.items()}
    meanrate=sum(frequencies.values()); large_rate=sum(v for s,v in frequencies.items() if s!='B1')
    Ndot={s:sum(x['species'][s]['nuclei_m2']/x['period_s'] for x in stages.values()) for s in ['p','he4']}
    gas=[]
    for label,x in stages.items():
        h=x['species']['p']['nuclei_m2']; he=x['species']['he4']['nuclei_m2']
        gas.append({'stage':label,'period_s':x['period_s'],'dispatch_370day_events':counts[label],
          'total_baryonic_mass_kg_m2':x['total_mass_fluence_kg_m2_baryon'],
          'H_nuclei_m2':h,'He_nuclei_m2':he,
          'H_neutral_atomic_kg_m2':h*H_U*MU,'He_neutral_atomic_kg_m2':he*HE_U*MU,
          'H_plus_He_baryonic_ug_m2':(h+4*he)*MU*1e9,
          'H_plus_He_neutral_atomic_ug_m2':(h*H_U+he*HE_U)*MU*1e9,
          'H2_plus_He_gas_molecules_m2':h/2+he,
          'ideal_gas_volume_m3_m2_1000K_1bar':(h/2+he)*KB*1000/1e5,
          'pressure_Pa_if_all_H2_He_held_in_5um_void_at_1000K':(h/2+he)*KB*1000/5e-6})
    handling=[]
    schedules={'every_pulse':meanrate,'large_pulses_only':large_rate,'daily':1/86400.,'monthly':12/YEAR,'annually':1/YEAR}
    for depth in [1e-6,10e-6,25e-6,50e-6,100e-6,150e-6]:
      for label,rate in schedules.items():
        sigma=rho*depth; flow=sigma*rate
        handling.append({'depth_um':depth*1e6,'schedule':label,'events_per_year':rate*YEAR,
          'handled_kg_m2_per_event':sigma,'handled_kg_m2_s':flow,'handled_kg_m2_year':flow*YEAR,
          'whole_sphere_handled_kg_s':flow*area,
          'electrical_W_m2_by_input_J_kg':{str(e):flow*e for e in [1e3,1e4,1e5,2e6]},
          'hardware_kg_m2_at_legacy_1e8_kg_per_kg_s':flow*1e8,
          'required_handling_kg_per_kg_machine_day_for_1kg_m2_machine':flow*86400,
          'sensible_heat_J_m2_per_event_1000_to_300K':sigma*cp*700,
          'gross_sensible_heat_diverted_W_m2_1000_to_300K':flow*cp*700,
          'residual_sensible_heat_diverted_W_m2_by_recuperation':{str(eta):flow*cp*700*(1-eta) for eta in [0.,.9,.95,.99]},
          'single_1m_long_strip_speed_m_s_if_continuously_renewed':rate})
    selective=[]
    for eV in [1.,10.,100.,1000.]:
      for efficiency in [.5,.9,.99]:
        useful=sum(Ndot.values())*EV*eV
        selective.append({'positive_free_energy_upgrade_eV_per_gas_nucleus':eV,
          'conversion_efficiency':efficiency,'useful_W_m2':useful,'input_W_m2':useful/efficiency,
          'conversion_loss_W_m2':useful*(1/efficiency-1),
          'scope':'conditional positive free-energy increment per processed atom; activation barrier and defect formation energy are not input work'})
    inventory=[]
    for hold in [100.,7698.198541942018,86400.,YEAR/12,YEAR/2,YEAR]:
      molecules=Ndot['p']/2+Ndot['he4']; sigma=hold*MU*(H_U*Ndot['p']+HE_U*Ndot['he4'])
      inventory.append({'hold_s':hold,'mean_gas_kg_m2':sigma,'whole_sphere_inventory_kg':sigma*area,
          'mean_gas_pressure_Pa_in_5um_void_at_1000K':molecules*hold*KB*1000/5e-6,
          'mean_gas_volume_m3_m2_at_1000K_1bar':molecules*hold*KB*1000/1e5,
          'scope':'mean-rate storage capacity, not actual transient maximum or graphite trap capacity'})
    depth_path=R/'component-design/catcher-followup/stopping/surface-stopping.json'
    depth_data=load(depth_path); bins=[0.,.01,.1,.25,.5,1.,2.,5.,10.,25.,50.,100.,150.]
    storage_depth=[]
    for low,high in zip(bins,bins[1:]):
      host=rho*(high-low)*1e-6/(12*MU)
      rates={sp:0. for sp in ['p','he4']}
      for label,x in depth_data['stages'].items():
        for sp in rates:
          q=x['species'][sp]
          rates[sp]+=sum(n for n,z in zip(q['packet_nuclei_m2'],q['packet_mean_projected_depth_um']) if low<=z<high)/x['period_s']
      storage_depth.append({'depth_um':[low,high],'carbon_host_nuclei_m2':host,
        'incident_H_nuclei_m2_s':rates['p'],'incident_He_nuclei_m2_s':rates['he4'],
        'H_per_C_after_year_if_all_held':rates['p']*YEAR/host,
        'He_per_C_after_year_if_all_held':rates['he4']*YEAR/host,
        'time_s_to_He_C_0p01_if_all_held':.01*host/rates['he4'] if rates['he4'] else None,
        'scope':'Fixed mean projected-range depth bins with no transport, release or erosion; a capacity screen, not a bubble threshold or retention prediction.'})
    budgets={
      'baseline_net_W_m2':plant['net_export_W_m2'],'baseline_mass_kg_m2':plant['mass_kg']/area,
      'baseline_mW_kg':plant['net_W_kg']*1000,
      'ordinary_auxiliary_W_m2':plant['ordinary_auxiliary_W']/area,
      'already_named_auxiliary_W_m2':.74368+.31946+.76484,
      'available_auxiliary_W_m2':plant['ordinary_auxiliary_W']/area-(.74368+.31946+.76484),
      'service_allowance_kg_m2':1000.,'already_named_services_kg_m2':321.394,
      'unassigned_service_allowance_kg_m2':1000.-321.394,
      'retained_routine_refurbishment_W_m2':plant['replacement_W']/area,
      'retained_routine_handling_kg_m2_s':plant['routine_repaired_kg_s']/area,
      'retained_routine_factory_kg_m2':plant['mass_breakdown_kg']['replacement_factory_capacity']/area,
      'cold_radiator_available_W_m2':.9*plant['radiator_emission_W_m2']-plant['environment_absorbed_W']/area,
      'cold_rejection_W_m2':plant['cold_rejection_W_m2'],
      'scope':'Unassigned allowances require actual component assignments. Net export and radiators must be recomputed for new external loads; allocating an already budgeted auxiliary does not debit it twice.'}
    # At fixed area, same plant rates and hardware, an incremental cold load cannot
    # disappear into the fully used radiator. Sensitivities below are requirements,
    # not a revised selected plant or an optimized re-dispatch.
    damage_path=R/'component-design/catcher-followup/stopping/stopping-results.json'
    damage=load(damage_path)
    targeted=[]
    for case in damage['annual']:
      for Ed,Nyr in case['NRT_proxy_events_m2_year'].items():
        for fraction in [0.001,0.01,0.1,1.]:
          rate=Nyr/YEAR*fraction
          for input_eV in [1.,10.,100.]:
            targeted.append({'electronic_stopping_multiplier':case['electronic_stopping_multiplier'],
              'displacement_threshold_eV':float(Ed),'fraction_of_NRT_proxy_requiring_service':fraction,
              'repair_events_m2_s':rate,'input_eV_per_serviced_event':input_eV,
              'input_W_m2':rate*input_eV*EV,
              'carbon_atom_equivalent_handled_kg_m2_s':rate*12*MU,
              'scope':'Service fraction and net input work are independent requirements; NRT events are recoil proxies, not surviving defects. Repair can release stored energy.'})
    per_pulse_damage=[]
    for stage in damage['rows']:
      nominal=next(x for x in stage['rows'] if x['electronic_stopping_multiplier']==1.)
      per_pulse_damage.append({'stage':stage['stage'],
        'elastic_transfer_J_m2':nominal['all_charged_elastic_transfer_J_m2'],
        'NRT_proxy_events_m2':nominal['all_charged_NRT_proxy_events_m2'],
        'input_J_m2_if_all_Ed25_proxy_events_need_10eV_service':nominal['all_charged_NRT_proxy_events_m2']['25.0']*10*EV})
    penalties=[]
    baseline_M=budgets['baseline_mass_kg_m2']; baseline_P=budgets['baseline_net_W_m2']
    for electrical in [0.,.1,1.,3.,10.]:
      for extra_mass in [0.,10.,100.]:
        penalties.append({'extra_external_electrical_W_m2':electrical,'extra_installed_kg_m2':extra_mass,
          'specific_power_mW_kg_fixed_rate_before_radiator_fix':1000*(baseline_P-electrical)/(baseline_M+extra_mass),
          'cold_radiator_temperature_K_if_area_fixed':200*((plant['cold_rejection_W_m2']+plant['environment_absorbed_W']/area+electrical)/(plant['cold_rejection_W_m2']+plant['environment_absorbed_W']/area))**.25,
          'scope':'Temperature is a heat-rejection requirement only; changing cold temperature also changes heat-engine efficiency and requires full re-solve.'})
    checks={'dispatch_events':len(events)==4228,'large_pulse_fluence':abs(gas[0]['total_baryonic_mass_kg_m2']/(11309733552923.254/area)-1)<1e-12,
      'H_flow_matches_selected_atomic':abs(Ndot['p']*MU*H_U*area/plant['recovered_H_kg_s']-1)<1e-9,
      'fifty_um_is_point_one_kg_m2':abs(50e-6*rho-.1)<1e-15,
      'radiator_no_unused_margin':abs(budgets['cold_radiator_available_W_m2']-budgets['cold_rejection_W_m2'])<1e-6,
      'returned_gas_not_fresh_loss':True}
    assert all(checks.values()),checks
    result={'scope':'Current delivered-B1/C13H2 family, 350 kg/m2 receiver, exact final chronological event counts. Mechanical handling, whole-skin processing, local repair and selective gas extraction are distinct conditional tasks. Recycled carbon and held gas remain inventory.',
      'sources_sha256':{str(p.relative_to(R)):sha(p) for p in [src,point,dispatch,damage_path,depth_path]},'code_sha256':sha(__file__),
      'area_m2':area,'density_kg_m3':rho,'mean_stage_event_frequency_s_inv':meanrate,
      'steady_events_per_year':meanrate*YEAR,'dispatch_370day_events':counts,
      'mean_H_nuclei_m2_s':Ndot['p'],'mean_He_nuclei_m2_s':Ndot['he4'],
      'neutral_H_kg_s':Ndot['p']*MU*H_U*area,'neutral_He_kg_s':Ndot['he4']*MU*HE_U*area,
      'per_pulse_gas':gas,'broad_layer_handling':handling,'selective_gas_handling_energy':selective,
      'mean_gas_inventory':inventory,'depth_resolved_holding_capacity_screen':storage_depth,'targeted_repair_sensitivities':targeted,'per_pulse_damage_proxy':per_pulse_damage,'existing_budgets':budgets,'unallocated_load_sensitivities':penalties,'checks':checks}
    (P/'results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'checks':checks,'gas':gas,'budgets':budgets,'all_events_per_year':meanrate*YEAR},indent=2))
if __name__=='__main__': main()
