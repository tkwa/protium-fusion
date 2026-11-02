#!/usr/bin/env python3
"""Conditional catcher service interfaces through the fixed-family plant ledger.

Actual pulse rate and dispatch remain unchanged. Radiator temperature, conversion,
support, spare, factory and heat-exchanger accounts are re-solved. Component mass,
idle power and heat-leak values are explicit acceptance targets, not a completed
hot-pump design. The archived selected plant and its solver remain unchanged.
"""
from pathlib import Path
import copy,hashlib,json,math,sys,types
from scipy.optimize import brentq
P=Path(__file__).resolve().parent; R=P.parents[2]
sys.path.insert(0,str(R/'mixture-design'))
import plant_compare as pc
import station_generalized as original
import station_catcher as extended
sys.path.insert(0,str(P.parent/'gas/cold-interface-followup'))
from bank_cycle import bank_cycle
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text())
point_path=R/'mixture-design/selected-plant/point.json';point=read(point_path);baseplant=point['plant']
family_path=R/'component-design/final-family/delivered-B1-C13H2/integration.json';family=read(family_path)
refurb_path=P.parent/'refurbishment/results.json';refurb=read(refurb_path)
gas_path=P.parent/'gas/plant-handoff.json';gas=read(gas_path)
geometry_path=P.parent/'geometry/results.json'; geometry=read(geometry_path)
mechanical_path=P.parent/'warm-mechanical-cleanup/results.json'
buffer_path=P.parent/'warm-mechanical-cleanup/buffer.json';buffer=read(buffer_path)
windage_path=P.parent/'warm-mechanical-cleanup/windage/results.json'
bank_path=P.parent/'gas/cold-interface-followup/bank_cycle.py'
cold_interface_path=P.parent/'warm-mechanical-cleanup/cold-interface-correction/results.json'
A=baseplant['area_m2'];selected=family['plant_selected'];cycle=point['cycle'];rate=baseplant['helium_per_s']
freeM=refurb['existing_budgets']['unassigned_service_allowance_kg_m2']
joint=gas['local_pump']['joint_port']; active=joint['required_active_panel_area_fraction_Hgetter_Hecryo'][1]
scenario=None

def loads(s,Tr):
    if s['architecture']=='direct_cold_local_port':
        port=joint['shared_mouth_area_fraction']*original.SB*(1000.**4-Tr**4)
    else:
        # Resolved hot-buffer molecular conduit: all admitted hot radiation
        # and solid conduction are charged to the cold account.
        port=s['warm_compressed_outlet_hot_to_cold_W_m2']
    # At most the stipulated fraction of warm-port light reaches the 5K bank.
    gas_enthalpy=1.380649e-23*(3.5*refurb['mean_H_nuclei_m2_s']/2+2.5*refurb['mean_He_nuclei_m2_s'])*(1000.-Tr)
    qcold=active*original.SB*(Tr**4-5.**4)+port*1e-4
    eta=s['refrigerator_fraction_Carnot']
    radiation_work=qcold*(Tr-5.)/(eta*5.)
    CA=gas['local_pump']['regeneration']['cycling_heat_capacity_J_per_active_m2_K']
    period=gas['local_pump']['regeneration']['active_bank_capture_period_s']
    regeneration_work=active*CA*(Tr*math.log(85./5.)-(85.-5.))/(period*eta)
    banks=bank_cycle(sink_K=Tr,fraction_Carnot=eta,active_area=active,heat_capacity_J_per_active_m2_K=CA,period_s=period,electrical_warming=True)
    regeneration_work+=banks['additional_to_old_85K_sensible_cycle_W_m2']
    # Includes getter sensible regeneration and the original compression duty.
    # Keeping the latter despite warm precompression is conservative and tiny.
    fixedgas=sum(gas['local_pump']['named_electrical_W_m2'][k] for k in [
      'H_getter_sensible_regeneration_if_electrically_heated',
      'gas_300K_isothermal_compression_to1bar_50percent'])
    gas_cooling=gas['local_pump']['named_electrical_W_m2']['He_gas200to5K_refrigeration_50percent']*.5/eta*Tr/200.
    inherited=.74368+.31946*(Tr/200.)**4*(Tr-65.)/135.+.76484*(Tr-65.)/135.
    warm=s['ambient_cleanup_electrical_W_m2']+s['local_warm_booster_electrical_W_m2']
    repair_extra=maintenance['all_Ed25_proxy_defects_at10eV_electrical_per_event']*(s['repair_electrical_work_multiplier']-1.)
    named=sum(maintenance.values())+repair_extra+warm+radiation_work+regeneration_work+fixedgas+gas_cooling+s['downstream_extra_electrical_reserve_W_m2']
    return {'port_hot_to_cold_W_m2':port,'bearing_support_hot_to_cold_W_m2':s['warm_bearing_support_bypass_W_m2'],
      'additional_internal_hot_to_cold_W_m2':port+s['warm_bearing_support_bypass_W_m2']+gas_enthalpy,
      'H2_He_sensible_enthalpy_1000_to_radiator_upper_W_m2':gas_enthalpy,
      'He5K_cold_duty_W_m2':qcold,'He5K_radiation_refrigerator_electrical_W_m2':radiation_work,
      'He_sorber_full_inactive_bank_electrical_W_m2':regeneration_work,
      'He_bank_cycle':banks,
      'gas_getter_compression_cooling_electrical_W_m2':fixedgas+gas_cooling,
      'new_named_electrical_W_m2':named,'already_named_electrical_W_m2':inherited,
      'warm_motor_electrical_W_m2':warm,'internal_hot_return_W_m2':warm*s['warm_motor_fraction_recovered_as_hot_heat']}

def augmented(base,cy,**pars):
    if scenario is None:return extended.evaluate(base,cy,**pars)
    b=copy.deepcopy(base);s=scenario;L=loads(s,pars['radiator_K'])
    b['mass_breakdown_kg']['processing_structure_thermal_service']+=A*max(0.,s['new_service_hardware_kg_m2']-freeM)
    b['mass_breakdown_kg']['peak_driver_hardware']*=s['source_peak_hardware_multiplier']
    pars['bypass_W_m2']=baseplant['parameters']['bypass_W_m2']+L['additional_internal_hot_to_cold_W_m2']
    pars['internal_hot_return_W_m2']=L['internal_hot_return_W_m2']
    captured=pars['energy_capture_fraction']*(cy['H_nonnu_MeV_per_He']+cy['U_prepared_MeV_per_He'])*original.MEV*b['helium_per_s']
    hot=pars['hot_eligible_fraction']*captured-pars['bypass_W_m2']*A+L['internal_hot_return_W_m2']*A
    margin=pars['heat_efficiency']*hot-cy['U_prepared_MeV_per_He']*original.MEV*b['helium_per_s']/pars['driver_efficiency']
    if margin<=0:return {'feasible':False,'reasons':['nonpositive_process_margin']}
    pars['auxiliary_fraction']=max(baseplant['parameters']['auxiliary_fraction'],(L['already_named_electrical_W_m2']+L['new_named_electrical_W_m2'])*A/margin)
    if pars['auxiliary_fraction']>=1:return {'feasible':False,'reasons':['named_auxiliaries_exceed_process_margin']}
    ans=extended.evaluate(b,cy,**pars)
    if 'net_W'in ans:
      ans['catcher_named_account']={**L,
        'unassigned_auxiliary_after_named_W_m2':ans['ordinary_auxiliary_W']/A-L['already_named_electrical_W_m2']-L['new_named_electrical_W_m2'],
        'new_service_hardware_kg_m2':s['new_service_hardware_kg_m2'],
        'explicit_rework_factory_requirement_kg_m2_at_legacy_coefficient':rework_mass,
        'unassigned_service_after_allocation_kg_m2':max(0.,freeM-s['new_service_hardware_kg_m2']),
        'additional_service_mass_outside_allowance_kg_m2':max(0.,s['new_service_hardware_kg_m2']-freeM)}
    return ans
pc.cold=types.SimpleNamespace(evaluate=augmented,MEV=original.MEV)

def at(T):
    pars=copy.deepcopy(baseplant['parameters']);pars['radiator_K']=T
    return pc.station(selected,cycle,rate,baseplant['radius_m'],cold_parameters=pars,absorber_sigma=350.)

baseline=at(200.)
checks={'baseline_'+k:abs(baseline[k]/baseplant[k]-1)<1e-10 for k in ['mass_kg','net_W','helium_per_s','cold_rejection_W','radiator_area_fraction']}
assert all(checks.values()),checks
cap=next(x for x in geometry['routine_full_layer_rework_sensitivities'] if x['depth_um']==.1 and x['schedule']=='each_pulse')['electric_W_m2_at2MJkg']
repair_row=next(x for x in refurb['targeted_repair_sensitivities'] if x['electronic_stopping_multiplier']==1 and x['displacement_threshold_eV']==25 and x['fraction_of_NRT_proxy_requiring_service']==1 and x['input_eV_per_serviced_event']==10)
monthly_row=next(x for x in refurb['broad_layer_handling'] if x['depth_um']==150 and x['schedule']=='monthly')
monthly=monthly_row['electrical_W_m2_by_input_J_kg']['2000000.0']
gas_work=next(x for x in refurb['selective_gas_handling_energy'] if x['positive_free_energy_upgrade_eV_per_gas_nucleus']==100 and x['conversion_efficiency']==.5)['input_W_m2']
maintenance={'full100nm_face_each_pulse_at2MJkg':cap,'all_Ed25_proxy_defects_at10eV_electrical_per_event':repair_row['input_W_m2'],
  'full150um_monthly_hot_service_at2MJkg':monthly,'selective_gas_100eV_useful_at50pct':gas_work}
rework_mass=(cap+monthly)/2e6*1e8+repair_row['carbon_atom_equivalent_handled_kg_m2_s']*1e8
reference={'architecture':'warm_compression_buffer', 'refrigerator_fraction_Carnot':.3,
  'ambient_cleanup_electrical_W_m2':1.,'local_warm_booster_electrical_W_m2':.13,
  'warm_bearing_support_bypass_W_m2':.5,'warm_compressed_outlet_hot_to_cold_W_m2':buffer['cold_conduit']['all_admitted_hot_radiation_W_m2_upper']+buffer['cold_conduit']['solid_conduction_W_m2'],
  'downstream_extra_electrical_reserve_W_m2':.005,'warm_motor_fraction_recovered_as_hot_heat':0.,
  'new_service_hardware_kg_m2':300.,'source_peak_hardware_multiplier':1.,'repair_electrical_work_multiplier':1.}

def solve(s):
    global scenario
    scenario=s
    q=at(200.)
    if 'radiator_area_fraction' not in q:return {'scenario':s,'feasible':False,'reasons':q['reasons']}
    if q['radiator_area_fraction']>.9:
        upper=at(249.)
        if 'radiator_area_fraction'not in upper or upper['radiator_area_fraction']>.9:
            return {'scenario':s,'feasible':False,'reasons':['no_fixed_rate_finite_engine_radiator_bracket']}
        T=brentq(lambda x:at(x)['radiator_area_fraction']-.9,200.,249.,xtol=1e-10)*(1+1e-12)
        q=at(T)
    if 'closure_relative_residuals'in q:
        assert max(abs(v) for v in q['closure_relative_residuals'].values())<1e-9
        assert q['catcher_named_account']['unassigned_auxiliary_after_named_W_m2']>-1e-9
    return {'scenario':s,'feasible':q['feasible'],'net_mW_kg':q['net_W_kg']*1000,'plant':q}

rows=[]
for architecture in ['warm_compression_buffer','direct_cold_local_port']:
    center={**reference,'architecture':architecture}
    if architecture=='direct_cold_local_port':center['local_warm_booster_electrical_W_m2']=0.
    rows.append({'label':architecture+'_reference',**solve(center)})
    for key,values in [
      ('refrigerator_fraction_Carnot',[.15,.5]),
      ('local_warm_booster_electrical_W_m2',[.52]),
      ('source_peak_hardware_multiplier',[1.5,2.]),
      ('ambient_cleanup_electrical_W_m2',[.1,3.,10.,30.,45.,100.,220.]),
      ('warm_bearing_support_bypass_W_m2',[.1,1.,5.,10.]),
      ('new_service_hardware_kg_m2',[600.,1000.,2000.]),
      ('repair_electrical_work_multiplier',[.1,10.]),
      ('warm_motor_fraction_recovered_as_hot_heat',[.5,.9])]:
        if architecture=='direct_cold_local_port' and key=='local_warm_booster_electrical_W_m2':continue
        for val in values:rows.append({'label':architecture+'_'+key+'_'+str(val),**solve({**center,key:val})})
frontiers=[]
for architecture in ['warm_compression_buffer','direct_cold_local_port']:
  for hotreturn in [0.,.5,.9]:
    for target in [20.,10.,1.]:
      s={**reference,'architecture':architecture,'warm_motor_fraction_recovered_as_hot_heat':hotreturn}
      if architecture=='direct_cold_local_port':s['local_warm_booster_electrical_W_m2']=0.
      def f(w):
        x=solve({**s,'ambient_cleanup_electrical_W_m2':w})
        return x.get('net_mW_kg',-1000.)-target
      if f(0.)<0:
        frontiers.append({'architecture':architecture,'warm_return_fraction':hotreturn,'target_mW_kg':target,'possible_with_other_reference_parameters':False});continue
      w=brentq(f,0.,300.,xtol=1e-9)
      q=solve({**s,'ambient_cleanup_electrical_W_m2':w})
      frontiers.append({'architecture':architecture,'warm_return_fraction':hotreturn,'target_mW_kg':target,'possible_with_other_reference_parameters':True,
        'maximum_ambient_pump_electrical_W_m2':w,'radiator_K':q['plant']['parameters']['radiator_K'],
        'improvement_over_221W_m2_small_module_areal_benchmark':221./w})
# Independently test that hot electrical dissipation re-enters the thermal loop.
x=solve({**reference,'ambient_cleanup_electrical_W_m2':10.,'warm_motor_fraction_recovered_as_hot_heat':0.})
y=solve({**reference,'ambient_cleanup_electrical_W_m2':10.,'warm_motor_fraction_recovered_as_hot_heat':.9})
checks['recovered_motor_heat_raises_net_power']=y['plant']['net_W']>x['plant']['net_W']
checks['rework_factory_is_explicit_and_within_service_allocation']=14.<rework_mass<15.<reference['new_service_hardware_kg_m2']
checks['all_reported_energy_closures']=all(max(abs(v)for v in r['plant']['closure_relative_residuals'].values())<1e-9 for r in rows if 'plant'in r)
assert all(checks.values()),checks
out={'scope':__doc__,'reference_selected_point':{'net_mW_kg':baseplant['net_W_kg']*1000,'mass_kg':baseplant['mass_kg'],'net_W':baseplant['net_W']},
 'checks':checks,'maintenance_electrical_scenario_W_m2':maintenance,'explicit_rework_hardware_kg_m2_at_legacy_factory_coefficient':rework_mass,
 'uncertainty':{'refrigerator_fraction_Carnot':gas['cryogenic_efficiency_uncertainty'],
  'other_parameters':'Performance requirements and deterministic sensitivities. No probability or narrow credible interval assigned to the full hot pump, the cap routing or source capsule hardware.'},
 'scenario_status':'Conditional fixed-rate interface sensitivities. 300kg/m2 is a new-service allocation, not a BOM. Warm buffer heat leakage uses the actual-dispatch conduit calculation; .13W/m2 booster input, 1W/m2 ambient pumping and .5W/m2 bearing/support bypass remain component acceptance targets. No row replaces selected-plant/point.json.',
 'rows':rows,'ambient_power_frontiers':frontiers,
 'sources_sha256':{str(p):sha(p) for p in [point_path,family_path,refurb_path,gas_path,geometry_path,mechanical_path,buffer_path,windage_path,bank_path,cold_interface_path,Path(pc.__file__),Path(original.__file__),Path(extended.__file__)]},'code_sha256':sha(__file__)}
(P/'results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
print(json.dumps({'checks':checks,'maintenance_W_m2':maintenance,'rework_kg_m2':rework_mass,
 'references':[{'label':r['label'],'mWkg':r.get('net_mW_kg'),'Tr':r.get('plant',{}).get('parameters',{}).get('radiator_K'),'account':r.get('plant',{}).get('catcher_named_account')}for r in rows if r['label'].endswith('_reference')],
 'frontiers':frontiers},indent=2))
