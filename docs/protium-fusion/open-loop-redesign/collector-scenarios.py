"""Refresh receiver cycle weighting for explicit fuel-carbon reserve cases.

Preserves the original collector/process baselines. No hydro, EOS or thermal
convolution rerun is required when the underlying pulse records are identical.
"""
from pathlib import Path
import copy,hashlib,json,math
HERE=Path(__file__).resolve().parent
YEAR=365.25*86400
AREA=4*math.pi*1e18

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def receiver_for_cycle(base,cycle):
    r=copy.deepcopy(base)
    for k,s in r['stages'].items():
        if s['source_sha256']['record']!=cycle['stage_pulse_data'][k]['sha256']:
            raise ValueError('Receiver pulse must be refreshed for '+k)
    r['occupied_s_per_He_at_R1e9']={mode:sum(s['timing_at_R1e9'][mode]['occupied_s']*cycle['stage_pulse_data'][k]['pulses_per_He'] for k,s in r['stages'].items()) for mode in ['mass999','energy99','full_resolved']}
    r['He_rate_at_80pct_occupancy_per_s']={mode:.8/x for mode,x in r['occupied_s_per_He_at_R1e9'].items()}
    rate=r['He_rate_at_80pct_occupancy_per_s']['mass999']
    reserve=cycle.get('net_structural_carbon_allowance_kg_per_He',0.)
    wall=reserve*rate*YEAR/AREA
    r.pop('wall_carbon_budget_at_R1e9_mass999_cap',None)
    r['selected_wall_carbon_budget_at_R1e9_mass999_cap']={
        'He_rate_s':rate,'fresh_H_kg_s':rate*cycle['fresh_hydrogen_kg_per_He'],
        'fuel_C_H_mass_ratio':cycle['fresh_C_H_mass_ratio'],
        'total_ongoing_C_H_mass_budget':cycle.get('total_ongoing_C_H_mass_budget'),
        'reserved_net_structural_C_H_mass_ratio':cycle.get('net_structural_C_H_allowance',0.),
        'maximum_net_wall_carbon_loss_kg_m2_year':wall,
        'maximum_fraction_lost_from_retired_Sigma5000_at_1e_4_per_year':wall/.5,
        'scope':'Entire structural reserve assigned here to wall loss only to obtain an upper allowance; coatings and other net carbon losses share the same reserve.'}
    r['annual_surface_fuel_loading_at_R1e9_mass999_cap_kg_m2']={k:rate*cycle['stage_pulse_data'][k]['pulses_per_He']*YEAR*cycle['stage_pulse_data'][k]['mass_kg']/AREA for k in r['stages']}
    r['scenario_input_sha256']={'collector.json':sha(HERE/'collector.json'),'process-scenarios.json':sha(HERE/'process-scenarios.json')}
    r['scenario_scope']='Pulse-level receiver calculation unchanged; per-helium timing and material flows refreshed for the selected process carbon reserve.'
    return r

def main():
    process=json.loads((HERE/'process-scenarios.json').read_text())
    base=json.loads((HERE/'collector.json').read_text())
    depth=json.loads((HERE/'collector-depth.json').read_text())
    photon_map={k:next(x['total_attenuation_depth']['prompt_peak_rise_K'] for x in v['monoenergetic_sensitivities'] if x['photon_energy_keV']==10) for k,v in depth['stages'].items()}
    cases={target:receiver_for_cycle(base,cycle) for target,cycle in process['families']['b2_64']['cases'].items()}
    out={'scope':__doc__,'families':{'b2_64':{'cases':cases}},
         'photon_depth_10keV_total_mu_K':photon_map,
         'photon_scope':'Sensitivity only: unknown escaping spectrum and prescribed exponential heat deposition; baseline remains all-surface.',
         'input_sha256':{p.name:sha(p) for p in [HERE/'process-scenarios.json',HERE/'collector.json',HERE/'collector-depth.json']},
         'script_sha256':sha(Path(__file__))}
    (HERE/'collector-scenarios.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    selected={'cycle':process['families']['b2_64']['cases']['0.0056'],'receiver':cases['0.0056'],
              'selected_family':'b2_64','selected_fuel_C_H_mass_ratio':.0056,
              'photon_depth_10keV_total_mu_K':photon_map,
              'source_sha256':{'process-scenarios.json':sha(HERE/'process-scenarios.json'),'collector-scenarios.json':sha(HERE/'collector-scenarios.json')}}
    (HERE/'station-inputs.json').write_text(json.dumps(selected,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'selected_occupied_s_per_He':cases['0.0056']['occupied_s_per_He_at_R1e9'],
        'rate_caps':cases['0.0056']['He_rate_at_80pct_occupancy_per_s'],
        'wall':cases['0.0056']['selected_wall_carbon_budget_at_R1e9_mass999_cap'],
        'photon_map':photon_map},indent=2))
if __name__=='__main__':main()
