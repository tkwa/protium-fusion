"""Separate receiver and station-input variant for a selected refined B2 family.

Reads completed process scenarios and raw B2 trace; never reruns hydro and
never overwrites the corrected64 collector or station-input baselines.
"""
from pathlib import Path
import argparse,copy,hashlib,importlib.util,json,math
import numpy as np
from scipy.optimize import brentq
import collector as c
HERE=Path(__file__).resolve().parent

def module(name,file):
    s=importlib.util.spec_from_file_location(name,HERE/file);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--family',required=True);ap.add_argument('--target',default='0.0056');args=ap.parse_args()
    scenarios=json.loads((HERE/'process-scenarios.json').read_text());family=scenarios['families'][args.family];cycle=family['cases'][args.target]
    tag=args.family.replace('_','-')
    path=Path(cycle['stage_pulse_data']['c12']['source']);d=json.loads(path.read_text())
    initial=Path(d['initial_profile_path']);trace=Path('/var/tmp/tkwa-chatgpt/protium-b2-completion')/path.with_suffix('.npz').name
    z=np.load(trace);t=z['time_s'];nz=int(z['zones']);nv=int(z['variables_per_zone']);state=z['state'].reshape(nz,nv,len(t))
    index=d['state_variables'].index('radiated_energy/uunit');w=np.asarray(d['profile']['mass_fraction_of_total'])
    cumulative=w@state[:,index,:];cumulative*=d['radiation_escaped_J']/cumulative[-1]
    assert np.min(np.diff(cumulative))>=-1e-7*cumulative[-1]
    item={'record':d,'initial':json.loads(initial.read_text()),'prompt_time_s':t.tolist(),'prompt_cumulative_J':cumulative.tolist(),
        'source_sha256':{'record':c.sha(path),'initial':c.sha(initial),'raw_trace':c.sha(trace)},
        'source_paths':{'record':str(path),'initial':str(initial),'raw_trace':str(trace)}}
    inputs=json.loads((HERE/'collector-inputs.json').read_text());inputs['stages']['c12']=item
    inputs.pop('B2_eight_sound_time_sensitivity',None)
    input_path=HERE/f'collector-inputs-{tag}.json';input_path.write_text(json.dumps(inputs,indent=2,allow_nan=False)+'\n')
    receiver=json.loads((HERE/'collector.json').read_text());receiver['stages']['c12']=c.stage_screen(item)
    receiver.pop('B2_eight_sound_time_sensitivity',None)
    receiver['resolution_status']=f'Completed B2{d["zones"]} with matched initial energy; A/B1/C13 remain64. Mesh convergence and pulse preparation remain separate questions.'
    receiver['scope']='Explicit refined B2 receiver sensitivity; all-surface baseline with unchanged material and scheduling assumptions.'
    receiver['peak_source_electric_W_at_R1e9']=max(s['source_rating']['peak_total_prepared_electric_W'] for s in receiver['stages'].values())
    receiver['peak_kinetic_source_electric_W_at_R1e9']=max(s['source_rating']['peak_kinetic_electric_W'] for s in receiver['stages'].values())
    receiver['thermal_radius_screens']=[]
    for rise in [500,1000,1500,2000]:
        rs={k:brentq(lambda r:s['material_peak_rise_K']*(1e9/r)**2.5+s['prompt_peak_rise_K']*(1e9/r)**2-rise,1e8,1e10) for k,s in receiver['stages'].items()}
        receiver['thermal_radius_screens'].append({'allowed_rise_K':rise,'minimum_R_m_in_surface_deposition_screen':max(rs.values()),'per_stage_R_m':rs})
    weighted=module('collector_scenarios','collector-scenarios.py').receiver_for_cycle(receiver,cycle)
    weighted['input_sha256']={p.name:c.sha(p) for p in [HERE/'process-scenarios.json',HERE/'collector.json',input_path]}
    weighted['scenario_input_sha256']=dict(weighted['input_sha256'])
    weighted['scenario_scope']='B2 pulse-level fields and per-helium stationary process weighting both refreshed for the specified refined family.'
    weighted['code_sha256'][str(Path(__file__).relative_to(c.ROOT))]=c.sha(Path(__file__))
    output_path=HERE/f'collector-{tag}.json';output_path.write_text(json.dumps(weighted,indent=2,allow_nan=False)+'\n')
    old_depth=json.loads((HERE/'collector-depth.json').read_text())
    photon_map={k:next(x['total_attenuation_depth']['prompt_peak_rise_K'] for x in v['monoenergetic_sensitivities'] if x['photon_energy_keV']==10) for k,v in old_depth['stages'].items()}
    depth=module('collector_depth','collector-depth.py');grid=np.linspace(0,t[-1],12001);dt=grid[1]-grid[0]
    q=np.diff(np.interp(grid,t,cumulative))/(dt*c.AREA)
    photon_map['c12']=float(max(depth.volume(q,dt,.2373*2000)))
    station={'cycle':cycle,'receiver':weighted,'selected_family':args.family,'selected_fuel_C_H_mass_ratio':float(args.target),
        'photon_depth_10keV_total_mu_K':photon_map,
        'source_sha256':{'process-scenarios.json':c.sha(HERE/'process-scenarios.json'),output_path.name:c.sha(output_path),'collector-depth.json':c.sha(HERE/'collector-depth.json')}}
    (HERE/f'station-inputs-{tag}.json').write_text(json.dumps(station,indent=2,allow_nan=False)+'\n')
    b=weighted['stages']['c12']
    print(json.dumps({'family':args.family,'B2':{k:b[k] for k in ['material_peak_rise_K','prompt_only_peak_rise_K','delayed_decay_peak_rise_K','sum_separate_maxima_K']},
      'source_peak_W':weighted['peak_source_electric_W_at_R1e9'],'occupied_s_per_He':weighted['occupied_s_per_He_at_R1e9'],
      'thermal_radius_screens':weighted['thermal_radius_screens'],'photon_map':photon_map},indent=2))
if __name__=='__main__':main()
