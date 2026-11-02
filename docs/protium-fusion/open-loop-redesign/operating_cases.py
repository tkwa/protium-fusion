"""Rounded operating cases, distinct from the finite search's boundary maxima."""
import hashlib
import json
from pathlib import Path
from station import evaluate, YEAR

P=Path(__file__).resolve().parent
inputs=json.loads((P/'station-inputs.json').read_text())
c=inputs['cycle']; r=inputs['receiver']
depth=json.loads((P/'collector-depth.json').read_text())
ten={k:next(x for x in v['monoenergetic_sensitivities'] if x['photon_energy_keV']==10)['total_attenuation_depth']['prompt_peak_rise_K'] for k,v in depth['stages'].items()}
rows={}
for name,R,rating,service,surface,soft in [
    ('balanced_surface',1.1e9,10.,5000.,'surface',0.),
    ('compact_service_surface',1.1e9,10.,1000.,'surface',0.),
    ('high_source_rating_surface',1.1e9,100.,1000.,'surface',0.),
    ('spectral_sensitivity',.7e9,100.,1000.,'finite_depth',.1),
]:
    f=.8/(r['occupied_s_per_He_at_R1e9']['mass999']*R/1e9)
    d=evaluate(f,R,c,r,absorber_sigma=5000,driver_peak_W_kg=rating,service_sigma=service,
               surface_model=surface,soft_surface_fraction=soft,photon_depth_K=ten)
    d.update({
        'source_and_service_status':'Specified hardware requirements; no component-resolved mass design.',
        'receiver_refurbishment_interpretation':'Minimum receiver carbon return if it alone uses the entire structural carbon allowance; other carbon losses tighten this requirement.',
        'net_electric_J_per_fresh_kg_H':d['net_W']/d['feed_H_kg_s'],
        'net_electric_J_per_total_kg_C_budget':d['net_W']/(.007*d['feed_H_kg_s']),
        'gross_manufacturing_energy_J_at100MJkg':d['mass_kg']*1e8,
        'one_way_shipment_kinetic_J_at30kms':.5*d['mass_kg']*(3e4)**2,
        'solar_30pct_1kgm2_equal_specific_power_AU':(1361*.3/d['net_W_kg'])**.5,
        'replacement_energy_budget_fraction_of_net':d['replacement_W']/d['net_W'],
        'pulse_intervals_days':{k:v/86400 for k,v in d['pulse_intervals_s'].items()},
        'thermal_mass_total_capacity_J_K':d['mass_breakdown_kg']['receiver']*1000,
        'mean_recovered_heat_flux_W_m2':d['gross_W']/.4/(4*3.141592653589793*R**2),
    })
    rows[name]=d
out={'scope':'Rounded conditional operating points. The spectral sensitivity assigns 10% of prompt energy to mathematical-surface deposition and90% to a10keV attenuation-depth model; the escaping spectrum is uncomputed.',
     'inputs_sha256':hashlib.sha256((P/'station-inputs.json').read_bytes()).hexdigest(),
     'station_code_sha256':hashlib.sha256((P/'station.py').read_bytes()).hexdigest(),
     'depth_input_sha256':hashlib.sha256((P/'collector-depth.json').read_bytes()).hexdigest(),
     'cases':rows}
(P/'operating-cases.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
print(json.dumps({k:{x:v[x] for x in ['feasible','net_W_kg','net_W','mass_kg','peak_surface_K']} for k,v in rows.items()},indent=2))
