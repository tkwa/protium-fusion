"""Manufacturing-energy sensitivity at fixed installed station hardware.

All 100 MJ/kg baseline manufacturing capacity, radiators, reserve and support
mass are retained. Only the recurring manufacturing electricity changes.
The rows are effective TOTAL process costs, not claimed chemical efficiencies.
"""
import hashlib
import json
from pathlib import Path
P=Path(__file__).resolve().parent
OLD=P.parent/'open-loop-redesign'/'operating-cases.json'
b=json.loads(OLD.read_text())['cases']['high_source_rating_surface']
YEAR=365.25*86400
mass=b['mass_kg']; flow=mass*b['replacement_fraction_per_year']/YEAR
rows=[]
for mj in [.1,1.,1.5,2.,3.1,10.,100.]:
    e=mj*1e6
    manufacturing=flow*e
    net=b['net_W']+b['replacement_W']-manufacturing
    rows.append({'effective_total_MJ_kg':mj,'routine_manufacturing_W':manufacturing,
        'net_W':net,'net_W_kg':net/mass,
        'construction_J_if_same_effective_total_cost':mass*e,
        'construction_energy_repayment_years_if_same_cost':mass*e/net/YEAR,
        'power_change_relative_residual':(net-b['net_W']-(b['replacement_W']-manufacturing))/net})
assert abs(rows[-1]['net_W']/b['net_W']-1)<1e-14
assert all(abs(r['power_change_relative_residual'])<1e-14 for r in rows)
out={'scope':__doc__,'mass_kg':mass,'refurbishment_mass_kg_s':flow,
    'fresh_C_H_mass_budget':.007,'carbon_budget_unchanged':True,
    'user_premise':{'minimum_efficiency_for_usable_energy_to_higher_energy_product':.5,
        'higher_efficiency_allowed':True,
        'definition':'Useful chemical free-energy increase divided by supplied electrical or usable supplementary chemical energy, with feed/product states specified.'},
    'rows':rows,'limit_if_routine_manufacturing_electricity_vanishes_W_kg':(b['net_W']+b['replacement_W'])/mass,
    'limitations':['No component manufacturing recipes or material-weighted chemical work have been solved.',
        'A dissipation estimate excludes energy retained in higher-energy products and is not automatically a total process cost.',
        'The old manufacturing hardware mass and throughput remain allocated; lower electricity does not grant a smaller factory.',
        'Other plant loads, source/receiver qualifications, transport, operating material return and membrane lifetime remain separate.'],
    'input_sha256':hashlib.sha256(OLD.read_bytes()).hexdigest(),
    'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(P/'sensitivity.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
print(json.dumps({'rows':rows,'zero_cost_limit_W_kg':out['limit_if_routine_manufacturing_electricity_vanishes_W_kg']},indent=2))
