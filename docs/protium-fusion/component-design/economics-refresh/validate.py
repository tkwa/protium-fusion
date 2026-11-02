"""Regression and physical scaling checks for the final-point economics interface."""
from pathlib import Path
import copy
import hashlib
import json
import math
import sys
from calculate import evaluate, render, select_pointer, sha

P=Path(__file__).resolve().parent
old_in=json.loads((P.parent/'economics-inputs.json').read_text())
old_out=json.loads((P.parent/'economics.json').read_text())
checks=[]
def check(name,value):
    checks.append({'name':name,'passed':bool(value)})
def close(a,b):return math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-15)
for name,record in old_in['cases'].items():
    new=evaluate(record['plant'],record['cycle']);old=old_out['cases'][name]
    check(name+': all ledger identities',new['validation_summary']['all_passed'])
    for value,reference in [(new['plant']['installed_kg_per_net_W'],old['installed_kg_per_net_W']),
       (new['embodied_energy']['payback_years_per_MJ_kg_at_declared_availability'],old['whole_plant_embodied_scenarios'][0]['equivalent_export_years']),
       (new['receiver_chemical_step']['receiver_C_kg'],old['receiver_C_kg']),
       (new['acquisition']['associated_resource_H_kg_at_complete_C_recovery'],old['resource_H_kg_accessed_for_receiver_at_chosen_C_H']),
       (new['receiver_chemical_step']['scenarios'][0]['equivalent_export_years'],old['methane_chemical_step'][0]['receiver_chemical_step_equivalent_export_years']),
       (new['solar_comparison']['mass_specific_power_parity_au'],old['direct_solar_mass_parity_distance_au'])]:
        check(name+': archived numerical regression '+str(len(checks)),close(value,reference))
    check(name+': generated report names exact power',f"{record['plant']['net_W']:.6g}" in render(new))
base=next(iter(old_in['cases'].values()))
original=copy.deepcopy(base)
x=evaluate(base['plant'],base['cycle'],{'availability':.5})
y=evaluate(base['plant'],base['cycle'])
check('half availability doubles energy payback',close(x['embodied_energy']['payback_years_per_MJ_kg_at_declared_availability'],2*y['embodied_energy']['payback_years_per_MJ_kg_at_declared_availability']))
check('half availability preserves installed mass and nominal export',x['plant']['mass_kg']==y['plant']['mass_kg'] and x['plant']['net_W']==y['plant']['net_W'])
check('evaluation leaves caller input intact',base==original)
p=copy.deepcopy(base['plant']);p['net_W']*=2;p['net_W_kg']*=2
z=evaluate(p)
check('double net export halves embodied payback',close(z['embodied_energy']['payback_years_per_MJ_kg_at_declared_availability'],.5*y['embodied_energy']['payback_years_per_MJ_kg_at_declared_availability']))
check('double net export leaves resource inventory unchanged',z['acquisition']==y['acquisition'])
w=evaluate(base['plant'],assumptions={'receiver_carbon_mass_fraction':.5})
check('half carbon fraction halves methane energy',close(w['receiver_chemical_step']['scenarios'][0]['supplied_total_J'],.5*y['receiver_chemical_step']['scenarios'][0]['supplied_total_J']))
check('half carbon fraction leaves reactor power unchanged',w['plant']['net_W']==y['plant']['net_W'])
minimal={'mass_kg':100.,'net_W':2.,'mass_breakdown_kg':{'receiver':40.,'services':60.}}
v=evaluate(minimal)
check('20 mW/kg scaling point',close(v['plant']['specific_power_W_kg'],.02))
check('minimal ledger does not invent feed/repair data',v['operating_materials'] is None and v['refurbishment'] is None and v['hydrogen_processing'] is None)
check('mass ranking updates to services',v['plant']['mass_items_descending'][0]=='services')
for label,point,assumptions in [('inconsistent mass',{**minimal,'mass_kg':101},None),('inconsistent specific power',{**minimal,'net_W_kg':.03},None),('infeasible plant',{**minimal,'feasible':False},None),('efficiency below user premise',minimal,{'methane_chemical_efficiencies':[.49]}),('invalid availability',minimal,{'availability':0})]:
    failed=False
    try:evaluate(point,assumptions=assumptions)
    except ValueError:failed=True
    check('reject '+label,failed)
check('explicit nested JSON selection',select_pointer({'rows':[{'plant':minimal}]},'/rows/0')['plant']==minimal)
check('JSON pointer escaping',select_pointer({'a/b':{'~x':2}},'/a~1b/~0x')==2)
report={'scope':'Arithmetic, physical-scaling and archived-point regression checks; no reactor simulation or money-price validation.',
        'code_sha256':sha(P/'calculate.py'),'validation_code_sha256':sha(__file__),
        'prior_input_sha256':sha(P.parent/'economics-inputs.json'),'prior_output_sha256':sha(P.parent/'economics.json'),
        'checks':checks,'passed':sum(c['passed'] for c in checks),'count':len(checks),'all_passed':all(c['passed'] for c in checks),
        '20mW_kg_scaling_only':{'scope':'Unit scaling illustration with arbitrary mass; not a selected reactor design.',
             'years_per_MJ_kg':v['embodied_energy']['payback_years_per_MJ_kg_at_declared_availability'],
             '20year_energy_MJ_kg':v['embodied_energy']['lifetime_acceptance'][0]['maximum_average_total_J_kg_for_energy_payback']/1e6,
             'solar_parity_au_at_declared_comparison':v['solar_comparison']['mass_specific_power_parity_au']}}
assert report['all_passed']
(P/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:report[k] for k in ['passed','count','all_passed','20mW_kg_scaling_only']},indent=2))
