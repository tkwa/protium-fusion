"""Check passive-tracer quadrature and fixed-neutron-background sensitivity."""
import argparse,json,hashlib
from pathlib import Path
ap=argparse.ArgumentParser();ap.add_argument('--four',type=Path,required=True);ap.add_argument('--eight',type=Path,required=True);ap.add_argument('--zero-neutron',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--cycle',type=Path);a=ap.parse_args();four=json.loads(a.four.read_text());eight=json.loads(a.eight.read_text());zero=json.loads(a.zero_neutron.read_text());rows={}
for tag,s in eight['stages'].items():
 row={}
 for n,b in s['cold_tracer'].items():
  f=four['stages'][tag]['cold_tracer'][n];z=zero['stages'][tag]['cold_tracer'][n]
  keys=set(b['cold_output_nuclei_per_tag'])|set(f['cold_output_nuclei_per_tag'])|set(z['cold_output_nuclei_per_tag'])
  row[n]={'max_probability_change_4_to_8':max(abs(b['cold_output_nuclei_per_tag'].get(k,0)-f['cold_output_nuclei_per_tag'].get(k,0)) for k in keys),
  'max_probability_change_zero_neutron':max(abs(b['cold_output_nuclei_per_tag'].get(k,0)-z['cold_output_nuclei_per_tag'].get(k,0)) for k in keys),
  'noncycle_probability':b['noncycle_probability'],'noncycle_probability_zero_neutron':z['noncycle_probability'],
  'baryon_residual_per_tag':b['baryon_residual_per_tag']}
 rows[tag]={'native_reconstruction_L1':s['native_hot_reconstruction_L1_per_initial_catalyst'],'tracers':row}
checks=[{'name':'all individual runs pass conservation checks','passed':all(x['all_checks_pass'] for x in [four,eight,zero])}, {'name':'zero-neutron control enabled','passed':zero['zero_neutron_background_control'] and not eight['zero_neutron_background_control']}, {'name':'same implementation and rate coefficients','passed':len({x['code_sha256'] for x in [four,eight,zero]})==1 and len({x['rate_coefficients_sha256'] for x in [four,eight,zero]})==1}]
assert all(c['passed'] for c in checks)
out={'scope':__doc__,'stages':rows,'checks':checks,'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [a.four,a.eight,a.zero_neutron]},'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
if a.cycle:
 cy=json.loads(a.cycle.read_text())['perfect_eventual_return_cycle'];attempts=sum(cy['stage_attempts_per_He'].values());maximum=max(v['noncycle_probability'] for stage in eight['stages'].values() for v in stage['cold_tracer'].values());rows_p=[]
 for impurity in [.001,.01,.02]:
  extra=attempts*impurity*maximum;dCH=extra*12.011*1.6605390666e-27/cy['fresh_hydrogen_kg_per_He']
  rows_p.append({'impurity_number_fraction':impurity,'extra_nonquartet_nuclei_per_He':extra,'extra_C_H_mass_ratio':dCH,'fraction_remaining_0056_fuel_allowance':dCH/(.0056-cy['fresh_C_H_mass_ratio'])})
 out['purity_budget_screen']={'reference':'Original pulse family with perfect eventual ordinary-material return','stage_attempts_per_He':attempts,'largest_observed_passive_nonquartet_probability':maximum,'cases':rows_p,'scope':'Assign the largest observed tagged nonquartet probability to every impurity in every original stage, without subtracting replaced native leakage. This is a conservative comparison within the fixed-history calculation, not a physical upper bound after energy/composition feedback. Purity is by catalyst number, not total fuel mass.'};out['source_sha256'][a.cycle.name]=hashlib.sha256(a.cycle.read_bytes()).hexdigest()
a.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({tag:{'native_L1':s['native_reconstruction_L1'],'max_subdivision_change':max(v['max_probability_change_4_to_8'] for v in s['tracers'].values()),'max_neutron_change':max(v['max_probability_change_zero_neutron'] for v in s['tracers'].values())} for tag,s in rows.items()},indent=2))
