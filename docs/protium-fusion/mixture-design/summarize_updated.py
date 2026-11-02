from pathlib import Path
import json,hashlib,argparse
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
a=argparse.ArgumentParser();a.add_argument('--runs',type=Path,required=True);a.add_argument('--output',type=Path,required=True);a=a.parse_args();rows=[]
for p in sorted(a.runs.glob('*.json')):
 if any(s in p.name for s in ('.inputs.','.initial.','.pretrace.')) or p.name.endswith('status.json'):continue
 d=json.loads(p.read_text())
 if 'rate_update' not in d:continue
 b=json.loads(Path(d['rate_update']['baseline_record']).read_text())
 rows.append({'case':p.stem,'source_record':str(p.resolve()),'record_sha256':sha(p),'initial_profile_sha256':sha(p.with_suffix('.initial.json')),'trace_sha256':sha(p.with_suffix('.npz')),'baseline_record':d['rate_update']['baseline_record'],'revision':d['rate_update']['revision'],'success':d['success'],'termination_reason':d['termination_reason'],'requested_duration_s':d['mixture_inputs']['duration_s'],'elapsed_s':d['elapsed_s'],'initial_energy_relative_change':d['initial_energy_J']/b['initial_energy_J']-1,'old_gain':b['deposited_nuclear_over_initial_energy'],'new_gain':d['deposited_nuclear_over_initial_energy'],'gain_relative_change':d['deposited_nuclear_over_initial_energy']/b['deposited_nuclear_over_initial_energy']-1,'old_target_yield':b['desired_product_per_initial_catalyst'],'new_target_yield':d['desired_product_per_initial_catalyst'],'energy_residual_fraction_initial':d['energy_residual_fraction_initial'],'output_mass_fractions':d['output_mass_fractions']})
a.output.write_text(json.dumps({'scope':'Exact-profile full hydro rate comparisons; endpoints retained for coupled-cycle routing.','cases':rows},indent=2)+'\n')
for d in rows:print(d['case'].split('-fixed')[0],f"Q {d['old_gain']:.9f} -> {d['new_gain']:.9f}; yield {d['old_target_yield']:.8f} -> {d['new_target_yield']:.8f}; residual {d['energy_residual_fraction_initial']:.3g}")
