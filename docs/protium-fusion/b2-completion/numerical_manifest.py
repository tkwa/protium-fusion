"""Hash and verify the numerical owner's research artifacts without rewriting runs."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

ap=argparse.ArgumentParser()
ap.add_argument('--packet',type=Path,required=True)
ap.add_argument('--output',type=Path,required=True)
a=ap.parse_args()
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
names=[
 'NUMERICAL-COMPLETION.md','radial_completion.py','radial_completion-v1.py',
 'run.py','run-v1.py','run-v2.py','rhs_audit.py','rhs-audit.json',
 'viscosity_checks.py','viscosity-checks.json','archive_outputs.py',
 'completion_summary.py','completion-summary.json','completion_plot.py',
 'completion.png','completion.svg','numerical_manifest.py','runtime.json',
]
checks=[]
cases=[]
for reference in sorted(a.packet.glob('B2-*-trace-reference.json')):
 r=json.loads(reference.read_text());name=r['result'].removesuffix('.json')
 result=a.packet/r['result'];compact=a.packet/r['compact_trace'];full=Path(r['full_trace_path'])
 checks.append({'case':name,'metadata_hash':sha(result)==r['result_sha256'],
  'compact_hash':sha(compact)==r['compact_trace_sha256'],
  'full_trace_hash':sha(full)==r['full_trace_sha256']})
 with np.load(compact) as small,np.load(full) as large:
  selected=small['full_trace_sample_indices']
  checks[-1]['exact_state_samples']=bool(np.array_equal(small['state'],large['state'][:,selected]))
  checks[-1]['exact_time_samples']=bool(np.array_equal(small['time_s'],large['time_s'][selected]))
  checks[-1]['endpoints_retained']=bool(selected[0]==0 and selected[-1]==len(large['time_s'])-1)
 for suffix in ['.json','.initial.json','.inputs.json','.trace.npz','-trace-reference.json']:
  if (a.packet/(name+suffix)).exists():names.append(name+suffix)
 rr=json.loads(result.read_text());cases.append({'name':name,'success':rr['success'],
  'termination_reason':rr['termination_reason'],'elapsed_sound_times':rr['elapsed_sound_times']})
all_pass=all(all(v for k,v in row.items() if k!='case') for row in checks)
manifest={'scope':'Numerical implementation/run artifacts only; independent audits and supporting-stage regressions have their own provenance.',
 'files':{n:{'sha256':sha(a.packet/n),'bytes':(a.packet/n).stat().st_size} for n in sorted(set(names)) if (a.packet/n).exists()},
 'cases':cases,'trace_integrity_checks':checks,'all_trace_integrity_checks_pass':all_pass}
a.output.write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps({'artifact_count':len(manifest['files']),'case_count':len(cases),'all_trace_integrity_checks_pass':all_pass}))
assert all_pass
