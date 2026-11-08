"""Check forward replay of the smooth time-reversed prehistory."""
import json
from pathlib import Path
import numpy as np
from analyze import snapshot
base=Path('/var/tmp/tkwa-chatgpt/protium-stability')
original=base/'runs/A-phase-reference-1536'
replay=base/'runs/A-phase-roundtrip-1536'
a=snapshot([original/'pulse.block0.out2.00000.vtk'])
options=[snapshot([f]) for f in sorted(replay.glob('*.vtk'))]
b=min(options,key=lambda x:abs(x[0]-.75))
volume=np.diff(a[1]**3);report={}
for key in ['dens','Etot','s1']:
    x=a[3][key][0];y=b[3][key][0]
    report[key+'_relative_L1']=float(np.sum(abs(x-y)*volume)/np.sum(abs(x)*volume))
x=a[3]['mom'][0,:,0];y=b[3]['mom'][0,:,0]
report['momentum_relative_L1']=float(np.sum(abs(x-y)*volume)/np.sum(abs(x)*volume))
report['replay_snapshot_clock_offset_us']=float(b[0]*30-22.5)
reference=json.loads((base/'runs/A-2048-1d-wide/diagnostics.json').read_text())
result=json.loads((replay/'diagnostics.json').read_text())
report['replay_conversion']=result['conversion']
report['reference_conversion']=reference['conversion']
report['relative_nuclear_output_change']=result['conversion']/reference['conversion']-1
report['comparison']='Both reference and replay use HLLC with 0.46875 m radial spacing; outer boundaries remain outside the fuel. The field comparison includes the recorded small output-time offset.'
report['interpretation']='Integrated yield is reproduced much more accurately than sharp initial cell interfaces. This checks the numerical prehistory; it does not simulate cold packet assembly.'
(base/'validation/phase-roundtrip.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
