"""Bounded 250 kg/m2 extension of the frozen source128 finite-backing screen.
No new hydro or source retuning. Preserve the original module and outputs.
The explicit substitutions change thickness grid, matching time-resolution
check thickness, and destination only. All physical input hashes are retained.
"""
from pathlib import Path
import hashlib,json
P=Path(__file__).parent
original=P/'source128_thin_thermal.py'
text=original.read_text()
replacements={
    'for sigma in[500.,750.,1000.]:':'for sigma in[250.]:',
    "front_step(ee,.25,50,1800,'N')":"front_step(ee,.125,50,1800,'N')",
    "P/'source128-thin-thermal.json'":"P/'source128-receiver250-thermal.json'",
}
for before,after in replacements.items():
    assert text.count(before)==1,(before,text.count(before))
    text=text.replace(before,after)
namespace={'__file__':str(original),'__name__':'receiver250_frozen_source_extension'}
exec(compile(text,str(original),'exec'),namespace)
output=P/'source128-receiver250-thermal.json'
data=json.loads(output.read_text())
data['extension']={'wrapper_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'frozen_parent_sha256':hashlib.sha256(original.read_bytes()).hexdigest(),'substitutions':replacements,'scope':'Only the backing column is changed to250 kg/m2. Exact source128 profiles and faster2000/.99 cadence are retained as conservative loading; all-source thermal retention remains distinct from the coupled plant capture account.'}
output.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
