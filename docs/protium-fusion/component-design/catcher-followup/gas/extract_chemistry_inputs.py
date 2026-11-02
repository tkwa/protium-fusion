"""Extract atomic masses and exact fully-decayed nitrogen flux from final350."""
from pathlib import Path
import sys,json,hashlib
P=Path(__file__).resolve().parent
R=Path('/home/tkwa/.codex/worktrees/protium-fusion/website/docs/protium-fusion')
sys.path.insert(0,str(R/'component-design'))
import material_cycle as mc
d=json.loads((P/'inputs.json').read_text())
mass={n:mc.atomic_kg(n) for n in ['p','he4','c12','c13','n14','n15']}
N={cat:{sp:s['cold_output_counts'].get(sp,0)*s['attempts_per_pulse']*mass[sp] for sp in ['n14','n15']} for cat,s in d['stages'].items()}
flow={sp:sum(z[sp]/d['plant']['pulse_interval_s'][cat] for cat,z in N.items()) for sp in ['n14','n15']}
out={'atomic_mass_kg':mass,'stage_fully_decayed_N_output_atomic_kg':N,'mean_fully_decayed_N_output_atomic_kg_s':flow,'source_sha256':{'inputs.json':hashlib.sha256((P/'inputs.json').read_bytes()).hexdigest(),'material_cycle.py':hashlib.sha256((R/'component-design/material_cycle.py').read_bytes()).hexdigest()},'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'Total material-ledger nitrogen after decays. A fraction released asgas is a separate parameter; this does not predict chemical liberation from the catcher.'}
(P/'chemistry-inputs.json').write_text(json.dumps(out,indent=2)+'\n')
print(flow)
