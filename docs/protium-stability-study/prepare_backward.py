"""Construct inputs for nonreacting expansion used to extend a smooth inflow backward."""
import argparse,hashlib,json,shutil
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('inputs',type=Path);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
a.out.mkdir(parents=True,exist_ok=True)
for name in ['eos.bin','rates.txt']:shutil.copy2(a.inputs/name,a.out/name)
lines=(a.inputs/'profile.txt').read_text().splitlines();rows=[lines[0]]
for line in lines[1:]:
    values=list(map(float,line.split()));values[-2]*=-1;values[-1]*=-1
    rows.append(' '.join(f'{v:.17g}' for v in values))
(a.out/'profile.txt').write_text('\n'.join(rows)+'\n')
manifest=dict(purpose='Velocity-reversed nonreacting expansion; reverse its output momenta and time labels to construct a smooth incoming prehistory',
              parent=str(a.inputs),profile_sha256=hashlib.sha256((a.out/'profile.txt').read_bytes()).hexdigest())
(a.out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
