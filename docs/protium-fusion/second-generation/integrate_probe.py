#!/usr/bin/env python3
"""Close actual quartet material/energy accounts for completed N15 probes."""
import hashlib
import json
from pathlib import Path
import sys
P=Path(__file__).resolve().parent; F=P.parent
sys.path.insert(0,str(F/'component-design'))
import material_cycle as mc

def main():
    base=json.loads((F/'component-design/final-family-dispatch/sigma350/input.json').read_text())
    rows=[]
    for path in sorted(P.glob('pulse-n15-*.json')):
        if '.pretrace.' in path.name:continue
        d=json.loads(path.read_text())
        if not d.get('success'):continue
        stages=dict(base['loaded_stages']); stages['n15']=mc.load_stage(path)
        new=mc.solve_cycle(stages)
        delta=(base['cycle']['U_prepared_MeV_per_He']-base['cycle']['U_hot_burn_start_MeV_per_He'])/base['cycle']['stage_attempts_per_He']['n15']
        new['U_hot_burn_start_MeV_per_He']=new['U_prepared_MeV_per_He']
        prep=delta*new['stage_attempts_per_He']['n15']
        new['U_prepared_MeV_per_He']+=prep;new['preparation_escaped_MeV_per_He']=prep
        rows.append(dict(record=path.name,record_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                         probe=d['probe'],cycle=new,stage=stages['n15']))
    (P/'probe-cycles.json').write_text(json.dumps(dict(scope=__doc__,cases=rows,
        qualification='Source128 A/C12/C13 unchanged. New N15 hot profiles are simulated; added precursor radiation is a linear per-attempt sensitivity. Full cold preparation and 3D assembly of the larger pulse remain unmodeled.'),indent=2)+'\n')
    for r in rows:
        c=r['cycle'];print(r['record'],'U',c['U_prepared_MeV_per_He'],'H',c['H_nonnu_MeV_per_He'],
            'N15 attempts',c['stage_attempts_per_He']['n15'],'C/H',c['fresh_C_H_mass_ratio'])

if __name__=='__main__':main()
