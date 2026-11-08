#!/usr/bin/env python3
"""Export pinned rates and compact source trajectories; existing research is read-only."""
import gzip
import hashlib
import json
from pathlib import Path
import numpy as np
import pynucastro as p

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RESEARCH = Path('/home/tkwa/src/protium-fusion/output')
sha = lambda f: hashlib.sha256(Path(f).read_bytes()).hexdigest()

def main():
    old = ['n','p','he4','c12','c13','n13','n14','n15','o14','o15','o16',
           'o17','o18','f17','f18','f19','ne18','ne19','ne20','ne21','ne22',
           'na21','na22','na23','mg22','mg23','mg24']
    extra = ['d','t','he3','li7','be7','be8','b8','b10','b11','c11','c14',
             'n12','n16','o13','o19','f16','f20','ne17','ne23','ne24',
             'na19','na20','na24','na25','mg20','mg21','mg25','mg26',
             'al22','al23','al24','al25','al26','al27','si24','si25','si26','si27','si28']
    names = old + extra
    nucs = [p.Nucleus(n) for n in names]
    lib = p.ReacLibLibrary()
    rates = lib.linking_nuclei(nucs, with_reverse=True).get_rates()
    out = dict(version=p.__version__, species=names, A=[n.A for n in nucs],
               Z=[n.Z for n in nucs], mass_MeV=[n.mass for n in nucs], rates=[])
    for r in rates:
        out['rates'].append(dict(name=r.fname, label=r.src, Q_MeV=r.Q,
            react=[nucs.index(n) for n in r.reactants],
            product=[nucs.index(n) for n in r.products], prefactor=r.prefactor,
            density_exponent=r.dens_exp, weak=bool(r.weak_type),
            ye_weighting=r.use_ye_weighting, coefficients=[s.a for s in r.sets]))
    (HERE/'rates.json').write_text(json.dumps(out,indent=1)+'\n')
    inp=ROOT/'docs/protium-fusion/component-design/final-family-dispatch/sigma350/input.json'
    family=json.loads(inp.read_text())
    stages={}
    for stage, s in family['loaded_stages'].items():
        path=Path(s['source']); d=json.loads(path.read_text()); z=np.load(d['trace_path'])
        m=z['cell_mass_kg']; time=z['time_s']; nv=d['state_variables']
        state=z['state'].reshape(len(m),int(z['variables_per_zone']),len(time))
        radii=np.vstack([np.zeros(len(time)),state[:,nv.index('r/R0')]*d['R0_m']])
        rho=m[:,None]/(4*np.pi/3*np.diff(radii**3,axis=0))
        stages[stage]=dict(time_s=time.tolist(), mass_fraction=(m/m.sum()).tolist(),
            temperature_K=(state[:,nv.index('T/T0')]*d['T0_K']).tolist(),
            density_g_cm3=(rho/1000).tolist(),
            light_Y=state[:,:3,:].tolist(),
            electron_Y=np.einsum('i,zit->zt',np.array(out['Z'][:27]),state[:,:27,:]).tolist(),
            initial_Y=d['initial_Y'],
            source_record=str(path),record_sha256=sha(path),
            trace_sha256=sha(d['trace_path']), mass_kg=d['mass_kg'],
            prepared_J=d['initial_energy_J'])
        recv=json.loads((path.parent/(stage+'.receiver.json')).read_text())
        stages[stage]['ballistic']=recv['ballistic']
        stages[stage]['source_rating']=recv['source_rating']
    with gzip.open(HERE/'trajectories.json.gz','wt') as f: json.dump(stages,f,separators=(',',':'))
    files=[RESEARCH/'second-generation-tamper/TAMPER-STUDY.md',
           RESEARCH/'larger-pulse-screen/README.md',
           RESEARCH/'fuel-assembly-run/radial-summary.json',
           RESEARCH/'reactive-b1/launcher_budget.json',
           RESEARCH/'reactive-b1/amr-benchmark-results.json']
    evidence={str(f):dict(sha256=sha(f),content=f.read_text()) for f in files}
    (HERE/'newer-research-snapshot.json').write_text(json.dumps(evidence,indent=2)+'\n')
    print(json.dumps(dict(species=len(names),rates=len(rates),stages=list(stages))))

if __name__=='__main__': main()
