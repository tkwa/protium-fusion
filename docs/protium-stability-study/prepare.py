"""Prepare immutable inputs for a bounded Athena++ stability experiment.

SI units externally; code units: 30 m, 1e8 kg/m3, 1e6 m/s, 30 microseconds.
The first model freezes the ion number per baryon, retains electron/positron
and LTE radiation thermodynamics, and traps all radiation. This explicitly
differs from the composition-dependent, diffusive reference model.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path
import numpy as np

p = argparse.ArgumentParser()
p.add_argument('--root', type=Path, required=True)
p.add_argument('--out', type=Path, required=True)
p.add_argument('--yi', type=float, default=2/15)
p.add_argument('--nrho', type=int, default=513)
p.add_argument('--ntemp', type=int, default=769)
a = p.parse_args()
a.out.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(a.root))
from fast_eos import FastEOS
import pulse_model as pm
eos = FastEOS()
# Stay inside the published source table, with a small endpoint margin.
lr = np.linspace(eos.ld[0]+1e-8-np.log(8/15)+np.log(1000/1e8),
                 eos.ld[-1]-1e-8-np.log(8/15)+np.log(1000/1e8), a.nrho)
lt = np.linspace(eos.lt[0]+1e-8, eos.lt[-1]-1e-8, a.ntemp)
r, t = np.meshgrid(np.exp(lr), np.exp(lt), indexing='ij')
s = eos.state(r.ravel()*1e5, t.ravel(), 8/15, a.yi)
# P and rho*e in code units; gamma1 is dimensionless.
data = np.stack([np.log(s.p/1e21), np.log(s.e*r.ravel()*1e-16), s.gamma1], -1)
with (a.out/'eos.bin').open('wb') as f:
    np.array([a.nrho,a.ntemp], dtype='<i4').tofile(f)
    lr.astype('<f8').tofile(f); lt.astype('<f8').tofile(f)
    data.astype('<f8').tofile(f)

source = a.root/'reassessment-2049/small-stages/scaled-hot-A-64-R30-initial.json'
initial = json.loads(source.read_text())
radii = np.r_[0, initial['outer_radius_m']]/30
vel = np.r_[0, initial['velocity_m_s']]/1e6
mass = 4*np.pi/3*np.diff(np.linspace(0,1,65)**3)
rho = mass/(4*np.pi/3*np.diff(radii**3))
temp = np.array(initial['T_K'])
with (a.out/'profile.txt').open('w') as f:
    f.write('64\n')
    for i in range(64):
        f.write(' '.join(f'{x:.17g}' for x in (radii[i],radii[i+1],rho[i],temp[i],vel[i],vel[i+1]))+'\n')
rates = [pm.LIB.get_rate_by_name(n) for n in ['n14(p,g)o15','o15(,p)n14']]
q = (pm.MASS[pm.SPECIES.index('p')]+pm.MASS[pm.SPECIES.index('n14')]-pm.MASS[pm.SPECIES.index('o15')])*pm.MEV/(15*pm.MU)
with (a.out/'rates.txt').open('w') as f:
    f.write(f'{q/1e12:.17g}\n')
    for rate in rates:
        f.write(f'{len(rate.sets)}\n')
        for rs in rate.sets: f.write(' '.join(f'{c:.17g}' for c in rs.a)+'\n')
manifest = dict(yi=a.yi, ye=8/15, q_J_kg=q, length_m=30, density_kg_m3=1e8,
    velocity_m_s=1e6, time_s=3e-5, source=str(source),
    source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
    eos_source_sha256=hashlib.sha256((a.root/'electron-eos.npz').read_bytes()).hexdigest(),
    eos_shape=[a.nrho,a.ntemp], input_hashes={x.name:hashlib.sha256(x.read_bytes()).hexdigest()
        for x in a.out.iterdir() if x.suffix in ['.bin','.txt']})
(a.out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps(manifest,indent=2))
