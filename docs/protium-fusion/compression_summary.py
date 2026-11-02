"""Summarize saved convergent-flow experiments and sampled density extrema."""
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).parent
rows=[]
archive=ROOT/'compression-traces.json'
compact=json.loads(archive.read_text())['traces'] if archive.exists() else {}
for path in sorted(ROOT.glob('compression-*-burning.json')):
    d=json.loads(path.read_text());trace=path.with_suffix('.npz')
    nz=d['zones']
    if path.stem in compact:
        peak_density=compact[path.stem]['sampled_peak_density_kg_m3']
    elif trace.exists():
        a=np.load(trace);nv=int(a['variables_per_zone']);nn=len(a['species'])
        z=a['state'].reshape(nz,nv,-1)
        r=np.vstack([np.zeros(z.shape[-1]),z[:,nn+1]*d['R0_m']])
        mass=d['mass_kg']*np.diff(np.linspace(0,1,nz+1)**3)
        rho=mass[:,None]/(4*np.pi/3*np.diff(r**3,axis=0))
        peak_density=float(rho.max())
    else: continue
    rows.append({'file':path.name,'zones':nz,
      'conversion':d['desired_product_per_initial_catalyst'],
      'nuclear_over_fuel_input':d['deposited_nuclear_over_initial_energy'],
      'initial_kinetic_fraction':d['initial_kinetic_energy_J']/d['initial_energy_J'],
      'sampled_peak_density_kg_m3':peak_density,
      'peak_temperature_K':d['peak_temperature_K'],
      'radiated_fraction_nuclear':d['radiation_escaped_J']/d['nuclear_energy_J'],
      'energy_residual_fraction_input':d['energy_residual_fraction_initial']})
out={'scope':'Specified convergent flow; input counts its kinetic and internal energy. Its preparation, driver efficiency, multidimensional stability and plant mass remain unmodeled.',
     'rows':rows}
(ROOT/'compression-summary.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
