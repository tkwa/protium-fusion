"""Compact trace diagnostics; sampled density extrema are not exact extrema."""
from pathlib import Path
import json
import sys
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1];sys.path.insert(0,str(ROOT))
from fast_eos import FastEOS

def main():
    eos=FastEOS(HERE/'electron-eos-cold-surrogate.npz')
    rows=[]
    for filename in sorted(HERE.glob('cold-driver-A-*.npz')):
        summary=json.loads(filename.with_suffix('.json').read_text())
        initial=json.loads(filename.with_name(filename.stem+'-initial.json').read_text())
        d=np.load(filename);nz=int(d['zones']);nv=int(d['variables_per_zone'])
        state=d['state'].reshape(nz,nv,-1)
        r=state[:,4,:]*300;T=state[:,3,:]*3e8
        m=summary['mass_kg']*np.diff(np.linspace(0,1,nz+1)**3)
        rho=m[:,None]/(4*np.pi/3*np.diff(np.vstack([np.zeros((1,r.shape[1])),r])**3,axis=0))
        where=np.unravel_index(np.argmax(rho),rho.shape)
        frac=np.sum((rho>=1e8)*m[:,None],axis=0)/m.sum()
        rho0=np.array(initial['rho_kg_m3']);temp0=np.array(initial['T_K'])
        vel0=np.r_[0.,initial['velocity_m_s']]
        s=eos.state(rho0/1000,temp0,8/15,2/15)
        cs=np.sqrt(s.gamma1*s.p/(rho0/1000))*.01
        dv=np.maximum(-np.diff(vel0),0.)
        q=rho0*(summary['artificial_viscosity_quadratic']*dv**2+
                 summary['artificial_viscosity_linear']*cs*dv)
        ratio=q/(s.p*.1)
        rows.append({'file':filename.name,'zones':nz,'initial_temperature_K':float(temp0[0]),
                     'initial_max_density_kg_m3':float(rho0.max()),
                     'nuclear_enabled':summary['nuclear_enabled'],
                     'viscosity_linear':summary['artificial_viscosity_linear'],
                     'viscosity_quadratic':summary['artificial_viscosity_quadratic'],
                     'termination':summary['termination_reason'],
                     'conversion_at_stop':summary['desired_product_per_initial_catalyst'],
                     'nuclear_over_input_at_stop':summary['deposited_nuclear_over_initial_energy'],
                     'exact_saved_peak_temperature_K':summary['peak_temperature_K'],
                     'sampled_peak_density_kg_m3':float(rho[where]),
                     'time_of_sampled_peak_density_s':float(d['time_s'][where[1]]),
                     'temperature_at_sampled_density_peak_K':float(T[where]),
                     'maximum_sampled_mass_fraction_above_1e8_kg_m3':float(frac.max()),
                     'minimum_sampled_outer_radius_m':float(np.min(r[-1])),
                     'peak_sampled_mass_weighted_density_kg_m3':float(np.max(np.sum(rho*m[:,None],axis=0)/m.sum())),
                     'initial_artificial_viscosity_to_pressure_mass_weighted':float(np.dot(m,ratio)/m.sum()),
                     'initial_artificial_viscosity_to_pressure_maximum':float(ratio.max()),
                     'energy_residual_fraction_initial':summary['energy_residual_fraction_initial'],
                     'wall_seconds':summary['wall_seconds']})
    failures=[json.loads(p.read_text()) for p in sorted(HERE.glob('*-failure.json'))]
    out={'scope':'Diagnostics of low-temperature surrogate compression. Nuclear values at temperature-guard termination are incomplete; no final yield is inferred.',
         'sampling':'At most 160 saved time samples; reported density and radius extrema may miss brief extrema.',
         'runs':rows,'numerical_failures':failures}
    (HERE/'cold-driver-summary.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    for r in rows:print(r['file'],r['termination'],'rho_max',r['sampled_peak_density_kg_m3'],'gain@stop',r['nuclear_over_input_at_stop'],'q/P initial',r['initial_artificial_viscosity_to_pressure_mass_weighted'])

if __name__=='__main__':main()
