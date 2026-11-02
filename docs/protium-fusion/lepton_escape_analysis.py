"""Analyze finite-length escape and parameterized electron-supply allowances.

The kinetic model retains unmagnetized encounter coefficients. The supply
table prescribes an all-in electrical cost per electron; it does not equate
radiation in the local drift frame with that cost or predict a gas pressure.
"""
from pathlib import Path
import argparse,json,gzip
import numpy as np

ROOT=Path(__file__).parent
KB=1.380649e-23
MH=1.6735575e-27


def analyze(root,outdir):
    names=['half','main','double','closed','closed-fine']
    data={name:json.loads((outdir/f'lepton-escape-{name}.json').read_text()) for name in names}
    L=data['main']['active_length_m']
    x=np.array([.5,1.,2.])
    times=np.array([data[n]['residence_s']['mean'] for n in names[:3]])
    se=np.array([data[n]['residence_s']['standard_error'] for n in names[:3]])
    # Three coefficients and three lengths: this is an interpolation with
    # propagated sampling uncertainty, not a goodness-of-fit test.
    A=np.column_stack([x*x,x,np.ones(3)])
    inv=np.linalg.inv(A);coef=inv@times;cov=inv@np.diag(se*se)@inv.T
    a,b,c=coef
    D=L*L/(12*a)
    dD=L*L/(12*a*a)*np.sqrt(cov[0,0])
    ell=2*D*b/L
    grad=np.array([-ell/a,2*D/L,0.])
    dell=np.sqrt(grad@cov@grad)
    closed=data['closed'];fine=data['closed-fine']
    dvalues=np.array([r['finite_time_axial_diffusion_m2_s']['mean'] for r in [closed,fine]])
    dses=np.array([r['finite_time_axial_diffusion_m2_s']['standard_error'] for r in [closed,fine]])
    weights=1/dses**2
    pooled=float(np.sum(weights*dvalues)/np.sum(weights));poolse=float(1/np.sqrt(np.sum(weights)))
    fit={'equation':'mean_time = a*(L/L0)^2 + b*(L/L0) + c',
         'reference_length_m':L,'coefficients_s':dict(zip(['a','b','c'],coef.tolist())),
         'coefficient_covariance_s2':cov.tolist(),
         'quadratic_diffusion_m2_s':{'mean':float(D),'standard_error':float(dD)},
         'effective_end_length_m':{'mean':float(ell),'standard_error':float(dell)},
         'end_fraction_at_reference_length':float(b/(a+b+c)),
         'free_space_pooled_diffusion_m2_s':{'mean':pooled,'standard_error':poolse},
         'qualification':'Three-point interpolation. Errors propagate independent trajectory sampling only. The end length is a descriptive transport parameter, not a designed electrode length.'}
    main=data['main'];local=main['local_state']
    displacement=main['diagnostic_radial_displacement_m']['mean']
    potential=local['radial_E_V_m']*displacement
    frame={'collision_work_eV':main['collision_energy_eV']['mean'],
           'radiated_eV':main['radiated_eV']['mean'],
           'exit_minus_birth_eV':main['exit_or_terminal_energy_eV']['mean']-main['birth_energy_eV'],
           'radial_potential_gain_eV':potential,
           'laboratory_collision_work_eV':main['collision_energy_eV']['mean']+potential,
           'electric_field_work_eV':-potential,
           'qualification':'Frozen local drift frame: laboratory collision work adds v_ExB times azimuthal collision impulse. The electric field takes the same increment into potential energy. Global charge return and voltage conversion remain unspecified.'}
    plant=next(p for p in json.loads((root/'network-plant.json').read_text())['cases']
        if p['assumptions']['B']==140 and p['assumptions']['heat']==.7 and p['assumptions']['end_in']==.999)
    source=plant['selected_source_profiles']['22049']
    transport=json.loads((root/'end-transport-screen.json').read_text())
    allowance=transport['power_allowances']['extra_electrical_MeV_per_primary_at_1mW_kg']
    F=source['reaction_rate_per_m_s']/(2*np.pi*source['radius_m'])
    rows=[]
    for cost in [.01,.1,1.,10.]:
        births=allowance/cost;flux=F*births
        rows.append({'assumed_all_in_MeV_per_electron':cost,
            'maximum_electrons_per_primary':births,
            'fresh_H_atom_supply_m2_s':flux,
            'fresh_H2_molecule_supply_m2_s':flux/2,
            'H2_outgassing_Pa_m3_s_per_m2_at_300K':flux/2*KB*300,
            'unit_flux_equivalent_H2_pressure_Pa_at_300K':flux/2*np.sqrt(2*np.pi*2*MH*KB*300)})
    supply={'primary_events_per_Ti_wall_area_m2_s':float(F),
        'total_Ti_inner_wall_area_m2':float(plant['primary_pn_rate_s']/F),
        'extra_electrical_allowance_MeV_per_primary':allowance,'cases':rows,
        'qualification':'Entire existing margin allocated to electrons from Ti-channel fresh neutral gas. Each H atom supplies at most one new electron without nuclear pair creation; H2 supplies two. All are assumed to enter the hot population. This is a sufficient supply bound under an assigned all-in cost, not a measured desorption yield or actual pressure. Other processors, impurities, electron sources and added hardware also consume the margin. Gross gas-wall recycling flux is not fresh gas supply.'}
    return {'scope':__doc__,'finite_length':fit,'local_energy_frames':frame,'neutral_supply':supply}


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=ROOT)
    ap.add_argument('--outdir',type=Path,default=ROOT);args=ap.parse_args()
    out=analyze(args.root,args.outdir)
    (args.outdir/'lepton-escape-analysis.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps(out,indent=2),flush=True)
