"""Conservation, independent closures, and thermodynamic benchmark checks."""
from pathlib import Path
import argparse,json,math
import numpy as np
from scipy.integrate import quad
from economic_comparison import YEAR,KB,opaque_replacement,discounted_payback
from nonneutral_screen import mass,MEV,C


def validate(root,data):
    d=json.loads((data/'economic-comparison.json').read_text());checks=[]
    def check(name,value,limit):checks.append(dict(name=name,metric=float(value),upper_limit=limit,passed=bool(value<limit)))
    base=json.loads((root/'plant-budget.json').read_text())['cases'][0]
    control=d['opaque_replacement_cases'][0]
    check('Zero-replacement closure recovers original opaque mass',abs(control['total_mass_kg']/base['total_mass_kg']-1),2e-12)
    check('Zero-replacement closure recovers original opaque throughput',abs(control['helium_atoms_s']/base['operations']['helium_atoms_per_s']-1),2e-12)
    check('Zero-replacement closure recovers original opaque support',abs(control['support_stored_energy_J']/base['operations']['support_stream_kinetic_energy_J']-1),4e-12)
    for row in d['opaque_replacement_cases']:
        label=f'replacement={row["replacement_fraction_year"]}'
        check('Mass closure '+label,abs(row['mass_relative_residual']),2e-12)
        check('Independent quadratic closure '+label,abs(row['independent_quadratic_mass_relative_error']),2e-12)
        check('Electric conservation '+label,abs(row['electric_balance_W'])/row['gross_electric_W'],2e-14)
        check('Heat conservation '+label,abs(row['thermal_balance_W'])/row['nuclear_captured_W'],2e-14)
        check('Stable lower mass fixed point '+label,row['mass_fixed_point_derivative'],1.)
    # Sum the measured-mass endpoints independently of the reaction-path Qs.
    direct=(6*mass('d')-2*mass('he4')-2*mass('p')-2*mass('n'))
    check('Six-deuteron channel sum against endpoint mass defect',abs(direct-d['DD_nuclear_ledger']['six_D_secondary_burn_Q_MeV']),2e-9)
    check('Neutron recapture closes D-to-helium mass defect',abs(d['DD_nuclear_ledger']['closure_residual_MeV']),2e-9)
    # Integrate the reversible chemical-potential cost as D is extracted from
    # an ideal binary atomic mixture; this independently recovers mixing entropy.
    for ratio in [1.9e-5,2.6e-5,3.3e-5,1e-7]:
        row=next(x for x in d['deuterium_scenarios'] if x['D_over_H_number_ratio']==ratio)
        x=row['D_atom_fraction'];NH=1-x
        value=quad(lambda t:2*x*t*math.log((NH+x*t*t)/(x*t*t)),0,1,epsabs=1e-17,epsrel=2e-11)[0]
        mbar=NH*d['constants']['hydrogen_atomic_mass_kg']+x*d['constants']['deuterium_atomic_mass_kg']
        work=value*KB*300/mbar
        check(f'Integrated chemical-potential separation work D/H={ratio}',abs(work/row['ideal_atomic_mixture_separation_J_per_kg_feed']-1),2e-9)
    # Discounted payback is checked by quadrature of the delivered energy.
    for row in d['capital_scenarios']:
        if row['construction_J_kg']!=1e8:continue
        for case in row['discounted_payback']:
            r=case['rate_year'];t=case['payback_year'];s=d['nonneutral_reference']['net_W_kg']
            if t is None:
                check(f'No finite discounted payback at r={r}, field={row["field_charge_included"]}',row['maximum_continuous_energy_growth_rate_year']/r,1.+1e-13)
            else:
                work=quad(lambda years:s*YEAR*math.exp(-r*years),0,t,epsrel=1e-12)[0]
                check(f'Discounted output repays capital at r={r}, field={row["field_charge_included"]}',abs(work/row['total_capital_energy_J_kg']-1),2e-12)
    out={'passed':sum(x['passed'] for x in checks),'total':len(checks),'checks':checks,
        'scope':'Arithmetic and physical bookkeeping checks; no validation of future component performance or market assumptions.'}
    (data/'economic-comparison-validation.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:out[k] for k in ['passed','total']}))
    for x in checks:
        if not x['passed']:print('FAILED',x)
    if out['passed']!=out['total']:raise SystemExit(1)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=Path(__file__).parent);ap.add_argument('--data',type=Path,default=Path(__file__).parent);a=ap.parse_args();validate(a.root,a.data)
