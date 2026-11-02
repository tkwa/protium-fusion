"""Check the gyro-averaged stochastic integrator against independent limits."""
from pathlib import Path
from math import sqrt,pi,exp
import argparse,json,gzip,time
import numpy as np
from scipy.special import gammainc
from scipy.integrate import quad
from numba import njit
from lepton_escape_mc import (ME,MP,EV,prepare_species,collision_coeff,
    rotate_direction,ensemble)
from lepton_collision_screen import ion_coefficients

ROOT=Path(__file__).parent


@njit(cache=True)
def direction_samples(q,n,seed):
    np.random.seed(seed)
    rows=np.empty((n,4))
    for i in range(n):
        x,y,z=rotate_direction(0.,0.,1.,q,.02)
        rows[i,0]=x;rows[i,1]=y;rows[i,2]=z
        rows[i,3]=x*x+y*y+z*z
    return rows


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=ROOT)
    ap.add_argument('--outdir',type=Path,default=ROOT);ap.add_argument('--without-cases',action='store_true')
    args=ap.parse_args();checks=[];details={};start=time.monotonic()
    def check(name,passed,metric=None):checks.append({'name':name,'passed':bool(passed),'metric':metric})
    def sample_check(name,samples,expected,extra=0.):
        mean=float(np.mean(samples));se=float(np.std(samples,ddof=1)/sqrt(len(samples)))
        check(name,abs(mean-expected)<5*se+extra,{'mean':mean,'expected':expected,'standard_error':se,'extra_discretization_tolerance':extra})
    speeds=np.r_[0.,np.geomspace(1.,1e9,301)]
    for mass,z in [(ME,1),(MP,1),(49*MP,22)]:
        table=prepare_species([(mass,1e19,z,1.6e6)])
        actual=np.array([collision_coeff(v,table) for v in speeds])
        reference=np.array(ion_coefficients(speeds,mass,1e19,z,1.6e6)).T
        error=float(np.max(np.abs(actual[:,:3]/reference-1)))
        check(f'Scalar/array collision tensor m/me={mass/ME:g}',error<1e-8,error)
        nu,dp,dt=reference.T;identity=-2*nu*speeds**2+dp+2*dt
        scale=2*nu*speeds**2+dp+2*dt
        error=float(np.max(np.abs(actual[:,3]-identity)/scale))
        check(f'Stable energy drift matches Ito moment m/me={mass/ME:g}',error<1e-10,error)
    for q in [.02,.2,2.,21.]:
        rows=direction_samples(q,20000,31119+int(q*100))
        sample_check(f'Spherical first harmonic q={q}',rows[:,2],exp(-q),1e-12)
        sample_check(f'Zero transverse angular drift q={q}',rows[:,0],0.,1e-12)
        check(f'Spherical norm q={q}',np.max(abs(rows[:,3]-1))<3e-14,float(np.max(abs(rows[:,3]-1))))
    # Exact Lorentz process: fixed energy and exponential directional memory.
    birth=1000.;v=sqrt(2*birth*EV/ME);nu=1e5;T=3/nu
    table=np.array([[nu*v**3,0.,1e-30]])
    rows,_,_=ensemble(8192,20173,table,0.,1e300,20235.,T,.02,birth,0.,False,137.)
    error=float(np.max(abs(rows[:,3]/birth-1)))
    check('Lorentz scattering conserves every particle energy',error<1e-10,error)
    sample_check('Lorentz velocity autocorrelation',rows[:,15]*rows[:,16],v*v/3*exp(-nu*T),v*v*1e-7)
    expected=2*v*v/(3*nu)*(T-(1-exp(-nu*T))/nu)
    sample_check('Lorentz exact finite-time mean-square displacement',rows[:,13]**2,expected,expected*2e-4)
    check('Lorentz nonnegative relative-energy update',rows[:,14].sum()==0,int(rows[:,14].sum()))
    # Pure perpendicular radiation has a closed-form trajectory for every
    # initial pitch, including the retained parallel component.
    tau=1.3712655035427434e-4;T=3*tau;birth=1000.
    rows,_,_=ensemble(4096,71117,np.array([[0.,0.,1.]]),0.,tau,20235.,T,.02,birth,0.,False,137.)
    initial_parallel=.5*ME*rows[:,15]**2/EV
    expected=initial_parallel+(birth-initial_parallel)*exp(-T/tau)
    error=float(np.max(abs(rows[:,3]-expected)))
    check('Exact pitch-dependent radiative cooling',error<1e-7,error)
    error=float(np.max(abs(rows[:,13]-rows[:,15]*T)))
    check('Radiation leaves axial free streaming unchanged',error<1e-6,error)
    # A Maxwellian electron bath has a Maxwellian stationary test distribution.
    bathT=100.;table=prepare_species([(ME,1e19,1,bathT)])
    v=sqrt(2*bathT*EV/ME);nu=collision_coeff(v,table)[0]
    rows,_,_=ensemble(8192,18133,table,0.,1e300,20235.,30/nu,.01,1.,0.,False,137.)
    sample_check('Maxwellian bath stationary mean energy',rows[:,3],1.5*bathT,.003*bathT)
    energies=np.sort(rows[:,3]);exact=gammainc(1.5,energies/bathT)
    empirical=(np.arange(len(rows))+.5)/len(rows)
    ks=float(np.max(abs(empirical-exact)))
    check('Maxwellian bath stationary energy CDF',ks<3/sqrt(len(rows))+.003,ks)
    sample_check('Maxwellian bath axial variance',rows[:,16]**2,bathT*EV/ME,.003*bathT*EV/ME)
    check('Maxwellian bath nonnegative energy update',rows[:,14].sum()==0,int(rows[:,14].sum()))
    if not args.without_cases:
        names=['main','fine','coarse','closed','half','double','closed-fine']
        data={name:json.loads((args.outdir/f'lepton-escape-{name}.json').read_text()) for name in names}
        for name,row in data.items():
            check('No relative-energy truncations '+name,row['negative_radial_energy_proposals']==0,row['negative_radial_energy_proposals'])
            check('Pathwise energy ledger '+name,row['pathwise_energy_closure_max_eV']<1e-5,row['pathwise_energy_closure_max_eV'])
            check('Normalized occupation histograms '+name,
                abs(sum(row['occupation_probability_per_energy_bin'])-1)<1e-10 and abs(sum(row['occupation_probability_per_pitch_bin'])-1)<1e-10)
            if row['absorbing_ends']:check('All sampled electrons reach an end '+name,row['censored']==0,row['censored'])
        for a,b in [('coarse','main'),('main','fine')]:
            for metric in ['residence_s','radiated_eV','exit_or_terminal_energy_eV']:
                aa=data[a][metric];bb=data[b][metric]
                se=sqrt(aa['standard_error']**2+bb['standard_error']**2)
                change=abs(aa['mean']-bb['mean'])
                check(f'Step refinement {a}/{b} {metric}',change<4*se+.01*bb['mean'],
                    {'absolute_change':change,'combined_standard_error':se,'relative_change':change/bb['mean']})
        aa=data['closed']['finite_time_axial_diffusion_m2_s']
        bb=data['closed-fine']['finite_time_axial_diffusion_m2_s']
        se=sqrt(aa['standard_error']**2+bb['standard_error']**2)
        check('Free-space step controls agree within sampling precision',abs(aa['mean']-bb['mean'])<4*se,
            {'change_in_combined_standard_errors':abs(aa['mean']-bb['mean'])/se})
        from lepton_escape_analysis import analyze,KB,MH
        analysis=analyze(args.root,args.outdir)
        fitted=analysis['finite_length']['quadratic_diffusion_m2_s']
        pooled=analysis['finite_length']['free_space_pooled_diffusion_m2_s']
        se=sqrt(fitted['standard_error']**2+pooled['standard_error']**2)
        check('Length coefficient agrees with independent free-space control',abs(fitted['mean']-pooled['mean'])<4*se,
            {'change_in_combined_standard_errors':abs(fitted['mean']-pooled['mean'])/se})
        # Independently check the Galilean energy ledger for arbitrary kicks.
        rng=np.random.default_rng(42879);v=rng.normal(size=(1000,3))*3e7
        dv=rng.normal(size=(1000,3))*2e5;flow=np.array([0.,-6.2e6,0.])
        delta=lambda a:.5*ME*np.sum(2*a*dv+dv*dv,axis=1)/EV
        discrepancy=delta(v+flow)-delta(v)-ME*np.sum(dv*flow,axis=1)/EV
        check('Independent laboratory/drift-frame work transform',np.max(abs(discrepancy))<1e-10,float(np.max(abs(discrepancy))))
        # The one-sided Maxwell flux is integrated in standard-normal units.
        moment=quad(lambda x:x*np.exp(-x*x/2)/sqrt(2*pi),0.,np.inf,epsabs=1e-12)[0]
        T=300.;mass=2*MH;pressure=1.
        numerical=pressure/(KB*T)*sqrt(KB*T/mass)*moment
        analytic=pressure/sqrt(2*pi*mass*KB*T)
        check('Neutral supply conversion: one-sided Maxwell flux',abs(numerical/analytic-1)<1e-12,float(numerical/analytic-1))
        supply=analysis['neutral_supply'];F=supply['primary_events_per_Ti_wall_area_m2_s']
        for row in supply['cases']:
            cost=row['assumed_all_in_MeV_per_electron'];births=row['maximum_electrons_per_primary']
            check(f'Parameterized supply energy budget {cost:g} MeV',
                abs(cost*births/supply['extra_electrical_allowance_MeV_per_primary']-1)<1e-13)
            flux=row['unit_flux_equivalent_H2_pressure_Pa_at_300K']*numerical
            check(f'Pressure-to-fresh-flux units {cost:g} MeV',
                abs(2*flux/(F*births)-1)<1e-12)
        details['production_cases']=data
    out={'scope':__doc__,'passed':sum(r['passed'] for r in checks),'total':len(checks),'checks':checks,
        'seconds':time.monotonic()-start,'details':details,
        'qualification':'Verifies this prescribed unmagnetized stochastic benchmark; not the strong-field collision operator or real collector geometry.'}
    name='lepton-escape-validation-preliminary.json' if args.without_cases else 'lepton-escape-validation.json'
    (args.outdir/name).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(out['passed'],'/',out['total'],'checks pass in',out['seconds'],'s',flush=True)
    for row in checks:
        if not row['passed']:print('FAILED',row,flush=True)
    if out['passed']!=out['total']:raise SystemExit(1)


if __name__=='__main__':main()
