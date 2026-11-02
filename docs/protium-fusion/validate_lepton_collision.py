"""Analytic and numerical checks of the diagnostic Landau energy closure."""
from pathlib import Path
from math import pi,sqrt
import argparse,json,gzip
import numpy as np
from scipy.special import gammainc
from lepton_collision_screen import (ME,MP,EV,E_CHARGE,EPS0,Cylinder,
    ion_coefficients,coefficients,stationary,case,green_kubo_diffusion)

ROOT=Path(__file__).parent


def maxwell(E,T):
    return 2/sqrt(pi)*np.sqrt(E)/T**1.5*np.exp(-E/T)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=ROOT)
    ap.add_argument('--outdir',type=Path,default=ROOT);args=ap.parse_args()
    checks=[]
    def check(name,passed,metric=None):
        checks.append({'name':name,'passed':bool(passed),'metric':metric})
    # At zero speed rotational symmetry requires equal diffusion components.
    for mass,z in [(ME,1),(MP,1),(49*MP,22)]:
        nu,dp,dt=ion_coefficients(np.array([0.,.01,1.]),mass,1e19,z,1000.)
        check(f'Zero-speed isotropic diffusion m/me={mass/ME:g}',
              np.isfinite(nu).all() and dp[0]==dt[0] and dp[0]>0)
    # An independent Ornstein-Uhlenbeck process has exact spatial diffusivity
    # kT/(m nu), although its energy and direction do not decorrelate separately.
    for T,nu in [(10.,1e3),(1000.,2e5)]:
        E=np.geomspace(T*1e-8,T*40,8193);f=maxwell(E,T)
        cc={'B_eV2_s':4*nu*T*E,'direction_rate_s-1':nu*T/E}
        out=green_kubo_diffusion(E,cc,f);actual=out['isotropized_green_kubo_diffusion_m2_s']
        error=abs(actual/(T*EV/ME/nu)-1)
        check(f'Green-Kubo Ornstein-Uhlenbeck diffusivity T={T:g} nu={nu:g}',error<3e-6,error)
    # Maxwellian detailed balance is an independent equilibrium of the
    # Coulomb operator, for any test/background mass ratio at common T.
    for mass in [ME,MP,49*MP]:
        for T in [100.,1000.]:
            E=np.geomspace(T*1e-8,T*40,16385)
            cc=coefficients(E,[(mass,1e19,1,T)],nangle=8)
            s=stationary(E,cc);exact=maxwell(E,T)
            exact/=np.trapezoid(exact,E)
            l1=float(np.trapezoid(np.abs(s['density_per_eV']-exact),E))
            check(f'Maxwellian detailed balance m/me={mass/ME:g} T={T:g}',
                  l1<3e-6 and abs(s['mean_energy_eV']/(1.5*T)-1)<3e-6,l1)
    # Independent Spitzer temperature-exchange moment for two Maxwellians.
    for mass,z in [(MP,1),(49*MP,22)]:
        Ti=1.6e6;n=1e19;logL=15.
        for Te in [100.,3000.,20000.]:
            E=np.geomspace(Te*1e-8,Te*40,16385);f=maxwell(E,Te)
            cc=coefficients(E,[(mass,n,z,Ti)],log_lambda=logL,nangle=8)
            actual=float(np.trapezoid(f*cc['A_eV_s'],E))
            exchange=8*sqrt(2*pi)/3*n*(z*EV**2/(4*pi*EPS0))**2*logL/(ME*mass)
            exchange/=(Te*EV/ME+Ti*EV/mass)**1.5
            expected=1.5*exchange*(Ti-Te)
            error=abs(actual/expected-1)
            check(f'Maxwellian temperature-exchange moment Z={z} Te={Te:g}',error<3e-6,error)
    # Analytic drag convolution for two drifting Maxwellians. This checks
    # angular quadrature, force normalization and frame/sign conventions.
    for mass,z in [(MP,1),(49*MP,22)]:
        Te=2000.;Ti=1.6e6;n=1e19;logL=15.
        E=np.geomspace(Te*1e-8,Te*40,8193);f=maxwell(E,Te)
        for U in [-6e6,1e6,6e6]:
            cc=coefficients(E,[(mass,n,z,Ti)],U,logL,nangle=128)
            actual=float(np.trapezoid(f*cc['force_N'],E))
            vrel2=2*(Te*EV/ME+Ti*EV/mass)
            G=4*pi*(z*EV**2/(4*pi*EPS0))**2*n*logL/ME**2
            expected=ME*G*(1+ME/mass)*U/abs(U)**3*gammainc(1.5,U*U/vrel2)
            error=abs(actual/expected-1)
            check(f'Maxwellian drag moment Z={z} U={U:g}',error<3e-6,error)
    saved=json.loads(gzip.decompress((args.outdir/'lepton-collision-screen.json.gz').read_bytes()))
    E=np.array(saved['energy_eV'])
    for row in saved['cases']:
        f=np.array(row['density_per_eV']);rel=abs(row['collision_eV_s']/row['radiation_eV_s']-1)
        name=f"f={row['reaction_fraction_inside']} logL={row['log_lambda']} flow={row['include_relative_flow']}"
        check('Stationary integrated energy balance '+name,rel<2e-5,rel)
        check('Probability normalization and nonnegative density '+name,
              np.min(f)>=0 and abs(np.trapezoid(f,E)-1)<1e-12)
    plant=next(p for p in json.loads((args.root/'network-plant.json').read_text())['cases']
        if p['assumptions']['B']==140 and p['assumptions']['heat']==.7 and p['assumptions']['end_in']==.999)
    c=Cylinder(plant['selected_source_profiles']['22049'])
    refined=[]
    for row in saved['cases']:
        if row['reaction_fraction_inside'] not in [.1,.5,.99] or row['log_lambda']!=15 or not row['include_relative_flow']:continue
        fine=case(c,row['local_state']['radius_m'],np.geomspace(.001,5e4,8193),15.,nangle=128)
        errors={k:abs(fine[k]/row[k]-1) for k in ['mean_energy_eV','force_N','fixed_speed_spatial_diffusion_proxy_m2_s','isotropized_green_kubo_diffusion_m2_s']}
        check(f"Energy-grid/angular refinement fraction={row['reaction_fraction_inside']}",max(errors.values())<2e-4,errors)
        refined.append({'reaction_fraction_inside':row['reaction_fraction_inside'],'relative_changes':errors})
    out={'scope':__doc__,'checks':checks,'passed':sum(c['passed'] for c in checks),'total':len(checks),'refinement':refined,
        'qualification':'Checks establish the stated unmagnetized isotropic mathematical closure only; strongly magnetized source kinetics remain unresolved.'}
    (args.outdir/'lepton-collision-validation.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(out['passed'],'/',out['total'],'checks pass')
    for c in checks:
        if not c['passed']:print('FAILED',c)
    if out['passed']!=out['total']:raise SystemExit(1)


if __name__=='__main__':main()
