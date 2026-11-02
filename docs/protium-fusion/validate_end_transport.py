"""Independent electrostatic and accounting checks of the end-transport screen.

The nonlinear planar problem checks transmitted positive-current solutions
between prescribed electrodes. It does not validate a finite collector, RF
extraction, transverse optics, or a beam neutralizer.
"""
from pathlib import Path
from math import pi, sqrt, log
import argparse, json
import numpy as np
from scipy.integrate import solve_bvp, quad
from scipy.optimize import brentq
from catalyst_catalog import nucmass

EPS=8.8541878128e-12
E=1.602176634e-19
C=299792458.
ME=9.1093837139e-31
MEV=1.602176634e-13
ROOT=Path(__file__).parent


def relativistic_velocity(za, potential):
    """Compute velocity from momentum, independently of the screen's formula."""
    rest=nucmass(*za)*MEV
    kinetic=za[0]*E*potential
    pc=np.sqrt(kinetic*(kinetic+2*rest))
    return C*pc/(rest+kinetic)


def planar_solution(products, current_density, voltage, gap, n=81, nonrel=False):
    total=sum(p['current_A'] for p in products)
    partial=[(p['ZA'],current_density*p['current_A']/total) for p in products]
    def rhs(x,y):
        remaining=voltage*(1-y[0])
        # The solver may probe an invalid trial during Newton iteration; a
        # reported result must separately satisfy the positive-energy bound.
        remaining=np.maximum(remaining,voltage*1e-3)
        rho=np.zeros_like(remaining)
        for za,j in partial:
            if nonrel:
                v=np.sqrt(2*za[0]*E*remaining/(nucmass(*za)*MEV/C**2))
            else:v=relativistic_velocity(za,remaining)
            rho+=j/v
        return np.vstack([y[1],-gap**2*rho/(EPS*voltage)])
    x=np.linspace(0,1,n)
    sol=solve_bvp(rhs,lambda a,b:np.array([a[0],b[0]]),x,
                  np.array([.36*x*(1-x),.36*(1-2*x)]),tol=2e-9,max_nodes=10000)
    peak=float(sol.sol(.5)[0])
    return sol,peak


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=ROOT)
    ap.add_argument('--screen-dir',type=Path,default=ROOT);args=ap.parse_args()
    screen=json.loads((args.screen_dir/'end-transport-screen.json').read_text())
    energy=json.loads((args.root/'network-ion-budget.json').read_text())
    source=next(p for p in json.loads((args.root/'network-plant.json').read_text())['cases']
                if p['assumptions']['B']==140 and p['assumptions']['heat']==.7 and p['assumptions']['end_in']==.999)
    checks=[];planar=[]
    def check(name,passed,metric=None):checks.append({'name':name,'passed':bool(passed),'metric':metric})
    rate=screen['reference_primary_rate_s'];streams=screen['processor_streams']
    expected=sum(s['pn_per_primary']*(s['ZA'][0]+1)*s['extraction_events_per_pn'] for s in energy['stages'])*rate*E
    actual=sum(s['current_A'] for s in streams)
    check('Network input-pair charge equals exclusive outgoing product current',abs(actual/expected-1)<2e-13,actual)
    check('Grouped and individual product currents agree',abs(sum(p['current_A'] for p in screen['aggregate_products'])/actual-1)<2e-13)
    for s in streams:
        check('Charge ledger '+s['processor'],abs(s['charge_ledger_error_per_pn'])<1e-12,s['charge_ledger_error_per_pn'])
    for case in screen['transport_cases']:
        U=case['energy_per_charge_V']
        # Independently integrate the electric field of a uniform cylinder and
        # its vacuum annulus. The charge is set by a velocity from momentum.
        line=sum(p['current_A']/relativistic_velocity(p['ZA'],U) for p in screen['aggregate_products'])
        for ratio in [1.,.5,.1]:
            a=.02;b=a/ratio
            integrated=quad(lambda r:line*r/(2*pi*EPS*a*a),0,a)[0]
            if b>a:integrated+=quad(lambda r:line/(2*pi*EPS*r),a,b)[0]
            analytic=line/(4*pi*EPS)*(1+2*log(b/a))
            check(f'Cylindrical Gauss-law field integral U={U:g}, a/b={ratio}',abs(integrated/analytic-1)<2e-11)
        check(f'Long merged beam potential U={U:g}',abs(case['long_merged_beam_axis_potential_V']/(line/(4*pi*EPS))-1)<2e-10)
        # At equal velocity, electron kinetic energy is electron/ion mass times
        # the ion kinetic energy. Count Z neutralization electrons per ion.
        k=sum(p['rate_s']/rate*p['ZA'][0]**2*U/1e6*(ME*C*C/(nucmass(*p['ZA'])*MEV))
              for p in screen['aggregate_products'])
        check(f'Comoving electron kinetic ledger U={U:g}',abs(k/case['co_moving_electron_kinetic_MeV_per_primary']-1)<2e-12)
        for detail,s in zip(case['channel_details'],streams):
            d=detail['sufficient_planar_gap_m_for_10percent_spacecharge']
            sol,peak=planar_solution(s['products'],s['aperture_current_density_A_m2'],U,d)
            check(f'Nonlinear mixed-ion planar bound {s["processor"]} U={U:g}',sol.success and 0<peak<.1+2e-8,peak)
            planar.append({'processor':s['processor'],'energy_per_charge_V':U,'gap_m':d,'peak_potential_fraction':peak,'nodes':len(sol.x)})
    # Independent first integral for the single-species, nonrelativistic slab:
    # a*d^2/U^(3/2) = 16/9*(1-t)*(1+2t)^2, t=sqrt(K_min/K_boundary).
    U=10000.;j=20.;beta=.1;m=nucmass(1,1)*MEV/C**2
    gap=sqrt(8*EPS*sqrt(2*E*(1-beta)*U/m)*beta*U/j)
    a=j/EPS*sqrt(m/(2*E));lhs=a*gap**2/U**1.5
    t=brentq(lambda t:16/9*(1-t)*(1+2*t)**2-lhs,.5,1.,xtol=1e-14)
    products=[{'ZA':[1,1],'current_A':j}]
    sol,peak=planar_solution(products,j,U,gap,nonrel=True)
    check('Planar solver matches independent exact first integral',sol.success and abs(peak-(1-t*t))<2e-8,{'bvp':peak,'first_integral':1-t*t})
    sol2,peak2=planar_solution(products,j,U,gap,n=161,nonrel=True)
    check('Planar grid refinement',sol2.success and abs(peak-peak2)<1e-9,abs(peak-peak2))
    # Numerically integrate the positive half of the Maxwell distribution.
    T=energy['T_MeV']*MEV;mp=nucmass(1,1)*MEV/C**2;vth=sqrt(T/mp)
    half_flux=quad(lambda u:u*np.exp(-u*u/2)/sqrt(2*pi),0,12)[0]*vth
    for gate in screen['axial_gate_screen']:
        src=source['selected_source_profiles'][gate['processor']]
        flux=2*src['proton_inventory_per_m']*half_flux
        check('Maxwell flux integral '+gate['processor'],abs(flux/gate['open_H_end_flux_s_per_channel']-1)<1e-12)
        check('Gate duty and closed-tail consistency '+gate['processor'],
              abs(gate['gate_duty_for_one_H_loss_per_processor_pn']*gate['open_H_losses_per_processor_pn']-1)<1e-13 and
              abs(gate['closed_H_losses_at_25kT_per_processor_pn']/gate['open_H_losses_per_processor_pn']/np.exp(-25)-1)<1e-13)
    budget=json.loads((args.screen_dir/'end-transport-budget.json').read_text())['cases']
    check('Additional-load sensitivity is monotonic',all(budget[i+1]['net_W_kg']<budget[i]['net_W_kg'] for i in range(5)))
    check('Independent extra-load limit reaches requested target',abs(budget[4]['net_W_kg']/.001-1)<1e-10,budget[4]['net_W_kg'])
    for row in budget:
        check('Plant accounting '+row['variant'],abs(row['mass_ledger_error_kg'])<2. and abs(row['energy_ledger_error_W'])<1.)
    out={'scope':__doc__,'checks':checks,'passed':sum(c['passed'] for c in checks),'total':len(checks),'planar_solutions':planar}
    (args.screen_dir/'end-transport-validation.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(out['passed'],'/',out['total'],'checks pass; nonlinear potential fractions',min(p['peak_potential_fraction'] for p in planar),max(p['peak_potential_fraction'] for p in planar))
    for c in checks:
        if not c['passed']:print('FAILED',c)
    if out['passed']!=out['total']:raise SystemExit(1)


if __name__=='__main__':main()
