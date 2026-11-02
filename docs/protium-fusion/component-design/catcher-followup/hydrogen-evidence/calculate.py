#!/usr/bin/env python3
"""Conditional hydrogen release kernels; no fitted reactor retention fractions.
Standard library only. Temperatures and histories are prescribed material fields.
Every released quantity uses H nuclei, including H bound in molecular products.
"""
import argparse, hashlib, json, math
from pathlib import Path
KB = 8.617333262145e-5  # eV/K
EV = 1.602176634e-19   # J
MU = 1.66053906660e-27 # kg; matches hydro baryon-mass convention
H_ATOMIC_U = 1.00782503223 # neutral protium atomic mass, for process-ledger comparison
YEAR = 365.25*86400
ROOT = Path(__file__).resolve().parent

def arrhenius(prefactor, energy_eV, temperature_K):
    assert prefactor >= 0 and energy_eV >= 0 and temperature_K > 0
    return prefactor * math.exp(-energy_eV / (KB*temperature_K))

def G_muto(n, q):
    """Positive form of Muto/Morita1992 Eq6, q=n0/C0."""
    assert 0 < n <= 1 and 0 < q <= 1
    return (1/n-1)/q + 2*math.log(n) + q*(1-n)

def G_hasebe(n, q):
    """Hasebe1992 Eq4 in K/SigmaT <<1 limit."""
    assert 0 < n <= 1 and 0 < q <= 1
    return (1/n-1)/q + math.log(n)

def inverse_survival(clock, q, equation=G_muto):
    assert clock >= 0
    if clock == 0: return 1.0
    lo, hi = 1e-300, 1.0
    for _ in range(120):
        mid = (lo+hi)/2
        if equation(mid,q) > clock: lo=mid
        else: hi=mid
    return (lo+hi)/2

def erfc_inverse(f):
    assert 0 < f < 1
    lo,hi=0.,32.
    for _ in range(100):
        mid=(lo+hi)/2
        if math.erfc(mid)>f:lo=mid
        else:hi=mid
    return (lo+hi)/2

def diffusion_release(clock_m2, distance_m):
    """1D half-space first passage to an absorbing plane, delta initial depth.
No retrapping, no finite rear surface. D(T) scalar, spatially uniform.
"""
    assert clock_m2 >= 0 and distance_m > 0
    return 0. if clock_m2 == 0 else math.erfc(distance_m/(2*math.sqrt(clock_m2)))

def required_diffusivity(distance_m, time_s, fraction):
    return distance_m**2/(4*erfc_inverse(fraction)**2*time_s)

def integrate_clock(times, temperatures, prefactor, energy_eV):
    """Cumulative Arrhenius clock: linear T between nodes, 256-panel Simpson.
Convergence is checked against1024 panels in validate.py. ConstantT is exact.
"""
    assert len(times)==len(temperatures) and len(times)>=2 and times[0]==0
    out=[0.]
    for t0,t1,T0,T1 in zip(times,times[1:],temperatures,temperatures[1:]):
        assert t1>t0 and min(T0,T1)>0
        n=256;dt=(t1-t0)/n
        vals=[arrhenius(prefactor,energy_eV,T0+(T1-T0)*i/n) for i in range(n+1)]
        area=dt/3*(vals[0]+vals[-1]+4*sum(vals[1:-1:2])+2*sum(vals[2:-1:2]))
        out.append(out[-1]+area)
    return out

def response(config):
    """One initial implanted cohort, independent disjoint populations.
Fractions/rates and gas routing must be supplied. Mixed populations are an
explicit scenario, never inferred from TDS peaks. Muto is nonlinear in occupancy:
repeated pulses require updating q and capacities, not naive convolution.
"""
    times=config['times_s']; temps=config['temperature_K']; N=config['initial_H_atoms_m2']
    pops=config['populations']; assert N>=0
    assert abs(sum(p['fraction'] for p in pops)-1)<1e-12
    total=[0.]*len(times); components={}
    for p in pops:
        assert 0<=p['fraction']<=1
        kind=p['model']; f=p['fraction']
        if kind=='prompt': vals=[1.]*len(times)
        elif kind=='held': vals=[0.]*len(times)
        elif kind=='diffusion':
            clock=integrate_clock(times,temps,p['D0_m2_s'],p['E_eV'])
            vals=[diffusion_release(c,p['distance_m']) for c in clock]
        elif kind=='first_order':
            clock=integrate_clock(times,temps,p['prefactor_s'],p['E_eV'])
            vals=[-math.expm1(-c) for c in clock]
        elif kind=='muto':
            clock=integrate_clock(times,temps,8.2,1.)
            vals=[1-inverse_survival(c,p['q']) for c in clock]
        else: raise ValueError(kind)
        components[p['name']]=[N*f*v for v in vals]
        total=[a+N*f*b for a,b in zip(total,vals)]
    cavity=config['cavity_fraction']; atom=config['atomic_H_fraction']; methane=config['methane_H_fraction']
    assert 0<=cavity<=1 and min(atom,methane)>=0 and atom+methane<=1
    return {'scope':'conditional independent-population response; prescribed T, no reactor-fitted fractions',
            'times_s':times,'released_H_atoms_m2':total,'retained_H_atoms_m2':[N-v for v in total],
            'released_cavity_H_atoms_m2':[cavity*v for v in total],
            'released_collector_H_atoms_m2':[(1-cavity)*v for v in total],
            'atomic_H_particles_m2':[atom*v for v in total],
            'H2_molecules_m2':[(1-atom-methane)*v/2 for v in total],
            'CH4_molecules_and_C_atoms_m2':[methane*v/4 for v in total],
            'population_released_H_atoms_m2':components}

def tables(snapshot):
    rows=[]
    for T in [723.15,773.15,823.15,873.15,1000.,1500.,1668.,1700.]:
        for q in [.001,.01,.1,.95]:
            a=arrhenius(8.2,1.,T)
            rows.append({'T_K':T,'q':q,'a_s_inv':a,'half_release_s':G_muto(.5,q)/a,
                         'ninety_percent_release_s':G_muto(.1,q)/a,
                         'released_after_B1_period':1-inverse_survival(a*snapshot['stages']['B1']['period_s'],q),
                         'status':'fit-range kernel; q is scenario' if T<=873.15 else 'high-temperature extrapolation; mechanism can change'})
    diff=[]
    for T in [1000.,1500.,1700.]:
        D=arrhenius(1.2e-4,.9,T)
        for L in [.1e-6,1e-6,10e-6,40e-6,100e-6]:
            diff.append({'T_K':T,'distance_m':L,'conditional_surface_D_m2_s':D,
                         't50_s':required_diffusivity(L,1,.5)/D,'t90_s':required_diffusivity(L,1,.9)/D,
                         'status':'Warrier2006 surface-channel coefficient; only accessible mobile population'})
    requirements=[{'distance_m':L,'time_s':t,'D_required_50_m2_s':required_diffusivity(L,t,.5),
                   'D_required_90_m2_s':required_diffusivity(L,t,.9)}
                  for L in [.1e-6,1e-6,10e-6,40e-6,100e-6]
                  for t in [100.,snapshot['stages']['B1']['period_s'],86400.]]
    flow=0.; cadence=0.; stages={}
    for s,x in snapshot['stages'].items():
        N=x['H_atoms_m2']; P=x['period_s']; flow+=N/P;cadence+=1/P
        clock=2*2.9e-22*x['He_atoms_m2']
        stages[s]={**x,'H_kg_m2_year_baryon':N*MU/P*YEAR,'H_kg_m2_year_neutral_atomic':N*MU*H_ATOMIC_U/P*YEAR,'H_to_C_if_uniform_10um':N/(2000*10e-6/(12*MU)),
                   'H_to_C_if_uniform_100um':N/(2000*100e-6/(12*MU)),
                   'He_induced_clock_if_1p5MeV':clock,
                   'He_induced_release_if_q_0p001':1-inverse_survival(clock,.001,G_hasebe),
                   'He_induced_release_if_q_0p95':1-inverse_survival(clock,.95,G_hasebe)}
    return {'scope':'bounds and sensitivity kernels, no asserted reactor release fraction',
            'ideal_laminate_paths':[{'T_K':T,'side_distance_m':L,'ideal_Dparallel_m2_s':arrhenius(2e-7,.38,T),'t90_s':required_diffusivity(L,1,.9)/arrhenius(2e-7,.38,T),'required_D90_over_B1_period_m2_s':required_diffusivity(L,snapshot['stages']['B1']['period_s'],.9),'scope':'ideal-crystal primary MD coefficient; effective damaged tensor unspecified'} for T in [1000.,1500.,1700.] for L in [10e-6,50e-6]],
            'muto_sensitivity':rows,'mobile_surface_channel':diff,'diffusivity_requirements':requirements,
            'stages':stages,'mean_H_atoms_m2_s':flow,'mean_H_kg_m2_year_baryon':flow*MU*YEAR,'mean_H_kg_m2_year_neutral_atomic':flow*MU*H_ATOMIC_U*YEAR,
            'whole_sphere_H_kg_s_baryon':flow*MU*snapshot['area_m2'],'whole_sphere_H_kg_s_neutral_atomic':flow*MU*H_ATOMIC_U*snapshot['area_m2'],
            'power_W_m2_per_eV_per_H':flow*EV,
            'input_power_W_m2_per_eV_positive_deltaG_per_H_at_eta_0p5':flow*EV/.5,
            'all_stage_pulses_per_year':cadence*YEAR,
            'broad_layer_refurbishment_sensitivity':[{'depth_m':d,'kg_m2_per_event':2000*d,
               'kg_m2_year_every_event':2000*d*cadence*YEAR,
               'W_m2_per_MJ_kg_input':2000*d*cadence*1e6} for d in [1e-6,10e-6,100e-6]],
            'source_input_sha256':snapshot['surface_spectra_sha256']}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--config',type=Path);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    result=response(json.loads(a.config.read_text())) if a.config else tables(json.loads((ROOT/'input-snapshot.json').read_text()))
    a.out.write_text(json.dumps(result,indent=2)+'\n')
