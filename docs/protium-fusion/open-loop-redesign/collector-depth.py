"""Monoenergetic finite-depth prompt deposition sensitivity, not transport.

Insulated semi-infinite constant-property carbon with prescribed volume heat
source q(t)*mu*exp(-mu*x). This normalizes ALL prompt energy to that deposition
profile. XCOM total attenuation and energy-absorption coefficients provide two
illustrative depths; neither substitution solves scattering/secondary escape.
"""
from pathlib import Path
import hashlib,json,math
import numpy as np
from scipy.signal import fftconvolve
from scipy.special import erfcx
from scipy.optimize import brentq
HERE=Path(__file__).resolve().parent
K=100.;RHO=2000.;CP=1000.;ALPHA=K/(RHO*CP);RREF=1e9;AREA=4*math.pi*RREF**2
# NIST graphite: energy keV, total mu/rho and mu_en/rho in cm2/g.
TABLE=[(1,2211,2209),(3,90.33,89.63),(5,19.12,18.66),(10,2.373,2.078),
       (30,.2562,.06614),(100,.1514,.02147)]

def integrated_kernel(t,mu):
    z=mu*np.sqrt(ALPHA*np.asarray(t))
    value=erfcx(z)-1+2*z/math.sqrt(math.pi)
    # Stable small-z expansion avoids cancellation between its two O(z) terms.
    small=z<.005;zs=z[small]
    value[small]=zs**2-4*zs**3/(3*math.sqrt(math.pi))+zs**4/2-8*zs**5/(15*math.sqrt(math.pi))+zs**6/6
    return value/(mu*K)

def surface(q,dt):
    n=np.arange(len(q));w=2*np.sqrt(dt)*(np.sqrt(n+1)-np.sqrt(n))/math.sqrt(math.pi*K*RHO*CP)
    return fftconvolve(q,w)[:len(q)]

def volume(q,dt,mu):
    edges=np.arange(len(q)+1)*dt;w=np.diff(integrated_kernel(edges,mu))
    assert np.min(w)>-1e-18
    return fftconvolve(q,w)[:len(q)]

def main():
    inp=json.loads((HERE/'collector-inputs.json').read_text());base=json.loads((HERE/'collector.json').read_text())
    stages={};checks=[]
    for cat,item in inp['stages'].items():
        time=np.asarray(item['prompt_time_s']);energy=np.asarray(item['prompt_cumulative_J'])
        grid=np.linspace(0,time[-1],12001);dt=grid[1]-grid[0]
        q=np.diff(np.interp(grid,time,energy))/(dt*AREA)
        s=float(max(surface(q,dt)));rows=[]
        for kev,mu_total,mu_energy in TABLE:
            row={'photon_energy_keV':kev}
            for name,coefficient in [('total_attenuation_depth',mu_total),('energy_absorption_depth',mu_energy)]:
                mu=coefficient*.1*RHO
                rise=float(max(volume(q,dt,mu)))
                row[name]={'mass_coefficient_m2_kg':coefficient*.1,'deposition_e_folding_depth_m':1/mu,
                    'prompt_peak_rise_K':rise,'ratio_to_surface_peak':rise/s,
                    'adiabatic_no_diffusion_peak_K':sum(q)*dt*mu/(RHO*CP)}
                assert rise<=row[name]['adiabatic_no_diffusion_peak_K']*(1+1e-8)
                assert rise<=s*(1+1e-8)
            rows.append(row)
        stages[cat]={'surface_peak_rise_K':s,'prompt_fluence_J_m2':sum(q)*dt,'monoenergetic_sensitivities':rows}
        checks.append({'test':'Prompt surface kernel agrees baseline '+cat,'error':s/base['stages'][cat]['prompt_only_peak_rise_K']-1})
    # Independent asymptotic checks: no-diffusion finite-depth limit and surface limit.
    q=np.ones(1000)*2e5;dt=.01
    exact=2*q[0]*math.sqrt(10/(math.pi*K*RHO*CP))
    checks.append({'test':'Surface rectangular exact','error':surface(q,dt)[-1]/exact-1})
    large=volume(q,dt,1e9)[-1]
    checks.append({'test':'Large mu approaches surface result','error':large/exact-1})
    small=volume(q,dt,.001)[-1];adiabatic=sum(q)*dt*.001/(RHO*CP)
    checks.append({'test':'Small mu approaches adiabatic depth result','error':small/adiabatic-1})
    assert max(abs(x['error']) for x in checks)<1e-4
    radius=[]
    for kev,_,_ in TABLE:
        for mode in ['total_attenuation_depth','energy_absorption_depth']:
            rvalues={}
            for cat,stage in stages.items():
                depth=next(x[mode]['prompt_peak_rise_K'] for x in stage['monoenergetic_sensitivities'] if x['photon_energy_keV']==kev)
                b=base['stages'][cat]
                rvalues[cat]=brentq(lambda r:b['material_peak_rise_K']*(RREF/r)**2.5+(b['delayed_decay_peak_rise_K']+depth)*(RREF/r)**2-1000,1e6,1e10)
            radius.append({'energy_keV':kev,'deposition_model':mode,'allowed_rise_K':1000.,'minimum_R_m':max(rvalues.values()),'stage_radii_m':rvalues})
    mixed=[]
    for radius_m in [6.5e8,7e8,7.5e8,8e8]:
        # Sufficient bound obtained by adding separate prompt maxima. A soft
        # fraction uses the original all-surface result; hard fraction at10keV
        # uses total-mu deposition, the more concentrated of our two closures.
        limits={}
        for cat,stage in stages.items():
            b=base['stages'][cat];x=RREF/radius_m
            hard=next(z['total_attenuation_depth']['prompt_peak_rise_K'] for z in stage['monoenergetic_sensitivities'] if z['photon_energy_keV']==10)
            remaining=1000-b['material_peak_rise_K']*x**2.5-b['delayed_decay_peak_rise_K']*x**2
            limits[cat]=(remaining/x**2-hard)/(stage['surface_peak_rise_K']-hard)
        mixed.append({'radius_m':radius_m,'allowed_rise_K':1000.,'maximum_soft_surface_component_energy_fraction_raw':min(limits.values()),'per_stage_raw_limits':limits})
    result={'scope':__doc__,'baseline_unchanged':True,'radius_reference_m':RREF,
      'material_properties':{'k_W_m_K':K,'density_kg_m3':RHO,'cp_J_kg_K':CP},
      'stages':stages,'radius_sensitivities':radius,'soft_component_limits_with10keV_hard_component':mixed,
      'checks':checks,'all_checks_passed':True,
      'source':'https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z06.html',
      'limitations':['No escaping photon spectrum is calculated by the grey pulse model; these cases cannot replace the baseline without that spectrum.',
       'Total-mu depth puts all heat at the primary interaction depth; mu_en depth puts all heat at an energy-absorption length. Both prescribe deposition and omit actual multiple scattering, electron transport and escaping energy.',
       'Total-mu is the more concentrated of these two profiles, not a universal upper bound for arbitrary material/spectral response.',
       'Constant properties, flat semi-infinite wall, no erosion or nonlinear radiation-hydrodynamic response.',
       'The bulk kinetic and delayed-decay terms are retained from the provisional64-zone baseline.'],
      'input_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [HERE/'collector.json',HERE/'collector-inputs.json']},
      'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (HERE/'collector-depth.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'C13':stages['c13'],'radii':radius,'soft_limits':mixed,'checks':checks},indent=2))
if __name__=='__main__':main()
