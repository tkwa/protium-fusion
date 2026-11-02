"""Optimistic fixed-state screens for two separated hot-CNO stages.

A: N14 + p -> O15, then allow O15 -> N15 outside the hot state.
B: N15 + 3p -> He4 + O14, then allow O14 -> N14 outside the hot state.

We integrate a finite-fuel reaction network at prescribed T and density.
EOS includes ideal ions, radiation and arbitrary-degeneracy/relativity e+-.
This script does NOT simulate an implosion, ignition or energy extraction.
The sound-speed length is a dimensional diagnostic, not a target design.
"""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
import pynucastro as p
from nuclear_screen import LIB, MU, MEV, MJ

SPECIES=['p','he4','c12','n13','c13','n14','o14','o15','n15','o16']
NUCS=[p.Nucleus(x) for x in SPECIES]
A=np.array([n.A for n in NUCS],dtype=float)
RATE_NAMES=['c12(p,g)n13','n13(p,g)o14','c13(p,g)n14',
            'n14(p,g)o15','n15(p,a)c12','n15(p,g)o16',
            'n13(,)c13','o14(,)n14','o15(,)n15']
RATES=[LIB.get_rate_by_name(x) for x in RATE_NAMES]
EOS=p.StellarEOS()


def stage(stage_name,temp,rho,goal=.5):
    composition=p.Composition(NUCS); composition.set_all(0)
    if stage_name=='A':
        composition['p']=1/15; composition['n14']=14/15
        destination=SPECIES.index('o15'); initial_catalysts=1/15
        time_scale=1/(rho/1000*(1/15)*RATES[3].eval(temp))
    else:
        composition['p']=3/18; composition['n15']=15/18
        destination=SPECIES.index('o14'); initial_catalysts=1/18
        time_scale=1/(rho/1000*(1/18)*min(RATES[0].eval(temp),RATES[1].eval(temp)))
    y0=np.array([composition[n]/n.A for n in NUCS])
    transitions=[]
    for rate in RATES:
        stoich=np.zeros(len(NUCS))
        for n in rate.reactants: stoich[NUCS.index(n)]-=1
        for n in rate.products: stoich[NUCS.index(n)]+=1
        coefficient=rate.eval(temp)*(rho/1000)**(len(rate.reactants)-1)
        transitions.append((stoich,[NUCS.index(n) for n in rate.reactants],coefficient))
    def rhs(t,y):
        dy=np.zeros_like(y)
        for stoich,indices,coefficient in transitions:
            flow=coefficient*np.prod(np.maximum(y[indices],0))
            dy+=stoich*flow
        return dy
    def reached(t,y):
        return y[destination]/initial_catalysts-goal
    reached.terminal=True
    reached.direction=1
    solution=solve_ivp(rhs,(0,20*time_scale),y0,method='BDF',rtol=2e-8,atol=1e-13,events=reached)
    y=solution.y[:,-1]
    elapsed=float(solution.t[-1])
    eos=EOS.pe_state(rho/1000,temp,composition)
    # Atomic masses include electrons; sum accounts for positron annihilation.
    # Beta decay during the hot stage should be negligible; retain a flag.
    masses=np.array([n.mass for n in NUCS])
    nuclear_total=float(-(y-y0)@masses*MEV/MU)
    cs=float(np.sqrt(eos.gamma1*eos.p/(rho/1000))*.01)
    length=cs*elapsed
    mass=4*np.pi/3*rho*length**3
    return {'stage':stage_name,'T_K':temp,'rho_kg_m3':rho,
            'destination_fraction_target':goal,'destination_fraction':float(y[destination]/initial_catalysts),
            'time_at_fixed_state_s':elapsed,'event_reached':bool(len(solution.t_events[0])),
            'initial_EOS_energy_J_kg':float(eos.e*1e-4),
            'initial_EOS_pressure_Pa':float(eos.p*.1),
            'initial_sound_speed_m_s':cs,
            'hot_stage_nuclear_total_J_kg':nuclear_total,
            'nuclear_over_initial_EOS_energy':float(nuclear_total/(eos.e*1e-4)),
            'sound_speed_times_reaction_time_m':length,
            'sphere_mass_at_that_length_kg':float(mass),
            'sphere_mass_at_that_length_over_Jupiter':float(mass/MJ),
            'mass_fraction_residual':float(A@y-1),
            'output_mass_fractions':dict(zip(SPECIES,(A*y).tolist())),
            'assumptions':['prescribed constant T and rho', 'finite initial fuel',
                           'no screening', 'limited CN network plus N15(p,g)O16 leak',
                           'no radiation loss, hydrodynamic expansion or driver model',
                           'EOS has no Coulomb corrections',
                           'hot nuclear energy includes any neutrinos; beta activity negligible only for short runs',
                           'cold decay and recycling energy not included',
                           'sound-speed length is not proof of sufficient confinement']}


if __name__=='__main__':
    records=[stage(s,t,rho) for rho in [1e6,1e7,1e8]
             for t in [2e8,3e8,5e8,8e8] for s in ['A','B']]
    Path(__file__).with_suffix('.json').write_text(json.dumps({'records':records},indent=2)+'\n')
    print('stage | rho SI | T MK | t s | nuclear/u | sound length m | mass/Mj')
    for r in records:
        print(f"{r['stage']} {r['rho_kg_m3']:8.1g} {r['T_K']/1e6:5.0f} {r['time_at_fixed_state_s']:10.3g} "
              f"{r['nuclear_over_initial_EOS_energy']:9.3g} {r['sound_speed_times_reaction_time_m']:12.3g} "
              f"{r['sphere_mass_at_that_length_over_Jupiter']:10.3g}")
