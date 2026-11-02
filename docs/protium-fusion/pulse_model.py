"""Homogeneous expanding-sphere screen for staged catalytic fusion.

This is a deliberately reduced hydrodynamic model, not an implosion design.
Initial conditions assume an already compressed, uniformly hot sphere at rest.
Expansion is homologous; a volume-average pressure accelerates its surface.
The scalar virial equation and PdV work conserve energy within this model.
Radiation loss is an approximate diffusion flux, not a solved boundary layer.
No external confinement, driver, shell, mixing or energy recovery is modeled.

Public arguments and results use SI; reaction network and EOS use cgs.
"""
import argparse
import json
import time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
import pynucastro as p
from nuclear_screen import LIB, MU, MEV, MJ, KB, C, ARAD

SPECIES=['n','p','he4','c12','c13','n13','n14','n15','o14','o15','o16',
         'o17','o18','f17','f18','f19','ne18','ne19','ne20','ne21','ne22',
         'na21','na22','na23','mg22','mg23','mg24']
NUCS=[p.Nucleus(n) for n in SPECIES]
AA=np.array([n.A for n in NUCS],dtype=float)
ZZ=np.array([n.Z for n in NUCS],dtype=float)
MASS=np.array([n.mass for n in NUCS])
RATES=LIB.linking_nuclei(NUCS,with_reverse=True).get_rates()
NR=len(RATES); NN=len(NUCS)
STOICH=np.zeros((NN,NR))
REACT=np.full((NR,3),-1,dtype=np.int64)
COEFF=[]; RATE_INDEX=[]
PREF=np.array([r.prefactor for r in RATES])
DEXP=np.array([r.dens_exp for r in RATES],dtype=np.int64)
for j,r in enumerate(RATES):
    for k,n in enumerate(r.reactants):
        STOICH[NUCS.index(n),j]-=1
        REACT[j,k]=NUCS.index(n)
    for n in r.products: STOICH[NUCS.index(n),j]+=1
    for s in r.sets:
        COEFF.append(s.a); RATE_INDEX.append(j)
COEFF=np.array(COEFF); RATE_INDEX=np.array(RATE_INDEX,dtype=np.int64)
Q=-MASS@STOICH
# Conservative hot-stage convention: all beta Q escapes as neutrinos.
# Actual cold-decay heat must be evaluated separately in the cycle ledger.
QNU=np.array([max(Q[j],0) if r.weak_type else 0 for j,r in enumerate(RATES)])
EOS=p.StellarEOS()


def flows(temp,rho,y):
    t=temp/1e9
    tf=np.array([1.,1/t,t**(-1/3),t**(1/3),t,t**(5/3),np.log(t)])
    vals=np.bincount(RATE_INDEX,weights=np.exp(COEFF@tf),minlength=NR)
    ypad=np.concatenate([np.maximum(y,0),[1.]])
    return vals*PREF*rho**DEXP*np.prod(ypad[REACT],axis=1)


def composition(y):
    comp=p.Composition(NUCS)
    for i,n in enumerate(NUCS): comp[n]=max(y[i]*AA[i],0.)
    comp.normalize()
    return comp


def simulate(stage,temp0,rho0_si,radius0,*,kappa=.001,rtol=2e-7,
             max_expansion=10.,max_steps=300,gravity=True):
    y0=np.zeros(NN)
    if stage=='A':
        y0[SPECIES.index('p')]=1/15
        y0[SPECIES.index('n14')]=1/15
        target='o15'; cats0=1/15
    elif stage=='B':
        y0[SPECIES.index('p')]=3/18
        y0[SPECIES.index('n15')]=1/18
        target='o14'; cats0=1/18
    elif stage=='B1':
        y0[SPECIES.index('p')]=1/16
        y0[SPECIES.index('n15')]=1/16
        target='c12';cats0=1/16
    elif stage=='B2':
        y0[SPECIES.index('p')]=2/14
        y0[SPECIES.index('c12')]=1/14
        target='o14';cats0=1/14
    elif stage=='C0':
        y0[SPECIES.index('p')]=2/18
        y0[SPECIES.index('o16')]=1/18
        target='ne18';cats0=1/18
    elif stage=='C1':
        y0[SPECIES.index('p')]=1/19
        y0[SPECIES.index('f18')]=1/19
        target='o15';cats0=1/19
    elif stage=='C2':
        y0[SPECIES.index('p')]=1/18
        y0[SPECIES.index('o17')]=1/18
        target='n14';cats0=1/18
    elif stage=='C3':
        y0[SPECIES.index('p')]=1/19
        y0[SPECIES.index('o18')]=1/19
        target='n15';cats0=1/19
    elif stage=='C13':
        y0[SPECIES.index('p')]=1/14
        y0[SPECIES.index('c13')]=1/14
        target='n14';cats0=1/14
    elif stage=='D3':
        y0[SPECIES.index('p')]=2/20
        y0[SPECIES.index('o18')]=1/20
        target='c12';cats0=1/20
    elif stage=='He':
        y0[SPECIES.index('he4')]=1/4
        # Normalize carbon output to the maximum possible 3-alpha yield.
        target='c12';cats0=1/12
    else: raise ValueError(stage)
    rho0=rho0_si/1000
    radius0_cm=radius0*100
    mass=4*np.pi/3*rho0_si*radius0**3
    s0=EOS.pe_state(rho0,temp0,composition(y0))
    u0=s0.e
    cs0=np.sqrt(s0.gamma1*s0.p/rho0)
    time_unit=radius0_cm/cs0
    grav=6.67430e-8*mass*1000/radius0_cm if gravity else 0
    # State: Y_i, T/T0, R/R0, V_surface/cs0, E_rad/u0, E_nu/u0.
    z0=np.concatenate([y0,[1.,1.,0.,0.,0.]])
    def rhs(t,z):
        y=z[:NN]
        temp=max(z[NN]*temp0,1e6)
        radius_ratio=max(z[NN+1],.1)
        speed=z[NN+2]*cs0
        rho=rho0/radius_ratio**3
        state=EOS.pe_state(rho,temp,composition(y))
        flow=flows(temp,rho,y)
        dy=STOICH@flow
        qtot=float(Q@flow)*MEV/MU*1e4
        qnu=float(QNU@flow)*MEV/MU*1e4
        r_cm=radius0_cm*radius_ratio
        # Flux limited approximate LTE diffusion; tau uses fixed opacity.
        r_si=r_cm/100; rho_si=rho*1000
        tau=kappa*rho_si*r_si
        flux=C*ARAD*(temp**4)/4/(1+3*tau/4)
        qrad=3*flux/(rho_si*r_si)*1e4
        rhodot=-3*rho*speed/r_cm
        yi=max(float(np.sum(y)),1e-30)
        ye=max(float(ZZ@y),1e-30)
        # Composition derivatives follow exactly from ideal ions plus an
        # electron EOS depending on rho*Ye, and composition-independent rad.
        de_dyi=1.5*KB*1e7/(MU*1000)*temp
        uion=de_dyi*yi
        de_dye=(state.e+rho*state.de_drho-uion)/ye
        tdot=(qtot-qnu-qrad+(state.p/rho**2-state.de_drho)*rhodot
              -de_dyi*np.sum(dy)-de_dye*(ZZ@dy))/state.c_v
        accel=5*state.p/(rho*r_cm)-grav/(radius0_cm*radius_ratio**2)
        result=np.empty_like(z)
        result[:NN]=dy*time_unit
        result[NN]=tdot*time_unit/temp0
        result[NN+1]=speed*time_unit/radius0_cm
        result[NN+2]=accel*time_unit/cs0
        result[NN+3]=qrad*time_unit/u0
        result[NN+4]=qnu*time_unit/u0
        return result
    def expanded(t,z): return z[NN+1]-max_expansion
    expanded.terminal=True; expanded.direction=1
    def cold(t,z): return z[NN]*temp0-3e7
    cold.terminal=True; cold.direction=-1
    def too_hot(t,z): return z[NN]*temp0-2e9
    too_hot.terminal=True; too_hot.direction=1
    begin=time.monotonic()
    sol=solve_ivp(rhs,(0,30),z0,method='BDF',rtol=rtol,atol=1e-12,
                  max_step=30/max_steps,events=[expanded,cold,too_hot])
    z=sol.y[:,-1]; y=z[:NN]; x=z[NN+1]
    end=EOS.pe_state(rho0/x**3,z[NN]*temp0,composition(y))
    nuclear=float(-(y-y0)@MASS)*MEV/MU*1e4
    kinetic=.3*(z[NN+2]*cs0)**2
    egrav0=-.6*grav; egrav=-.6*grav/x
    residual=(end.e+kinetic+egrav+(z[NN+3]+z[NN+4])*u0-u0-egrav0-nuclear)
    ymax=np.max(sol.y[:NN],axis=1)
    heavy=(ZZ>8)
    return {
      'stage':stage,'T0_K':temp0,'rho0_kg_m3':rho0_si,'R0_m':radius0,
      'mass_kg':mass,'mass_over_Jupiter':mass/MJ,
      'initial_EOS_energy_J_kg':u0*1e-4,
      'initial_sound_crossing_s':time_unit,'elapsed_s':float(sol.t[-1]*time_unit),
      'initial_electron_eta':float(s0.eta) if np.isfinite(s0.eta) else None,
      'peak_T_K':float(np.max(sol.y[NN])*temp0),
      'final_T_K':float(z[NN]*temp0),'expansion_factor':float(x),
      'initial_optical_depth':kappa*rho0_si*radius0,
      'final_optical_depth':kappa*rho0_si*radius0/x**2,
      'nuclear_total_J_kg':nuclear*1e-4,'radiated_J_kg':float(z[NN+3]*u0*1e-4),
      'weak_energy_escaped_upper_J_kg':float(z[NN+4]*u0*1e-4),
      'final_thermal_J_kg':float(end.e*1e-4),'final_kinetic_J_kg':float(kinetic*1e-4),
      'deposited_nuclear_over_initial_EOS_energy':float(nuclear/u0-z[NN+4]),
      'energy_residual_fraction_initial':float(residual/u0),
      'mass_fraction_residual':float(AA@y-1),
      'desired_product_per_initial_catalyst':float(y[SPECIES.index(target)]/cats0),
      'final_heavier_than_oxygen_number_per_initial_catalyst':float(y[heavy].sum()/cats0),
      'output_mass_fractions':dict(zip(SPECIES,(AA*y).tolist())),
      'peak_mass_fractions':dict(zip(SPECIES,(AA*ymax).tolist())),
      'integrator_success':sol.success,'integrator_message':sol.message,
      'termination_reason':('expanded' if len(sol.t_events[0]) else
                            'cooled' if len(sol.t_events[1]) else
                            'temperature_guard' if len(sol.t_events[2]) else 'time_limit'),
      'expansion_completed':bool(len(sol.t_events[0]) or len(sol.t_events[1])),
      'rhs_evaluations':sol.nfev,'jacobian_evaluations':sol.njev,
      'wall_seconds':time.monotonic()-begin,'network_nuclei':NN,'network_rates':NR,
      'settings':{'fixed_opacity_m2_kg':kappa,'rtol':rtol,'self_gravity':gravity},
      'limitations':['already-compressed initial state; no driver or implosion',
        'homologous one-zone expansion; no spatial burn or rarefaction front',
        'no screening or Coulomb EOS correction', 'network truncated at magnesium',
        'constant opacity and approximate diffusion; no surface boundary layer',
        'no thermal conduction, mixing, instability or external confinement',
        'all hot beta Q treated as escaping, conservatively',
        'no cold-decay, separation, repeat-cycle, recovery or plant-mass budget']}


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--stage',choices=['A','B','B1','B2','C0','C1','C2','C3','C13','D3','He'],default='A')
    parser.add_argument('--temp',type=float,default=3e8)
    parser.add_argument('--rho',type=float,default=1e8)
    parser.add_argument('--radius',type=float,default=1000)
    parser.add_argument('--rtol',type=float,default=2e-7)
    parser.add_argument('--fast-eos',action='store_true')
    args=parser.parse_args()
    if args.fast_eos:
        from fast_eos import FastEOS
        EOS=FastEOS()
    record=simulate(args.stage,args.temp,args.rho,args.radius,rtol=args.rtol)
    record['interpolated_EOS']=args.fast_eos
    name=f'pulse-{args.stage}-{args.temp:g}-{args.rho:g}-{args.radius:g}'+('-fast' if args.fast_eos else '')+'.json'
    Path(__file__).with_name(name).write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in record.items() if k not in
                     ['output_mass_fractions','peak_mass_fractions','limitations']},indent=2))
