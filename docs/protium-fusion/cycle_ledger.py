"""Steady isotope/energy ledger for separated, repeating staged pulses.

The default data use specified convergent-flow stages, including auxiliaries.
Separation is ideal in the base ledger; optional catalyst collection losses
are charged explicitly. This is a process balance, not a complete plant.
"""
import json
from pathlib import Path
import numpy as np
import pynucastro as p
from nuclear_screen import MU,MEV

ROOT=Path(__file__).parent
STATES=['c12','c13','n14','n15']
DECAY={'n':'p','n13':'c13','o14':'n14','o15':'n15','f17':'o17',
       'f18':'o18','ne18':'o18','ne19':'f19','na21':'ne21','na22':'ne22',
       'mg22':'ne22','mg23':'na23'}
MEAN_NU={'n13':.706,'o14':1.05,'o15':1.00}
FILES={'n14':('compression-A-128-full-burning.json',1),
       'n15':('compression-B1-128-burning.json',1),
       'c12':('compression-B2-128-burning.json',2),
       'c13':('compression-C13-64-burning.json',1)}
FACTORY='compression-He-64-burning.json'

def nuc(name): return p.Nucleus(name)
def add(dst,src,scale=1.):
    for k,v in src.items(): dst[k]=dst.get(k,0.)+scale*v

def sanitize(record):
    raw=record['output_mass_fractions']
    correction=sum(-v for v in raw.values() if v<0)
    if correction>1e-7: raise ValueError('Trace negativity exceeds ledger tolerance')
    norm=sum(max(v,0.) for v in raw.values())
    return {k:max(v,0.)/norm for k,v in raw.items()},correction

def cold(counts):
    final={};q=0.;nu=0.
    for k,v in counts.items():
        dest=DECAY.get(k,k);final[dest]=final.get(dest,0.)+v
        qq=nuc(k).mass-nuc(dest).mass
        q+=v*qq
        # Conservatively discard all decay energy for rare unmodeled spectra.
        nu+=v*MEAN_NU.get(k,max(qq,0.))
    return final,q,nu

def specific(record,key_radial,key_uniform):
    if key_radial in record: return record[key_radial]/record['mass_kg']
    return record[key_uniform]

def stage(cat,file,protons):
    d=json.loads((ROOT/file).read_text());x,neg=sanitize(d)
    mix=d.get('initial_carbon13_number_fraction',0.) if cat=='c12' else 0.
    fuel={cat:1.,'p':protons}
    if mix: fuel={'c12':1-mix,'c13':mix,'p':protons}
    a=sum(nuc(k).A*v for k,v in fuel.items())
    counts={k:v*a/nuc(k).A for k,v in x.items()}
    final,qcold,nucold=cold(counts)
    qhot=sum(nuc(k).mass*v for k,v in fuel.items())-sum(nuc(k).mass*v for k,v in counts.items())
    input_e=specific(d,'initial_energy_J','initial_EOS_energy_J_kg')*a*MU/MEV
    nuhot=specific(d,'weak_energy_escaped_upper_J','weak_energy_escaped_upper_J_kg')*a*MU/MEV
    escaped=specific(d,'radiation_escaped_J','radiated_J_kg')*a*MU/MEV
    kinetic_fraction=d.get('final_kinetic_J',0.)/(d.get('final_kinetic_J',0.)+d.get('final_internal_J',1.))
    return {'source':file,'cat':cat,'protons_input':protons,'fuel_input':fuel,'cold_output':final,
            'Q_hot_MeV':qhot,'Q_cold_MeV':qcold,'Q_total_MeV':qhot+qcold,
            'nu_MeV':nuhot+nucold,'nu_hot_MeV':nuhot,'nu_cold_MeV':nucold,
            'final_kinetic_fraction':kinetic_fraction,
            'radiated_MeV':escaped,'input_MeV':input_e,
            'clipped_negative_mass_fraction':neg}

def factory(file):
    d=json.loads((ROOT/file).read_text());x,neg=sanitize(d)
    a=12/x['c12']
    counts={k:v*a/nuc(k).A for k,v in x.items()}
    final,qcold,nucold=cold(counts)
    input_he=a/4
    qhot=input_he*nuc('he4').mass-sum(nuc(k).mass*v for k,v in counts.items())
    return {'source':file,'input_he4_atoms':input_he,'cold_output':final,
            'Q_total_MeV':qhot+qcold,'Q_hot_MeV':qhot,'Q_cold_MeV':qcold,
            'nu_hot_MeV':specific(d,'weak_energy_escaped_upper_J','weak_energy_escaped_upper_J_kg')*a*MU/MEV,
            'nu_cold_MeV':nucold,
            'final_kinetic_fraction':d.get('final_kinetic_J',0.)/(d.get('final_kinetic_J',0.)+d.get('final_internal_J',1.)),
            'nu_MeV':specific(d,'weak_energy_escaped_upper_J','weak_energy_escaped_upper_J_kg')*a*MU/MEV+nucold,
            'radiated_MeV':specific(d,'radiation_escaped_J','radiated_J_kg')*a*MU/MEV,
            'input_MeV':specific(d,'initial_energy_J','initial_EOS_energy_J_kg')*a*MU/MEV,
            'clipped_negative_mass_fraction':neg}

def ledger(collection_loss=0.,aux_input_multiplier=1.,files=None,factory_file=None,mixed_carbon=False):
    files=files or FILES
    states=['c12','n14','n15'] if mixed_carbon else STATES
    stages=[stage(cat,*files[cat]) for cat in states]
    fac=factory(factory_file or FACTORY)
    def circulating(out,k):
        return out.get(k,0.)+(out.get('c13',0.) if mixed_carbon and k=='c12' else 0.)
    matrix=np.zeros((len(states),len(states)));b=np.array([circulating(fac['cold_output'],k) for k in states])
    for j,s in enumerate(stages):
        for i,k in enumerate(states): matrix[i,j]=circulating(s['cold_output'],k)*(1-collection_loss)
    attempts=np.linalg.solve(np.eye(len(states))-matrix,b)
    if np.any(attempts<0): raise ValueError('Nonphysical stationary solution')
    net={};q=fac['Q_total_MeV'];nu=fac['nu_MeV'];rad=fac['radiated_MeV']
    u=fac['input_MeV']*aux_input_multiplier
    add(net,fac['cold_output']);add(net,{'he4':-fac['input_he4_atoms']})
    for j,(s,f) in enumerate(zip(stages,attempts)):
        add(net,s['cold_output'],f);add(net,s['fuel_input'],-f)
        q+=f*s['Q_total_MeV'];nu+=f*s['nu_MeV'];rad+=f*s['radiated_MeV']
        u+=f*s['input_MeV']*(aux_input_multiplier if s['cat']=='c13' else 1.)
    # Residual circulating catalysts are the explicitly lost output inventory.
    # They enter ash rather than being silently removed from the mass balance.
    helium=net['he4'];hydrogen=-net['p']
    if helium<=0 or hydrogen<=0: raise ValueError('No net hydrogen-to-helium operation')
    baryon=sum(nuc(k).A*v for k,v in net.items())
    qcheck=-sum(nuc(k).mass*v for k,v in net.items())
    assert abs(baryon)/hydrogen<2e-10,baryon
    assert abs(q-qcheck)/q<1e-7,(q,qcheck)
    heat=(q-nu-rad)/helium;input_e=u/helium
    # Use exact mass-defect energy for totals and the final modeled K/U ratio
    # only for the partition. Additional expansion could transfer more U to K.
    kinetic=0.;thermal=0.;cold_heat=0.
    for s,f in [(fac,1.)]+list(zip(stages,attempts)):
        e=s['input_MeV']+s['Q_hot_MeV']-s['nu_hot_MeV']-s['radiated_MeV']
        kinetic+=f*e*s['final_kinetic_fraction']
        thermal+=f*e*(1-s['final_kinetic_fraction'])
        cold_heat+=f*(s['Q_cold_MeV']-s['nu_cold_MeV'])
    # The multiplier penalizes preparation; its extra energy is conservatively
    # assumed to become heat rather than usable directed motion.
    thermal+=u-(fac['input_MeV']+sum(f*s['input_MeV'] for s,f in zip(stages,attempts)))
    kinetic/=helium;thermal/=helium;cold_heat/=helium
    assert abs(kinetic+thermal+cold_heat-(heat+input_e))<1e-7
    ash={k:v/helium for k,v in net.items() if k not in ['p','he4'] and v>1e-8}
    efficiency=[]
    for eta_heat,eta_driver in [(.4,.15),(.4,.2),(.4,.3),(.7,.3),(.85,.5),(.9,.5)]:
        w=eta_heat*(heat+input_e)-input_e/eta_driver
        efficiency.append({'recovery_efficiency':eta_heat,'driver_efficiency':eta_driver,
          'net_MeV_per_helium_before_other_auxiliaries':w,
          'net_fraction_hydrogen_rest_energy':w/(hydrogen/helium*nuc('p').mass)})
    component_cases=[]
    for eta_k,eta_th,eta_driver in [(.8,.4,.3),(.9,.5,.5),(.9,.5,.8),(.95,.6,.8)]:
        w=eta_k*kinetic+eta_th*(thermal+cold_heat)-input_e/eta_driver
        component_cases.append({'kinetic_recovery_efficiency':eta_k,'thermal_recovery_efficiency':eta_th,
          'driver_efficiency':eta_driver,'net_MeV_per_helium_before_other_auxiliaries':w,
          'net_fraction_hydrogen_rest_energy':w/(hydrogen/helium*nuc('p').mass)})
    carbon13_mismatch=net.get('c13',0.)/attempts[0] if mixed_carbon else 0.
    return {'collection_loss_per_stage':collection_loss,'auxiliary_input_multiplier':aux_input_multiplier,
       'mixed_carbon_stage':mixed_carbon,'carbon13_stationarity_residual_per_carbon_stage':carbon13_mismatch,
       'stage_attempts_per_net_helium':dict(zip(states,(attempts/helium).tolist())),
       'new_carbon_per_net_helium':1/helium,'hydrogen_atoms_per_net_helium':hydrogen/helium,
       'ash_nuclei_per_net_helium':ash,
       'ash_baryon_fraction_of_net_fuel':sum(nuc(k).A*v for k,v in ash.items())/(hydrogen/helium),
       'nuclear_Q_MeV_per_helium':q/helium,'neutrino_MeV_per_helium':nu/helium,
       'escaped_radiation_MeV_per_helium':rad/helium,'available_heat_MeV_per_helium':heat,
       'fuel_preparation_MeV_per_helium':input_e,
       'factory_input_MeV_per_helium':fac['input_MeV']*aux_input_multiplier/helium,
       'final_bulk_kinetic_MeV_per_helium':kinetic,
       'final_hot_internal_MeV_per_helium':thermal,
       'cold_decay_heat_MeV_per_helium':cold_heat,
       'baryon_residual_per_hydrogen':baryon/hydrogen,
       'mass_defect_energy_relative_residual':(q-qcheck)/q,
       'net_energy_cases':efficiency,'component_recovery_cases':component_cases,
       'stage_data':stages,'factory_data':fac}

if __name__=='__main__':
    cases=[ledger(loss,mult) for loss,mult in [(0,1),(0,2),(0,10),(1e-5,1),(1e-4,1),(1e-3,1)]]
    out={'scope':'Steady material/energy screen; ideal separation, specified convergent inflow in all stages, no complete plant model.',
         'beta_neutrino_approximation_MeV':MEAN_NU,
         'cases':cases}
    (ROOT/'cycle-ledger.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    for d in cases:
        print({k:d[k] for k in ['collection_loss_per_stage','auxiliary_input_multiplier',
          'hydrogen_atoms_per_net_helium','new_carbon_per_net_helium',
          'available_heat_MeV_per_helium','fuel_preparation_MeV_per_helium']})
        print('Thermal-40% / driver-30% net:',d['net_energy_cases'][2]['net_MeV_per_helium_before_other_auxiliaries'])
