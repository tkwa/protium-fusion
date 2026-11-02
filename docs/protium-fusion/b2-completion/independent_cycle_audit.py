"""Independently normalized closed-cycle balance; raw trace records preferred.

Uses a five-variable linear system normalized to one net He4 rather than
cycle_ledger's factory-normalized four-variable fixed-point solve. Atomic mass
excess provides Q independently of the full atomic mass sums in cycle_ledger.
"""
import argparse,json,sys,hashlib
from pathlib import Path
import numpy as np
import pynucastro as pn
BASE=Path(__file__).resolve().parent.parent
if not (BASE/'fast_eos.py').exists():BASE=Path('/home/tkwa/.codex/worktrees/protium-fusion/website/docs/protium-fusion')
sys.path.insert(0,str(BASE))
MU=1.66053906660e-27;MEV=1.602176634e-13
CATS=['c12','c13','n14','n15']
STAGES={'c12':('B2',2),'c13':('C13',1),'n14':('A',1),'n15':('B1',1)}
DECAY={'n':'p','n13':'c13','o14':'n14','o15':'n15','f17':'o17','f18':'o18','ne18':'o18','ne19':'f19','na21':'ne21','na22':'ne22','mg22':'ne22','mg23':'na23'}
NU={'n13':.706,'o14':1.05,'o15':1.00}
def mass(n):return pn.Nucleus(n).dm
def A(n):return pn.Nucleus(n).A
def Z(n):return pn.Nucleus(n).Z
def add(a,b,s=1):
    for k,v in b.items():a[k]=a.get(k,0)+s*v

def one_stage(record,fuel):
    raw=record['output_mass_fractions'];negs=sum(-x for x in raw.values() if x<0)
    if negs>1e-7:raise ValueError('Abundance negativity exceeds allowed roundoff')
    norm=sum(max(x,0) for x in raw.values());x={k:max(v,0)/norm for k,v in raw.items()}
    baryons=sum(A(n)*v for n,v in fuel.items())
    hot={k:v*baryons/A(k) for k,v in x.items()}
    cold={};coldnu=0
    for k,v in hot.items():
        dest=DECAY.get(k,k);cold[dest]=cold.get(dest,0)+v
        coldnu+=v*NU.get(k,max(mass(k)-mass(dest),0))
    hotq=sum(mass(k)*v for k,v in fuel.items())-sum(mass(k)*v for k,v in hot.items())
    coldq=sum(mass(k)*v for k,v in hot.items())-sum(mass(k)*v for k,v in cold.items())
    factor=baryons*MU/MEV/record['mass_kg']
    inp=record['initial_energy_J']*factor;hotnu=record['weak_energy_escaped_upper_J']*factor;rad=record['radiation_escaped_J']*factor
    kf=record['final_kinetic_J']/(record['final_kinetic_J']+record['final_internal_J'])
    return dict(fuel=fuel,hot=hot,cold=cold,qhot=hotq,qcold=coldq,nu=hotnu+coldnu,hotnu=hotnu,coldnu=coldnu,input=inp,rad=rad,kf=kf,negative_clipping=negs)

def solve(records,loss):
    stages=[]
    for cat in CATS:
        stage,protons=STAGES[cat];fuel={cat:1,'p':protons}
        mix=records[stage].get('initial_carbon13_number_fraction',0)
        if mix:fuel={'c12':1-mix,'c13':mix,'p':protons}
        stages.append(one_stage(records[stage],fuel))
    stages.append(one_stage(records['He'],{'he4':1}))
    matrix=np.zeros((5,5))
    for j,s in enumerate(stages):
        for i,cat in enumerate(CATS):matrix[i,j]=(1-loss if j<4 else 1)*s['cold'].get(cat,0)-s['fuel'].get(cat,0)
        matrix[4,j]=s['cold'].get('he4',0)-s['fuel'].get('he4',0)
    rhs=np.array([0,0,0,0,1.]);factors=np.linalg.solve(matrix,rhs)
    if min(factors)<0:raise ValueError('Nonpositive operating cycle')
    net={};q=nu=rad=prep=kin=hotinternal=coldheat=0.
    for count,s in zip(factors,stages):
        add(net,s['cold'],count);add(net,s['fuel'],-count)
        q+=count*(s['qhot']+s['qcold']);nu+=count*s['nu'];rad+=count*s['rad'];prep+=count*s['input']
        hotenergy=s['input']+s['qhot']-s['hotnu']-s['rad']
        kin+=count*hotenergy*s['kf'];hotinternal+=count*hotenergy*(1-s['kf']);coldheat+=count*(s['qcold']-s['coldnu'])
    baryon=sum(A(k)*v for k,v in net.items());qcheck=-sum(mass(k)*v for k,v in net.items());charge=sum(Z(k)*v for k,v in net.items())
    heat=q-nu-rad;closed=q-nu
    return {'loss_per_stage':loss,'stage_attempts_per_net_helium':dict(zip(CATS,factors[:4].tolist())),
      'factory_incoming_helium_per_net_helium':float(factors[4]),'new_carbon_per_net_helium':factors[4]*stages[4]['hot'].get('c12',0),
      'hydrogen_atoms_per_net_helium':-net['p'],'nuclear_Q_MeV_per_helium':q,'neutrino_MeV_per_helium':nu,
      'escaped_radiation_MeV_per_helium':rad,'available_heat_MeV_per_helium':heat,'closed_station_heat_MeV_per_helium':closed,
      'fuel_preparation_MeV_per_helium':prep,'factory_input_MeV_per_helium':factors[4]*stages[4]['input'],
      'final_bulk_kinetic_MeV_per_helium':kin,'final_hot_internal_MeV_per_helium':hotinternal,'cold_decay_heat_MeV_per_helium':coldheat,
      'selected_60driver_40heat_99capture_net_MeV_per_helium':.4*.99*(closed+prep)-prep/.6,
      'net_species_per_helium':net,'baryon_residual':baryon,'nuclear_charge_change_per_net_helium':charge,
      'mass_excess_Q_sum_minus_net_MeV':q-qcheck,'energy_partition_residual_MeV':kin+hotinternal+coldheat-heat-prep,
      'balance_matrix_residual':float(np.max(abs(matrix@factors-rhs))), 'balance_matrix_condition':float(np.linalg.cond(matrix))}

def open_external_carbon(records,recovery_c,recovery_h,carbon13_atom_fraction=.0107):
    """Six balances normalized to one helium; unlike factory/catalyst-normalized solve."""
    stages=[one_stage(records[STAGES[c][0]],{c:1,'p':STAGES[c][1]}) for c in CATS]
    # Four catalyst balances, one fresh-H balance, and one product-He balance.
    # Unknowns are stage attempts[4], fresh carbon atoms, and fresh H atoms.
    matrix=np.zeros((6,6));rhs=np.zeros(6);rhs[5]=1
    makeup={'c12':1-carbon13_atom_fraction,'c13':carbon13_atom_fraction}
    for i,c in enumerate(CATS):
        for j,s in enumerate(stages):matrix[i,j]=s['fuel'].get(c,0)-recovery_c*s['cold'].get(c,0)
        matrix[i,4]=-makeup.get(c,0)
    for j,s in enumerate(stages):
        matrix[4,j]=s['fuel'].get('p',0)-recovery_h*s['cold'].get('p',0)
        matrix[5,j]=s['cold'].get('he4',0)
    matrix[4,5]=-1
    solution=np.linalg.solve(matrix,rhs)
    if np.min(solution)<0:raise ValueError('Unphysical open-loop solution')
    attempts,freshc,freshh=solution[:4],solution[4],solution[5]
    incoming={k:freshc*v for k,v in makeup.items()};incoming['p']=freshh
    outgoing={};q=nu=rad=prep=kin=hotinternal=coldheat=0.
    for count,s in zip(attempts,stages):
        for k,v in s['cold'].items():
            r=recovery_c if k in CATS else recovery_h if k=='p' else 0
            outgoing[k]=outgoing.get(k,0)+count*v*(1-r)
        q+=count*(s['qhot']+s['qcold']);nu+=count*s['nu'];rad+=count*s['rad'];prep+=count*s['input']
        e=s['input']+s['qhot']-s['hotnu']-s['rad']
        kin+=count*e*s['kf'];hotinternal+=count*e*(1-s['kf']);coldheat+=count*(s['qcold']-s['coldnu'])
    net=outgoing.copy();add(net,incoming,-1)
    qnet=-sum(mass(k)*v for k,v in net.items());baryon=sum(A(k)*v for k,v in net.items())
    cmass=freshc*sum(pn.Nucleus(k).mass*v for k,v in makeup.items());hmass=freshh*pn.Nucleus('p').mass
    closed=q-nu;w=.4*.99*(closed+prep)-prep/.6
    return {'recovery_c':recovery_c,'recovery_h':recovery_h,'carbon13_atom_fraction':carbon13_atom_fraction,
      'attempts_per_helium':dict(zip(CATS,attempts.tolist())),
      'fresh_carbon_atoms_per_helium':freshc,'fresh_hydrogen_atoms_per_helium':freshh,'fresh_C_to_H_mass_ratio':cmass/hmass,
      'fresh_carbon_kg_per_helium':cmass*MEV/299792458.**2,'fresh_hydrogen_kg_per_helium':hmass*MEV/299792458.**2,
      'closed_station_heat_MeV_per_helium':closed,'preparation_MeV_per_helium':prep,'nuclear_Q_MeV_per_helium':q,
      'neutrino_MeV_per_helium':nu,'prompt_photon_MeV_per_helium':rad,
      'kinetic_MeV_per_helium':kin,'hot_internal_MeV_per_helium':hotinternal,'cold_heat_MeV_per_helium':coldheat,
      'selected_60driver_40heat_99capture_net_MeV_per_helium':w,
      'fresh_C_kg_per_preaux_electric_J':cmass/(299792458.**2*w) if w>0 else None,
      'fresh_H_kg_per_preaux_electric_J':hmass/(299792458.**2*w) if w>0 else None,
      'mass_excess_Q_sum_minus_net_MeV':q-qnet,'baryon_residual':baryon,
      'matrix_residual':float(np.max(abs(matrix@solution-rhs))),
      'external_incoming_atoms_per_helium':incoming,'unreturned_outgoing_atoms_per_helium':outgoing}

def load_records(files,trace_audits=None):
    records={k:json.loads(Path(v).read_text()) for k,v in files.items()}
    for k,p in (trace_audits or {}).items():
        audited=json.loads(Path(p).read_text())['reconstructed_result'];records[k].update(audited)
    return records

def compare(files,trace_audits=None):
    import cycle_ledger as cl
    records=load_records(files,trace_audits);out=[]
    comparison_files={cat:(str(Path(files[st]).resolve()),protons) for cat,(st,protons) in STAGES.items()}
    for loss in [0,1e-5,1e-4,1e-3]:
        ours=solve(records,loss)
        theirs=cl.ledger(loss,files=comparison_files,factory_file=str(Path(files['He']).resolve()))
        comparison={k:ours[k]-theirs[k] for k in ours.keys() & theirs.keys() if isinstance(ours[k],(int,float))}
        for cat in CATS:comparison['attempts_'+cat]=ours['stage_attempts_per_net_helium'][cat]-theirs['stage_attempts_per_net_helium'][cat]
        out.append({'independent':ours,'difference_from_cycle_ledger':comparison})
    return {'scope':'Independent normalization and mass-excess process balance; shared atomic data and explicitly inherited cold-decay neutrino approximation.',
      'files':{k:{'path':str(v),'sha256':hashlib.sha256(Path(v).read_bytes()).hexdigest()} for k,v in files.items()},
      'trace_audits':trace_audits or {},'cases':out,
      'limits':['Cycle arithmetic only; does not establish acceptable viscosity, completed trajectories, separation, collection or plant operation.',
        'Comparisons use original JSON through cycle_ledger; independent results can substitute raw-trace reconstruction.',
        'Hot beta Q assigned entirely to neutrinos; cold O14/N13/O15 mean-neutrino convention inherited and not revalidated here.']}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--files',required=True,help='JSON stage-key to result path mapping');ap.add_argument('--trace-audits');ap.add_argument('--output',required=True);a=ap.parse_args()
    files=json.loads(Path(a.files).read_text());tr=json.loads(Path(a.trace_audits).read_text()) if a.trace_audits else None
    out=compare(files,tr);Path(a.output).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps([{'loss':c['independent']['loss_per_stage'],'H':c['independent']['closed_station_heat_MeV_per_helium'],'U':c['independent']['fuel_preparation_MeV_per_helium'],'max_abs_compare':max(abs(x) for x in c['difference_from_cycle_ledger'].values())} for c in out['cases']],indent=2))
if __name__=='__main__':main()
