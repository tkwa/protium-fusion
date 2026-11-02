"""Stationary open-loop CNO process balance, with explicit fresh C and H.

No helium-to-carbon factory. All quantities normalize to one net He-4 atom.
This module supplies process numbers, not apparatus masses or attainable yields.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from cycle_ledger import stage, nuc, MU, MEV

STATES = ['c12', 'c13', 'n14', 'n15']
PROTONS = {'c12': 2, 'c13': 1, 'n14': 1, 'n15': 1}
LEGACY_PATHS = {
    'c12': 'inputs/pulse-profiles/scaled-hot-B2-64-R30-tight.json',
    'c13': 'inputs/pulse-profiles/scaled-hot-C13-64-R30.json',
    'n14': 'inputs/pulse-profiles/scaled-hot-A-64-R30.json',
    'n15': 'inputs/pulse-profiles/scaled-hot-B1-64-R30.json',
}
MASS_MEV_PER_U = 931.49410242
M_H_U = 1.00782503223
M_C13_U = 13.00335483507


def load_stages(root=ROOT, stage_paths=None):
    """Read explicit c12/c13/n14/n15 paths (relative to root or absolute).

    Without explicit paths these are historically invalid legacy diagnostics;
    consumers must inspect trajectory_status before claiming a design result.
    """
    root = Path(root)
    paths = stage_paths or LEGACY_PATHS
    result = {}
    for cat in STATES:
        path = Path(paths[cat])
        if not path.is_absolute():
            path = root / path
        raw = json.loads(path.read_text())
        s = stage(cat, str(path), PROTONS[cat])
        fuel_baryons = sum(nuc(k).A * v for k, v in s['fuel_input'].items())
        if abs(s['fuel_input'].get(cat, 0.) - 1.) > 1e-12:
            raise ValueError('This four-state interface requires isotope-separated pulse inputs')
        effective_end = raw.get('effective_termination_reason', raw.get('termination_reason'))
        viscosity_mode = raw.get('viscosity_mode')
        if viscosity_mode not in ['compression_cap', 'divergence']:
            status = 'legacy_or_unqualified_viscosity'
        elif not raw.get('success') or effective_end != 'final_time' or raw.get('elapsed_sound_times', 0.) < 4.-1e-8:
            status = 'incomplete_or_nonfinal'
        else:
            status = 'completed_corrected_viscosity'
        s['trajectory_status'] = status
        s['pulse_data'] = {
            'source': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
            'trajectory_status': status, 'success': raw.get('success'),
            'viscosity_mode': viscosity_mode,
            'clipped_negative_mass_fraction': s['clipped_negative_mass_fraction'],
            'termination_reason': raw.get('effective_termination_reason', raw.get('termination_reason')),
            'stage': raw.get('stage'), 'mass_kg': raw['mass_kg'],
            'initial_energy_J': raw['initial_energy_J'],
            'attempts_per_pulse': raw['mass_kg'] / (fuel_baryons * MU),
            'fuel_baryons_per_attempt': fuel_baryons,
            'elapsed_s': raw.get('elapsed_s'), 'elapsed_sound_times': raw.get('elapsed_sound_times'),
            'energy_residual_fraction_initial': raw.get('energy_residual_fraction_initial'),
            'initial_surface_radius_m': raw.get('initial_surface_radius_m'),
            'reference_radius_m': raw.get('R0_m'),
        }
        result[cat] = s
    return result


def atomic_kg(name):
    return nuc(name).mass / MASS_MEV_PER_U * MU


def case(stages, C13fraction=.0107, carbon_loss=0., hydrogen_loss=None,
         eta_driver=.6, eta_heat=.4, capture=.99):
    """Stagewise return fractions apply after cold decay and isotope sorting.

    hydrogen_loss=None means coupled H/C losses. Carbon loss applies to all
    four circulating CNO species, not merely elemental carbon. Energy capture
    is independent: unreturned atoms may be collected as heat and then ash.
    """
    if hydrogen_loss is None:
        hydrogen_loss = carbon_loss
    if not all(0 <= x <= 1 for x in [C13fraction, carbon_loss, hydrogen_loss, capture, eta_heat]):
        raise ValueError('Fractions outside [0,1]')
    if not 0 < eta_driver <= 1:
        raise ValueError('Driver efficiency outside (0,1]')
    rC, rH = 1-carbon_loss, 1-hydrogen_loss
    feed = np.array([1-C13fraction, C13fraction, 0., 0.])
    matrix = np.array([[stages[j]['cold_output'].get(i, 0.) for j in STATES] for i in STATES])
    attempts = np.linalg.solve(np.eye(4) - rC*matrix, feed)
    if np.min(attempts) < -1e-9:
        raise ValueError('Nonphysical stationary solution')
    outgoing = {}
    def add(dst, src, factor):
        for name, count in src.items():
            dst[name] = dst.get(name, 0.) + factor*count
    for k, a in zip(STATES, attempts):
        add(outgoing, stages[k]['cold_output'], a)
    gross_H_input = sum(a*stages[k]['protons_input'] for k,a in zip(STATES,attempts))
    gross_H_output = outgoing.get('p',0.)
    fresh_H = gross_H_input - rH*gross_H_output
    helium = outgoing.get('he4',0.)
    if helium <= 0 or fresh_H <= 0:
        raise ValueError('No positive He output / hydrogen supply')
    external_inputs = {'c12': 1-C13fraction, 'c13': C13fraction, 'p': fresh_H}
    returned = {k: outgoing.get(k,0.)*rC for k in STATES}
    returned['p'] = rH*gross_H_output
    ash = {k: v-returned.get(k,0.) for k,v in outgoing.items()}
    ash = {k: v for k,v in ash.items() if v > 0}
    net = external_inputs.copy()
    add(net, ash, -1.)
    def weighted(field):
        return sum(a*stages[k][field] for k,a in zip(STATES,attempts)) / helium
    q, nu, rad, u = (weighted(k) for k in ['Q_total_MeV','nu_MeV','radiated_MeV','input_MeV'])
    h = q-nu
    baryon_residual = sum(nuc(k).A*v for k,v in net.items()) / helium
    mass_q = sum(nuc(k).mass*v for k,v in net.items()) / helium
    input_baryons_per_He = sum(nuc(k).A*v for k,v in external_inputs.items()) / helium
    assert abs(baryon_residual) / input_baryons_per_He < 2e-10, baryon_residual
    assert abs(mass_q-q) / max(1.,abs(q)) < 1e-7, (mass_q,q)
    kinetic, thermal, delayed = 0., 0., 0.
    for k,a in zip(STATES,attempts):
        s=stages[k]
        hot=s['input_MeV']+s['Q_hot_MeV']-s['nu_hot_MeV']-s['radiated_MeV']
        kinetic += a*hot*s['final_kinetic_fraction']/helium
        thermal += a*hot*(1-s['final_kinetic_fraction'])/helium
        delayed += a*(s['Q_cold_MeV']-s['nu_cold_MeV'])/helium
    assert abs(kinetic+thermal+delayed+rad-h-u) < 1e-7*max(1.,h+u)
    mC = ((1-C13fraction)*12 + C13fraction*M_C13_U)*MU/helium
    mH = fresh_H*M_H_U*MU/helium
    preaux = eta_heat*capture*(h+u)-u/eta_driver
    pulse_data = {k: dict(stages[k]['pulse_data']) for k in STATES}
    for k,a in zip(STATES,attempts):
        pulse_data[k]['pulses_per_He'] = float(a/helium/pulse_data[k]['attempts_per_pulse'])
        pulse_data[k]['gross_fuel_kg_per_He'] = float(a/helium*sum(nuc(n).A*v for n,v in stages[k]['fuel_input'].items())*MU)
    return {
        'C13_number_fraction_in_fresh_carbon': C13fraction,
        'carbon_loss_per_cold_stage_output': carbon_loss,
        'hydrogen_loss_per_cold_stage_output': hydrogen_loss,
        'eta_driver': eta_driver, 'eta_heat': eta_heat, 'energy_capture_fraction': capture,
        'stage_attempts_per_He': dict(zip(STATES,(attempts/helium).tolist())),
        'stage_pulse_data': pulse_data,
        'fresh_carbon_atoms_per_He': 1/helium,
        'fresh_hydrogen_atoms_per_He': fresh_H/helium,
        'fresh_carbon_kg_per_He': mC, 'fresh_hydrogen_kg_per_He': mH,
        'fresh_C_H_mass_ratio': mC/mH,
        'fresh_C_H_atom_ratio': 1/fresh_H,
        'nuclear_H_consumed_atoms_per_He': (gross_H_input-gross_H_output)/helium,
        'gross_H_input_atoms_per_He': gross_H_input/helium,
        'unreturned_H_atoms_per_He': hydrogen_loss*gross_H_output/helium,
        'gross_fuel_kg_per_He': sum(d['gross_fuel_kg_per_He'] for d in pulse_data.values()),
        'returned_fuel_atoms_per_He': {k:v/helium for k,v in returned.items()},
        'returned_fuel_kg_per_He': sum(atomic_kg(k)*v for k,v in returned.items())/helium,
        'ash_atoms_per_He': {k:v/helium for k,v in ash.items()},
        'ash_kg_per_He': sum(atomic_kg(k)*v for k,v in ash.items())/helium,
        'nonhelium_ash_kg_per_He': sum(atomic_kg(k)*v for k,v in ash.items() if k!='he4')/helium,
        'nuclear_Q_MeV_per_He': q, 'neutrino_MeV_per_He': nu,
        'H_nonnu_MeV_per_He': h, 'U_prepared_MeV_per_He': u,
        'escaped_prompt_MeV_per_He': rad,
        'bulk_kinetic_MeV_per_He': kinetic, 'hot_internal_MeV_per_He': thermal,
        'delayed_decay_heat_MeV_per_He': delayed,
        'preaux_net_MeV_per_He': preaux,
        'preaux_electric_J_per_He': preaux*MEV,
        'carbon_kg_per_preaux_J': mC/(preaux*MEV) if preaux>0 else None,
        'hydrogen_kg_per_preaux_J': mH/(preaux*MEV) if preaux>0 else None,
        'returned_fuel_kg_per_preaux_J': sum(atomic_kg(k)*v for k,v in returned.items())/helium/(preaux*MEV) if preaux>0 else None,
        'baryon_balance_residual_per_He': baryon_residual,
        'baryon_balance_relative_residual': baryon_residual/input_baryons_per_He,
        'mass_defect_Q_residual_MeV_per_He': mass_q-q,
        'stationary_matrix_condition_number': float(np.linalg.cond(np.eye(4)-rC*matrix)),
        'trajectory_status': {k:stages[k]['trajectory_status'] for k in STATES},
    }


def solve_supply_target(stages, target=.007, C13fraction=.0107, hydrogen_loss=None, **kwargs):
    """Largest stagewise CNO loss meeting a fresh C:H mass-ratio upper bound.

    Fixed hydrogen_loss is independent; None couples H loss to candidate C loss.
    This optimizes no machinery and gives no credit for arbitrary bypass H.
    """
    def evaluate(loss):
        return case(stages,C13fraction,loss,hydrogen_loss,**kwargs)
    low=evaluate(0.)
    if low['fresh_C_H_mass_ratio']>target:
        return {'feasible':False,'reason':'Nuclear catalyst leakage exceeds target at perfect process recovery','perfect_return_case':low}
    if evaluate(1.)['fresh_C_H_mass_ratio']<=target:
        answer=evaluate(1.)
    else:
        loss=brentq(lambda x:evaluate(x)['fresh_C_H_mass_ratio']-target,0.,1.,xtol=1e-14,rtol=1e-13)
        answer=evaluate(loss)
    answer['feasible']=True
    answer['target_fresh_C_H_mass_ratio']=target
    answer['minimum_CNO_return_fraction']=1-answer['carbon_loss_per_cold_stage_output']
    return answer


def build_report(stages):
    selected=solve_supply_target(stages)
    return {
        'scope':'Stationary external-carbon process; no local helium-to-carbon factory. Hardware, supply extraction, separation energy, scheduling and multidimensional inflow are outside this balance.',
        'supply_constraint':'Fresh carbon / fresh protium MASS ratio <=0.007 kg/kg, including actual unreturned hydrogen; no hydrogen bypass dilution.',
        'stage_paths': {k:s['source'] for k,s in stages.items()},
        'stage_sha256': {k:s['pulse_data']['sha256'] for k,s in stages.items()},
        'trajectory_status': {k:s['trajectory_status'] for k,s in stages.items()},
        'selected': selected,
        'recovery_scan': [case(stages,carbon_loss=1-r) for r in [0.,.9,.99,.999,.9999,.99999,1.]],
        'H_recovery_sensitivity': [solve_supply_target(stages,hydrogen_loss=1-r) for r in [0.,.5,.9,.99,.999,.9999,.99999,1.]],
        'C13feed_sensitivity': [solve_supply_target(stages,C13fraction=x) for x in [0.,.0096,.0107,.0116,1.]],
        'atomic_masses_u': {'protium':M_H_U,'c12':12.,'c13':M_C13_U},
        'sources': ['https://www.ciaaw.org/carbon.htm','https://www.ciaaw.org/hydrogen.htm',
                    'https://pml.nist.gov/cgi-bin/Compositions/stand_alone.pl?ascii=ascii&ele=C',
                    'https://pml.nist.gov/cgi-bin/Compositions/stand_alone.pl?ascii=ascii&ele=H'],
    }


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--stage-map',type=Path,help='JSON mapping c12/c13/n14/n15 to corrected stage JSON paths')
    ap.add_argument('--output',type=Path,default=Path(__file__).with_name('process.json'))
    args=ap.parse_args()
    stages=load_stages(stage_paths=json.loads(args.stage_map.read_text()) if args.stage_map else None)
    report=build_report(stages)
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    d=report['selected']
    print(json.dumps({k:d[k] for k in ['minimum_CNO_return_fraction','fresh_C_H_mass_ratio','H_nonnu_MeV_per_He','U_prepared_MeV_per_He','preaux_net_MeV_per_He','carbon_kg_per_preaux_J']},indent=2))
