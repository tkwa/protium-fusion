"""Explicit B2 spatial-grid and structural-carbon-reserve process scenarios.

Does not overwrite the frozen corrected64 process.json. Other stages remain64.
"""
import argparse
import hashlib
import json
from pathlib import Path
from process import ROOT, load_stages, solve_supply_target

HERE=Path(__file__).resolve().parent
BASE_MAP={
    'c12':'b2-completion/B2-64-compression-cap.json',
    'c13':'b2-completion/regression-C13-compression_cap.json',
    'n14':'b2-completion/regression-A-compression_cap.json',
    'n15':'b2-completion/regression-B1-compression_cap.json',
}
DEFAULT_B2_MAP={
    'b2_64':BASE_MAP['c12'],
    'b2_128_energy_matched':'b2-completion/B2-128-compression-cap-energy-matched.json',
}
FUEL_TARGETS=[.0056,.0063,.007]


def scenario_case(stages,fuel_C_H_target=.0056,total_C_H_budget=.007,**kwargs):
    """The usual flat process case plus explicit net structural-carbon reserve.

    Supply ratios use the SAME fresh H denominator. The structural allowance is
    net unrecoverable carbon, after actual material recovery and reuse. It is
    not permission to equate all replaced receiver material with fresh import.
    Efficiencies, H recovery and fresh C13 fraction pass through to the solver.
    """
    if not 0<=fuel_C_H_target<=total_C_H_budget:
        raise ValueError('Fuel carbon target must fit total carbon budget')
    result=solve_supply_target(stages,target=fuel_C_H_target,**kwargs)
    if not result.get('feasible'):
        return result
    ratio_allowance=total_C_H_budget-fuel_C_H_target
    reserve_kg_per_He=ratio_allowance*result['fresh_hydrogen_kg_per_He']
    total_kg_per_He=total_C_H_budget*result['fresh_hydrogen_kg_per_He']
    electric=result['preaux_electric_J_per_He']
    result.update({
        'fuel_C_H_target':fuel_C_H_target,
        'total_ongoing_C_H_mass_budget':total_C_H_budget,
        'carbon_budget_reserved_fraction':ratio_allowance/total_C_H_budget if total_C_H_budget else 0.,
        'net_structural_C_H_allowance':ratio_allowance,
        'net_structural_carbon_allowance_kg_per_He':reserve_kg_per_He,
        'net_structural_carbon_allowance_kg_per_preaux_J':reserve_kg_per_He/electric if electric>0 else None,
        'total_fresh_carbon_budget_kg_per_He':total_kg_per_He,
        'total_fresh_carbon_budget_kg_per_preaux_J':total_kg_per_He/electric if electric>0 else None,
        'carbon_reserve_boundary':'Steady net unrecoverable carbon from receiver, coatings, substrates and all other material makeup; includes salvage/recycling credits only when physically routed. Startup construction carbon is separate.',
    })
    # Both budgets use the actual new fresh H, not fixed baseline H or gross H.
    assert abs((result['fresh_carbon_kg_per_He']+reserve_kg_per_He)/result['fresh_hydrogen_kg_per_He']-total_C_H_budget)<1e-11
    return result


def build_scenarios(b2_map=None,targets=FUEL_TARGETS):
    families={}
    if b2_map is None:
        b2_map=json.loads((HERE/'b2-process-scenario-map.json').read_text())
    for name,b2_path in b2_map.items():
        paths=dict(BASE_MAP,c12=b2_path)
        stages=load_stages(stage_paths=paths)
        if any(s['trajectory_status']!='completed_corrected_viscosity' for s in stages.values()):
            raise ValueError('Grid scenario requires completed corrected stages: '+name)
        raw=json.loads((ROOT/b2_path).read_text())
        families[name]={
            'scope':'Only B2 spatial resolution varies; A/B1/C13 stay corrected64. Completed time integration does not establish spatial convergence.',
            'stage_map':paths,
            'stage_sha256':{k:s['pulse_data']['sha256'] for k,s in stages.items()},
            'B2_zones':raw['zones'],
            'B2_initial_energy_J':raw['initial_energy_J'],
            'B2_nuclear_gain':raw['deposited_nuclear_over_initial_energy'],
            'B2_desired_product_per_initial_catalyst':raw['desired_product_per_initial_catalyst'],
            'cases':{f'{t:.4f}':scenario_case(stages,t) for t in targets},
        }
    baseline=families['b2_64']
    for name,f in families.items():
        f['comparison_with_b2_64']={}
        for target,point in f['cases'].items():
            ref=baseline['cases'][target]
            f['comparison_with_b2_64'][target]={
                'delta_H_nonnu_MeV_per_He':point['H_nonnu_MeV_per_He']-ref['H_nonnu_MeV_per_He'],
                'delta_U_prepared_MeV_per_He':point['U_prepared_MeV_per_He']-ref['U_prepared_MeV_per_He'],
                'delta_preaux_net_MeV_per_He':point['preaux_net_MeV_per_He']-ref['preaux_net_MeV_per_He'],
                'relative_preaux_net_change':point['preaux_net_MeV_per_He']/ref['preaux_net_MeV_per_He']-1.,
                'delta_required_CNO_return':point['minimum_CNO_return_fraction']-ref['minimum_CNO_return_fraction'],
            }
    if 'b2_256_energy_matched' in families and 'b2_128_energy_matched' in families:
        fine=families['b2_256_energy_matched']
        mid=families['b2_128_energy_matched']
        fine['comparison_with_b2_128_energy_matched']={}
        for target,point in fine['cases'].items():
            ref=mid['cases'][target]
            fine['comparison_with_b2_128_energy_matched'][target]={
                'delta_H_nonnu_MeV_per_He':point['H_nonnu_MeV_per_He']-ref['H_nonnu_MeV_per_He'],
                'delta_U_prepared_MeV_per_He':point['U_prepared_MeV_per_He']-ref['U_prepared_MeV_per_He'],
                'delta_preaux_net_MeV_per_He':point['preaux_net_MeV_per_He']-ref['preaux_net_MeV_per_He'],
                'relative_preaux_net_change':point['preaux_net_MeV_per_He']/ref['preaux_net_MeV_per_He']-1.,
                'delta_required_CNO_return':point['minimum_CNO_return_fraction']-ref['minimum_CNO_return_fraction'],
            }
    return {
        'scope':'Explicit process scenarios; no new apparatus-mass assumptions. Original process.json is frozen as corrected64 provisional family.',
        'scenario_module_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'source_module_sha256':hashlib.sha256((HERE/'process.py').read_bytes()).hexdigest(),
        'base_process_json_sha256':hashlib.sha256((HERE/'process.json').read_bytes()).hexdigest(),
        'total_ongoing_C_H_mass_budget':.007,
        'fuel_C_H_targets':list(targets),
        'families':families,
    }


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--b2-map',type=Path,help='Explicit family name -> B2 JSON map; must include b2_64')
    ap.add_argument('--output',type=Path,default=HERE/'process-scenarios.json')
    args=ap.parse_args()
    result=build_scenarios(json.loads(args.b2_map.read_text()) if args.b2_map else None)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    for name,f in result['families'].items():
        for t,c in f['cases'].items():
            print(name,t,'return',c['minimum_CNO_return_fraction'],'W',c['preaux_net_MeV_per_He'],
                  'structural_C_kg_per_J',c['net_structural_carbon_allowance_kg_per_preaux_J'])
