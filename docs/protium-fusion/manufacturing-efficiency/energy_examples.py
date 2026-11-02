"""Thermodynamic examples and declared factory-cost sensitivities.

No reactor inputs are modified. The user's >=50% premise applies to a specified
positive free-energy upgrade including that conversion's dissipation. The
zero-net-chemistry examples instead expose physical operations and losses.
"""
import hashlib
import importlib.util
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
YEAR = 365.25 * 86400
M_C = 0.0120107  # kg/mol; same natural-carbon bookkeeping convention as prior study
RGAS = 8.31446261815324


def hashfile(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    inputs_path = ROOT / 'open-loop-redesign/station-inputs.json'
    cases_path = ROOT / 'open-loop-redesign/operating-cases.json'
    station_path = ROOT / 'open-loop-redesign/station.py'
    inp = json.loads(inputs_path.read_text())
    case = json.loads(cases_path.read_text())['cases']['high_source_rating_surface']
    spec = importlib.util.spec_from_file_location('factory_station_sensitivity', station_path)
    station = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(station)
    # Explicit transcription of NIST-JANAF rows, standard pressure 1 bar.
    raw = {
        'carbon_gas_298': {'dHf_J_mol':716670., 'dGf_J_mol':671244.,
                         'source':'https://janaf.nist.gov/tables/C-003.html'},
        'methane_gas_298': {'dHf_J_mol':-74873., 'dGf_J_mol':-50768.,
                           'source':'https://janaf.nist.gov/tables/C-067.html'},
        'graphite_heating': {'T_initial_K':298.15, 'T_final_K':1000.,
                             'dH_J_mol':11795., 'S_initial_J_mol_K':5.740,
                             'S_final_J_mol_K':24.457,
                             'source':'https://janaf.nist.gov/tables/C-002.html'},
        'graphite_basal_cleavage': {'cleavage_J_m2':.39, 'uncertainty_J_m2':.02,
                                  'selected_density_kg_m3':2000.,
                                  'source':'https://arxiv.org/abs/1506.00536',
                                  'scope':'Ideal AB stacking estimate from measured incommensurate cleavage plus modeled twist boundary; room-temperature basal surfaces.'},
    }
    ch4 = raw['methane_gas_298']
    dg = -ch4['dGf_J_mol']/M_C
    methane = {
        'reaction':'CH4(g) -> C(graphite) + 2 H2(g), all standard states at 298.15 K',
        'dG_J_kg_C':dg, 'dH_J_kg_C':-ch4['dHf_J_mol']/M_C,
        'input_usable_J_kg_C_at_efficiency': {str(e):dg/e for e in [.5,.9,.99,1.]},
        'alternative_fixed_dissipation_J_kg_C': {
            str(d): {'input_J_kg_C':dg+d, 'efficiency':dg/(dg+d)} for d in [1e5,2e6,3.1e6]},
        'note':'Efficiency and fixed dissipation are alternative bookkeeping choices for the same conversion, not additive penalties. H2 is retained as a useful product; oxidizing it requires a separately counted oxidant and changes the reaction.'
    }
    atom = raw['carbon_gas_298']
    atomization = {'dG_J_kg_C':atom['dGf_J_mol']/M_C,
                   'dH_J_kg_C':atom['dHf_J_mol']/M_C,
                   'input_at_50pct_J_kg_C':2*atom['dGf_J_mol']/M_C,
                   'note':'A specified high-energy monatomic gas endpoint, not a compulsory graphite recycling step.'}
    heat = raw['graphite_heating']
    heating = {'heat_J_kg':heat['dH_J_mol']/M_C,
               'availability_relative_298K_J_kg':(heat['dH_J_mol']-heat['T_initial_K']*(heat['S_final_J_mol_K']-heat['S_initial_J_mol_K']))/M_C,
               'note':'Availability and heat are distinct. A loop retaining hot graphite does not repeat the full cold-to-hot step every pass.'}
    cleavage = raw['graphite_basal_cleavage']
    surfaces = [{'plate_thickness_m':d, 'retained_surface_energy_J_kg':cleavage['cleavage_J_m2']/(cleavage['selected_density_kg_m3']*d)} for d in [2.5,1e-3,1e-6,1e-8]]
    x13=.0107
    mixing = [{'T_K':t,'ideal_complete_C12_C13_separation_J_kg':-RGAS*t*((1-x13)*math.log(1-x13)+x13*math.log(x13))/M_C} for t in [298.15,1000.]]
    rows=[]
    for cost in [0.,1e5,1e6,2e6,3.1e6,1e7,1e8]:
        out=station.evaluate(case['helium_per_s'],case['radius_m'],inp['cycle'],inp['receiver'],
            radiator_T=case['radiator_K'],hot_T=case['receiver_base_K'],
            absorber_sigma=case['absorber_sigma_kg_m2'],driver_peak_W_kg=case['driver_peak_W_kg'],
            service_sigma=case['service_sigma_kg_m2'],timing=case['timing_basis'],replacement_J_kg=cost)
        assert out['feasible']
        assert abs(out['mass_relative_residual'])<1e-10
        assert abs(out['electrical_relative_residual'])<1e-10
        assert abs(out['global_energy_relative_residual'])<1e-10
        if cost==1e8:
            assert abs(out['mass_kg']/case['mass_kg']-1)<1e-12
            assert abs(out['net_W']/case['net_W']-1)<1e-12
        rows.append({'replacement_energy_J_kg':cost, 'mass_kg':out['mass_kg'],
                     'net_W':out['net_W'], 'net_W_kg':out['net_W_kg'],
                     'replacement_W':out['replacement_W'],
                     'radiator_area_fraction':out['radiator_area_fraction'],
                     'net_increase_fraction':out['net_W']/case['net_W']-1,
                     'fabrication_arithmetic_years_if_delivered_feed_is_already_suitable':cost*out['mass_kg']/out['net_W']/YEAR})
    receiver_mass=case['mass_breakdown_kg']['receiver']
    receiver_flow=receiver_mass*case['replacement_fraction_per_year']/YEAR
    receiver_only={
        'scope':'Fixed-installed-mass diagnostic only; other component maintenance remains 100 MJ/kg.',
        'refurbished_receiver_kg_s':receiver_flow,
        'old_receiver_refurbishment_W':receiver_flow*1e8,
        'new_receiver_refurbishment_W_at_2MJkg':receiver_flow*2e6,
        'saved_W':receiver_flow*(1e8-2e6),
        'net_W_at_fixed_mass':case['net_W']+receiver_flow*(1e8-2e6),
        'new_same_material_processing_allowance_J_kg':2e6,
        'receiver_startup_methane_upgrade_J_at_50pct':receiver_mass*2*dg,
        'receiver_startup_methane_upgrade_years_of_old_net_output_at_50pct':receiver_mass*2*dg/case['net_W']/YEAR,
        'net_material_C_allowance_kg_s':case['net_material_C_allowance_kg_s'],
        'makeup_from_methane_at_50pct_W_for_net_material_allowance':case['net_material_C_allowance_kg_s']*2*dg,
        'note':'Fresh synthesis applies to net new material, whereas refurbishment applies to gross returned material. 2 MJ/kg is a sensitivity, not a graphite-specific measured cost.'}
    transport={'selected_speed_m_s':30000.,'kinetic_J_kg':.5*30000.**2,
               'source_electric_J_kg_if_separately_selected_50pct_acceleration':.5*30000.**2/.5,
               'note':'Separate transport example; manufacturing premise does not fix propulsion efficiency. Net journey cost also depends on braking recovery and orbital changes.'}
    output={'scope':'Independent arithmetic and conditional cost sensitivity; no future manufacturing hardware is demonstrated by these calculations.',
            'premise':'For a declared positive useful free-energy increase, supplied usable energy = deltaG/eta with eta >=0.5; all losses inside that conversion boundary are already included. Chemical energy consumed must be counted. Zero-net graphite chemistry needs explicit dissipation rather than this ratio.',
            'input_sha256':{str(p.relative_to(ROOT)):hashfile(p) for p in [inputs_path,cases_path,station_path]},
            'code_sha256':hashfile(Path(__file__)), 'transcribed_data':raw,
            'methane_to_graphite':methane,'carbon_atomization':atomization,
            'graphite_heating':heating,'basal_surface_examples':surfaces,'ideal_isotope_separation_examples':mixing,
            'sensitivity_scope':'Secondary mass-reallocation calculation using the old 1 kg/W manufacturing mass and resized radiators. Prefer sensitivity.json for fixed installed hardware; lower energy does not independently establish lower factory mass.',
            'initial_cost_scope':'The same-cost fabrication arithmetic assumes delivered feed already has suitable chemical free energy. In canon, methane-to-graphite upgrading must be added once for fresh receiver carbon; a 2 MJ/kg refurbishment cost cannot replace that initial chemistry charge.',
            'secondary_mass_reallocation_sensitivity':rows,'receiver_only_fixed_mass_sensitivity':receiver_only,
            'transport_example':transport,
            'checks':{'baseline_reproduced':True,'all_mass_and_energy_residuals_below_1e10_inverse':True}}
    path=HERE/'energy-examples.json'
    path.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'methane':methane,'heating':heating,'surfaces':surfaces,'mixing':mixing,'station':rows,'receiver':receiver_only},indent=2))

if __name__=='__main__':
    main()
