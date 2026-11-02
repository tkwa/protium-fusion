"""Energy-only membrane renewal sensitivities with installed factory mass retained.

A same-material recycling cost is explicit dissipation/handling work. It is not
inferred from a positive-free-energy manufacturing efficiency. No frame, cooler,
extra factory throughput, lifetime or material-return performance is qualified.
"""
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / 'open-loop-redesign/operating-cases.json'
b = json.loads(SOURCE.read_text())['cases']['high_source_rating_surface']
YEAR = 365.25 * 86400
area = 4 * math.pi * b['radius_m']**2
routine_flow = b['mass_kg'] * b['replacement_fraction_per_year'] / YEAR
cadences = {
    'large_pulses': sum(1/b['pulse_intervals_s'][k] for k in ['c12', 'c13', 'n14']),
    'all_pulses': sum(1/t for t in b['pulse_intervals_s'].values()),
}
rows = []
for material, density, carbon_fraction in [('silica', 2200., 0.), ('carbon', 2000., 1.)]:
    for thickness in [1e-8, 1e-7, 1e-6]:
        film_mass = area * density * thickness
        for cadence, frequency in cadences.items():
            flow = film_mass * frequency
            for routine_MJ in [100., 2.]:
                before_film = b['net_W'] + b['replacement_W'] - routine_flow * routine_MJ * 1e6
                row = {
                    'material': material, 'density_kg_m3': density, 'thickness_m': thickness,
                    'cadence': cadence, 'renewals_per_year': frequency * YEAR,
                    'installed_film_mass_kg': film_mass, 'extra_recycling_kg_s': flow,
                    'one_film_processing_time_s_at_original_total_factory_flow': film_mass / routine_flow,
                    'B1_interval_s_for_deadline_comparison': b['pulse_intervals_s']['n15'],
                    'routine_effective_total_MJ_kg': routine_MJ,
                    'net_before_extra_film_W': before_film,
                    'max_film_recycling_MJ_kg_for_1mWkg_energy_only':
                        (before_film - .001 * (b['mass_kg'] + film_mass)) / flow / 1e6,
                    'required_total_factory_productivity_kg_per_kg_s_if_original_factory_mass':
                        (routine_flow + flow) / b['mass_breakdown_kg']['replacement_manufacturing'],
                    'total_factory_own_mass_processing_time_days':
                        b['mass_breakdown_kg']['replacement_manufacturing'] / (routine_flow + flow) / 86400,
                    'carbon_escape_fraction_ceiling_if_entire_structural_reserve_assigned_to_film':
                        b['net_material_C_allowance_kg_s'] / flow if carbon_fraction else None,
                    'noncarbon_makeup_kg_s_if_unreturned_fraction_1e8_inverse':
                        flow * (1 - carbon_fraction) * 1e-8,
                    'cost_rows': [],
                }
                for film_MJ in [1., 2., 3.1, 10., 100.]:
                    renewal_W = flow * film_MJ * 1e6
                    net = before_film - renewal_W
                    residual = (net + renewal_W - before_film) / b['gross_W']
                    assert abs(residual) < 1e-14
                    row['cost_rows'].append({
                        'film_effective_total_MJ_kg': film_MJ,
                        'extra_recycling_W': renewal_W,
                        'extra_recycling_W_m2': renewal_W / area,
                        'energy_only_net_W': net,
                        'energy_only_net_W_kg_including_film': net / (b['mass_kg'] + film_mass),
                        'power_balance_relative_residual': residual,
                    })
                rows.append(row)
# The 100 MJ/kg routine case reproduces the previous output before extra renewal.
assert abs((b['net_W'] + b['replacement_W'] - routine_flow * 1e8)/b['net_W']-1) < 1e-14
out = {
    'scope': __doc__, 'source_sha256': hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
    'code_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'radius_m': b['radius_m'], 'surface_area_m2': area,
    'baseline_mass_kg': b['mass_kg'], 'routine_recycling_kg_s': routine_flow,
    'retained_installed_factory_mass_kg': b['mass_breakdown_kg']['replacement_manufacturing'],
    'structural_carbon_makeup_allowance_kg_s': b['net_material_C_allowance_kg_s'],
    'fresh_C_H_budget': .007,
    'rows': rows,
    'limitations': [
        'All original hardware is retained; only the bare film is added to the denominator. Support, cooling, additional heat rejection and any larger factory remain unallocated.',
        'Productivity values are requirements imposed on the retained factory, not demonstrated factory ratings.',
        'Mean throughput does not ensure completion between pulses; renewal deadlines, batch buffers and peak handling capacity require a schedule.',
        'Effective recycling costs are explicit sensitivity assumptions; no complete process recipe selects one.',
        'A film recycled after a pulse must function during that pulse. Reprocessing does not establish irradiation lifetime or leak tightness.',
        'Material removed toward the central cavity could worsen contamination. Receiver-side collection is a separate design requirement.',
        'The entire carbon allowance is shared with all other equipment. Silicon and oxygen losses need their own supply account.',
    ],
}
(HERE/'renewal-energy.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
print(json.dumps({'case_count':len(rows),'power_balances_passed':True}))
