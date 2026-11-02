"""Final-point energy/material economics, with no reactor operating-point mutation.

Library: evaluate(plant, cycle=None, assumptions=None) -> JSON-safe dict.
CLI selects one explicitly named object from a hashed JSON source. JSON pointers
use RFC 6901 escaping. The selected object is {'plant':..., 'cycle':...} or a plant.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math

HERE = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def select_pointer(document, pointer):
    value = document
    if pointer and not pointer.startswith('/'):
        raise ValueError('Selection must be empty or an explicit /JSON/pointer.')
    for part in pointer.split('/')[1:]:
        token = part.replace('~1', '/').replace('~0', '~')
        value = value[int(token)] if isinstance(value, list) else value[token]
    return value


def default_assumptions():
    return json.loads((HERE / 'assumptions.json').read_text())['assumptions']


def evaluate(plant, cycle=None, assumptions=None):
    """Read a complete plant mass/power ledger; return conditional economic scales.

    Required plant fields: mass_kg, net_W, mass_breakdown_kg (with receiver).
    area_m2, receiver_areal_mass_kg_m2 and net_W_kg are validated when present.
    Flow, refurbishment and driver sections are calculated only if their fields
    are present; unavailable sections remain null instead of inventing values.
    """
    a = default_assumptions()
    if assumptions is not None:
        a.update(assumptions)
    p = plant
    m, power = float(p['mass_kg']), float(p['net_W'])
    if not (m > 0 and power > 0 and math.isfinite(m + power)):
        raise ValueError('Economics requires positive finite installed mass and net export.')
    if p.get('feasible') is False:
        raise ValueError('Selected plant is marked infeasible; select a feasible point explicitly.')
    if not (0 < a['availability'] <= 1):
        raise ValueError('Availability must be in (0,1].')
    if not (0 <= a['receiver_carbon_mass_fraction'] <= 1):
        raise ValueError('Receiver carbon fraction must be in [0,1].')
    for eta in a['methane_chemical_efficiencies']:
        if not (.5 <= eta <= 1):
            raise ValueError('Manufacturing scenarios must honor the specified >=50% efficiency premise.')
    masses = p['mass_breakdown_kg']
    if any(v < 0 or not math.isfinite(v) for v in masses.values()):
        raise ValueError('Mass ledger contains a negative or nonfinite component.')
    checks = {'mass_ledger': math.isclose(sum(masses.values()), m, rel_tol=1e-9)}
    q = power / m
    if 'net_W_kg' in p:
        checks['specific_power'] = math.isclose(q, p['net_W_kg'], rel_tol=1e-9)
    if 'area_m2' in p and 'receiver_areal_mass_kg_m2' in p:
        checks['receiver_area_mass'] = math.isclose(masses['receiver'], p['area_m2'] * p['receiver_areal_mass_kg_m2'], rel_tol=1e-9)
    if not all(checks.values()):
        raise ValueError('Inconsistent selected plant ledger: ' + str(checks))
    yr = a['year_s']
    effective_power = power * a['availability']
    effective_q = effective_power / m
    receiver_C = masses['receiver'] * a['receiver_carbon_mass_fraction']
    resource_H = receiver_C / a['resource_C_H_mass_ratio']
    dg = a['standard_methane_decomposition_G_J_mol'] / a['carbon_molar_mass_kg_mol']
    methane_H = receiver_C * 4 * a['hydrogen_atomic_mass_g_mol'] / 1000 / a['carbon_molar_mass_kg_mol']
    out = {
        'scope': 'Economic evaluation of the explicitly selected plant; all net operating debits and capture/cold-heat routing remain those of the upstream plant record. Chemical scenarios and acquisition/throughput scales are separate.',
        'assumptions': a,
        'plant': {'mass_kg': m, 'net_W': power, 'specific_power_W_kg': q,
                  'installed_kg_per_net_W': 1 / q, 'availability_adjusted_export_W': effective_power,
                  'receiver_areal_mass_kg_m2': p.get('receiver_areal_mass_kg_m2'),
                  'selected_cases': p.get('selected_cases'),
                  'mass_breakdown_kg': masses,
                  'mass_fraction': {k: v / m for k, v in masses.items()},
                  'mass_items_descending': sorted(masses, key=masses.get, reverse=True),
                  'prepared_pulse_W': p.get('prepared_pulse_W'), 'driver_W': p.get('driver_W'),
                  'driver_fraction_of_net_export': p.get('driver_W', 0) / power if 'driver_W' in p else None},
        'embodied_energy': {
            'payback_years_per_MJ_kg_at_declared_availability': 1e6 / effective_q / yr,
            'scenarios': [{'actual_total_J_kg': e, 'total_J': e * m,
                           'equivalent_export_years': e / effective_q / yr}
                          for e in a['whole_plant_embodied_J_kg_scenarios']],
            'lifetime_acceptance': [{'years': t, 'maximum_average_total_J_kg_for_energy_payback': effective_q * yr * t,
                                     'net_kWh_per_installed_kg': effective_q * yr * t / 3.6e6}
                                    for t in a['lifetimes_years']],
            'interpretation': 'Actual total embodied energy includes the chosen construction boundary. Efficiency applies to specified useful positive free-energy changes; it does not assign an embodied-energy cost to every reused kilogram.'},
        'receiver_chemical_step': {
            'reaction': 'CH4(g) -> graphite + 2 H2(g), standard states 298.15 K, 1 bar',
            'receiver_C_kg': receiver_C, 'useful_G_J_kg_C': dg,
            'methane_H_coproduct_kg': methane_H,
            'scenarios': [{'efficiency': eta, 'supplied_J_kg_C': dg / eta,
                           'supplied_total_J': receiver_C * dg / eta,
                           'equivalent_export_years': receiver_C * dg / eta / effective_power / yr}
                          for eta in a['methane_chemical_efficiencies']],
            'interpretation': 'One specified chemical endpoint upgrade, applicable when receiver carbon is obtained from methane. Extraction, delivery, installation and chosen finished structure retain their own accounts. Existing condensed carbon has a different starting endpoint.'},
        'acquisition': {
            'receiver_C_kg': receiver_C, 'receiver_only_carbon_stock': True,
            'resource_C_H_mass_ratio': a['resource_C_H_mass_ratio'],
            'associated_resource_H_kg_at_complete_C_recovery': resource_H,
            'interpretation': 'Associated hydrogen belongs to the accessed resource; it is not automatically transported, stored, consumed or manufactured. Services/reserves require their own elemental bill.'},
        'operating_materials': None, 'refurbishment': None, 'hydrogen_processing': None,
        'construction_throughput': [],
        'solar_comparison': {
            'collector_efficiency': a['comparison_solar_efficiency'],
            'complete_collector_kg_m2': a['comparison_solar_complete_kg_m2'],
            'delivery_efficiency': a['comparison_solar_delivery_efficiency'],
            'mass_specific_power_parity_au': math.sqrt(a['nominal_solar_irradiance_W_m2_at_au'] * a['comparison_solar_efficiency'] * a['comparison_solar_delivery_efficiency'] / a['comparison_solar_complete_kg_m2'] / q),
            'scope': 'Illuminated collector comparison at full output. Siting, duty cycle, storage and delivery apparatus need an equal-service comparison.'},
        'cycle_reference': cycle,
        'validation': checks,
    }
    if 'fresh_H_kg_s' in p and 'fresh_C_kg_s' in p:
        h, c = p['fresh_H_kg_s'], p['fresh_C_kg_s']
        if h <= 0 or c < 0:
            raise ValueError('Fresh feed rates must have H > 0 and C >= 0.')
        ratio = c / h
        if 'fresh_C_H_mass_ratio' in p:
            checks['fresh_C_H_identity'] = math.isclose(ratio, p['fresh_C_H_mass_ratio'], rel_tol=1e-8, abs_tol=1e-15)
        out['operating_materials'] = {
            'fresh_H_kg_s': h, 'fresh_C_kg_s': c, 'fresh_C_H_mass_ratio': ratio,
            'total_007_C_allowance_kg_s': a['operating_fresh_C_H_ceiling'] * h,
            'remaining_C_allowance_kg_s_after_listed_fuel_makeup': a['operating_fresh_C_H_ceiling'] * h - c,
            'original_0056_fuel_allowance_kg_s': a['operating_fuel_C_H_allocation'] * h,
            'original_0014_structural_allowance_kg_s': a['operating_structural_C_H_allocation'] * h,
            'net_J_per_fresh_H_kg': power / h, 'net_GWh_per_fresh_H_kg': power / h / 3.6e12,
            'methane_H_coproduct_years_of_full_rate_H_feed': methane_H / h / yr,
            'fuel_C_methane_chemical_W_scenarios': [{'efficiency': eta, 'W': c * dg / eta} for eta in a['methane_chemical_efficiencies']],
            'lifetime_flows': [{'years': t, 'fresh_H_kg': h * yr * t * a['availability'],
                                'fresh_C_kg': c * yr * t * a['availability'],
                                'receiver_C_over_operating_007_budget': receiver_C / (a['operating_fresh_C_H_ceiling'] * h * yr * t * a['availability'])}
                               for t in a['lifetimes_years']],
            'scope': 'Listed fresh C is the plant material-ledger input. The residual budget is available for separately identified terminal structural makeup. Repair and retained sorting hold-up do not automatically consume that budget.'}
    if 'routine_repaired_kg_s' in p and p['routine_repaired_kg_s'] > 0:
        repaired = p['routine_repaired_kg_s']; included = p.get('replacement_W')
        out['refurbishment'] = {
            'processed_kg_s': repaired, 'annual_processed_fraction_of_installed_mass': repaired * yr / m,
            'included_W': included,
            'included_J_kg': included / repaired if included is not None else None,
            'scenarios': [{'actual_process_J_kg': e, 'W': e * repaired,
                           'fraction_of_reference_export': e * repaired / power,
                           'incremental_W_above_included': e * repaired - included if included is not None else None}
                          for e in a['whole_plant_refurbishment_J_kg_scenarios']],
            'J_kg_consuming_10pct_reference_export': .1 * power / repaired,
            'scope': 'Fixed-throughput sensitivities, not a re-solved plant. Included operating repair is already charged. Processed material can remain in the closed inventory.'}
    if 'recovered_H_kg_s' in p:
        out['hydrogen_processing'] = {
            'recovered_H_kg_s': p['recovered_H_kg_s'],
            'incremental_scenarios': [{'additional_J_kg': e, 'additional_W': e * p['recovered_H_kg_s'],
                                       'fraction_of_reference_export': e * p['recovered_H_kg_s'] / power}
                                      for e in a['incremental_recovered_H_processing_J_kg_scenarios']],
            'scope': 'Use additional process energy above any named existing auxiliary allocation; no per-pass atomization is assumed.'}
    for duration in a['construction_years']:
        flow = m / yr / duration
        out['construction_throughput'].append({
            'construction_years': duration, 'installed_kg_s': flow,
            'receiver_C_kg_s': receiver_C / yr / duration,
            'associated_resource_H_kg_s': resource_H / yr / duration,
            'throughput_multiple_of_installed_repair': flow / p['routine_repaired_kg_s'] if p.get('routine_repaired_kg_s', 0) > 0 else None,
            'external_factory_kg_if_repair_productivity_applies': flow * p['factory_mass_per_repaired_kg_s'] if 'factory_mass_per_repaired_kg_s' in p else None,
            'scope': 'Average construction throughput, conditional on this schedule. Repair productivity is an explicit comparison, not an established new-build specification.'})
    checks['payback_identity'] = all(math.isclose(row['total_J'] / effective_power / yr, row['equivalent_export_years'], rel_tol=1e-12) for row in out['embodied_energy']['scenarios'])
    if not all(checks.values()):
        raise ValueError('Economic validation failed: ' + str(checks))
    out['validation_summary'] = {'passed': sum(checks.values()), 'count': len(checks), 'all_passed': all(checks.values())}
    return out


def render(out):
    p=out['plant']; e=out['embodied_energy']; r=out['receiver_chemical_step']; ac=out['acquisition']; op=out['operating_materials']
    ranking=p['mass_items_descending']; top=ranking[0]
    lines=['# Selected-plant manufacturing and economics ✨','',
           f"This explicitly selected point exports **{p['net_W']:.6g} W** from **{p['mass_kg']:.6g} kg**, or **{1000*p['specific_power_W_kg']:.4g} mW/kg**. It requires **{p['installed_kg_per_net_W']:.4g} installed kg/W**. The largest named mass is `{top}` at **{100*p['mass_fraction'][top]:.3g}%** of the total.",'',
           '## Manufacturing energy','',
           f"At the declared availability {out['assumptions']['availability']:.4g}, every **1 MJ/kg** of actual whole-plant embodied energy takes **{e['payback_years_per_MJ_kg_at_declared_availability']:.4g} years** of export to repay. Manufacturing efficiency applies to a specified useful free-energy change, E = ΔG/η, with η ≥ 0.5. It does not set a universal fabrication cost per kilogram.",'',
           '| Actual total embodied energy | Equivalent export years |','|---:|---:|']
    for row in e['scenarios']:
        lines.append(f"| {row['actual_total_J_kg']/1e6:.4g} MJ/kg | {row['equivalent_export_years']:.4g} |")
    lines+=['',f"The receiver contains **{r['receiver_C_kg']:.6g} kg C** under the stated carbon-fraction assumption. For CH₄ → graphite + 2 H₂ at 298.15 K and 1 bar, the standard useful chemical free-energy increase is **{r['useful_G_J_kg_C']/1e6:.4g} MJ/kg C**. [NIST-JANAF methane](https://janaf.nist.gov/tables/C-067.html).",'']
    for row in r['scenarios']:
        lines.append(f"At {100*row['efficiency']:.4g}% efficiency, this specified chemical step requires {row['supplied_J_kg_C']/1e6:.4g} MJ/kg C and equals {row['equivalent_export_years']:.4g} years of this plant's export.")
    lines+=['','Existing condensed carbon changes the starting chemical endpoint. Acquisition, selected material structure, construction machinery and delivery have separate accounts. Recycled material returning to the same bulk state can have a small net chemical free-energy change.','',
            '## Material access and operating supply','',
            f"At resource C:H = {ac['resource_C_H_mass_ratio']:.4g}, the receiver carbon corresponds to **{ac['associated_resource_H_kg_at_complete_C_recovery']:.6g} kg of associated resource hydrogen**, at complete carbon recovery. That is a resource-access scale; only the selected products need be delivered. Services and reserve stocks need an elemental bill before total construction carbon can be stated."]
    if op:
        lines+=['',f"Operation uses **{op['fresh_H_kg_s']:.6g} kg/s fresh H** and **{op['fresh_C_kg_s']:.6g} kg/s fresh C**, giving **C:H = {op['fresh_C_H_mass_ratio']:.6g}**. The unused part of the 0.007 feed allowance is **{op['remaining_C_allowance_kg_s_after_listed_fuel_makeup']:.6g} kg/s C** before separately identified terminal structural makeup. The former 0.0056 fuel / 0.0014 structural split is an allocation, while the complete material balance determines actual makeup.",'',
                f"Each kilogram of fresh hydrogen supplies **{op['net_GWh_per_fresh_H_kg']:.4g} GWh** of net electricity. Retained skins, sorting tails and material in processing remain inventory; refurbishment throughput alone does not imply fresh carbon consumption."]
    lines+=['','## Factories and recurring work','',
            '| Construction duration | Average installed kg/s | Multiple of installed repair throughput |','|---:|---:|---:|']
    for row in out['construction_throughput']:
        x=row['throughput_multiple_of_installed_repair'];txt=f'{x:.4g}' if x is not None else 'unspecified'
        lines.append(f"| {row['construction_years']:.4g} years | {row['installed_kg_s']:.6g} | {txt} |")
    lines+=['','These are construction schedules supplied as scenarios. Energy efficiency does not establish a factory growth rate or deployment date. Delivery depends on site and trajectory, momentum exchange and braking recovery; a kinetic-energy scale alone is not a dissipative transport charge.']
    if out['refurbishment']:
        f=out['refurbishment'];lines+=['',f"Routine repair handles **{f['processed_kg_s']:.6g} kg/s**, equivalent to **{100*f['annual_processed_fraction_of_installed_mass']:.4g}%** of installed mass annually. The selected plant already charges its included refurbishment load. Fixed-throughput alternatives in the data report only the increment above that load; thermal resizing would require a new plant solve."]
    if out['hydrogen_processing']:
        f=out['hydrogen_processing'];row=next((x for x in f['incremental_scenarios'] if x['additional_J_kg']==1e9),None)
        if row:lines+=['',f"An additional 1 GJ/kg applied to all recovered hydrogen would draw **{100*row['fraction_of_reference_export']:.4g}%** of reference export. This sensitivity does not replace geometric capture or processing-residence requirements."]
    sol=out['solar_comparison'];lines+=['','## Economic comparison','',
        f"With an illuminated collector at {100*sol['collector_efficiency']:.4g}% efficiency, {sol['complete_collector_kg_m2']:.4g} kg/m² complete areal mass and {100*sol['delivery_efficiency']:.4g}% delivery efficiency, collector mass-specific output equals this plant at **{sol['mass_specific_power_parity_au']:.4g} AU**. The irradiance normalization is the IAU nominal 1,361 W/m² at 1 AU. [IAU 2015 Resolution B3](https://www.iau.org/common/Uploaded%20files/IAUGA2015-Resolution-B3-recommended-nominal-conversion.pdf). Equal-service comparison still needs siting, availability, storage and delivery apparatus.",'',
        'Useful roles depend on the value of independent long-lived power and local resource access. Thermodynamic bounds establish energy scales; a money price requires the chosen supply chain, useful life, industrial throughput and opportunity cost. Numerical scenarios here are conditional inputs, not subjective probability intervals.','',
        '## Provenance','',
        'The companion JSON retains the exact plant selection, source hash, economic assumptions and arithmetic checks. Receiver capture, cold rejection, isotope flows and source preparation come unchanged from the selected upstream plant record. This calculation adds no operating debit or mass to that record.','']
    return '\n'.join(lines)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--select',default='',help='Explicit JSON pointer, e.g. /top_detailed/0')
    parser.add_argument('--assumptions',type=Path,default=HERE/'assumptions.json')
    parser.add_argument('--out-dir',type=Path,required=True)
    args=parser.parse_args()
    doc=json.loads(args.input.read_text()); selected=select_pointer(doc,args.select)
    wrapped='plant' in selected;plant=selected['plant'] if wrapped else selected
    cycle=selected.get('cycle') if wrapped else None
    adoc=json.loads(args.assumptions.read_text());out=evaluate(plant,cycle,adoc.get('assumptions',adoc))
    out['provenance']={'source':str(args.input.resolve()),'source_sha256':sha(args.input),'selection_JSON_pointer':args.select,'assumptions_file':str(args.assumptions.resolve()),'assumptions_sha256':sha(args.assumptions),'code_sha256':sha(__file__),'source_scope':doc.get('scope') if isinstance(doc,dict) else None,'primary_sources':adoc.get('primary_sources',{})}
    args.out_dir.mkdir(parents=True,exist_ok=True)
    for fn in ['summary.json','README.md','selected-input.json']:
        if (args.out_dir/fn).exists():
            raise FileExistsError(f'Refusing to overwrite existing research output {args.out_dir/fn}; use a new output directory.')
    (args.out_dir/'selected-input.json').write_text(json.dumps(selected,indent=2,allow_nan=False)+'\n')
    (args.out_dir/'summary.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    (args.out_dir/'README.md').write_text(render(out))
    print(json.dumps({'specific_power_mW_kg':1000*out['plant']['specific_power_W_kg'],'years_per_MJ_kg':out['embodied_energy']['payback_years_per_MJ_kg_at_declared_availability'],'validation':out['validation_summary'],'out_dir':str(args.out_dir)},indent=2))


if __name__=='__main__':
    main()
