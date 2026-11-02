"""Energy and material economics from an explicitly selected plant ledger.

All cost ranges are declared scenarios, not inferred probability intervals.
The script does not alter reactor output or double-charge existing auxiliaries.
"""
from pathlib import Path
import argparse, hashlib, json, math

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def calculate(data):
    a=data['assumptions'];year=a['year_s'];rows={}
    for name,case in data['cases'].items():
        p=case['plant'];m=p['mass_kg'];power=p['net_W'];q=power/m
        c=p['mass_breakdown_kg']['receiver'];h=p['fresh_H_kg_s']
        assert math.isclose(q,p['net_W_kg'],rel_tol=1e-12)
        assert math.isclose(c,p['area_m2']*p['receiver_areal_mass_kg_m2'],rel_tol=1e-12)
        assert math.isclose(sum(p['mass_breakdown_kg'].values()),m,rel_tol=1e-12)
        dg=a['standard_methane_decomposition_G_J_mol']/a['carbon_molar_mass_kg_mol']
        byproduct=4*a['hydrogen_atomic_mass_g_mol']/1000/a['carbon_molar_mass_kg_mol']
        chemicals=[]
        for eta in a['methane_chemical_efficiencies']:
            e=dg/eta*c
            chemicals.append({'efficiency':eta,'supplied_J_kg_C':dg/eta,'receiver_chemical_step_J':e,'receiver_chemical_step_equivalent_export_years':e/power/year,'operating_fresh_C_chemical_W':p['fresh_C_kg_s']*dg/eta,'example_days_if_external_construction_power_1e24W':e/1e24/86400})
        lifetime=[]
        for life in a['lifetimes_years']:
            mh=h*year*life
            lifetime.append({'years':life,'fresh_H_kg':mh,'fresh_C_kg':p['fresh_C_kg_s']*year*life,'receiver_C_per_consumed_H':c/mh,'receiver_C_over_operating_007_budget':c/(a['operating_fresh_C_H_ceiling']*mh),'maximum_average_embodied_J_kg_for_energy_payback':q*year*life,'delivered_kWh_per_installed_kg':q*year*life/3.6e6})
        refurbishment=[]
        for e in a['whole_plant_refurbishment_J_kg_scenarios']:
            load=e*p['routine_repaired_kg_s']
            refurbishment.append({'J_kg':e,'W':load,'fraction_of_reference_export':load/power,'incremental_W_above_included_refurbishment':0. if math.isclose(load,p['replacement_W'],rel_tol=1e-12) else load-p['replacement_W']})
        build=[]
        for duration in a['construction_years']:
            rate=m/year/duration
            build.append({'years':duration,'mean_installed_kg_s':rate,'mean_receiver_C_kg_s':c/year/duration,'mean_H_resource_access_kg_s_at_chosen_abundance':c/a['resource_C_H_mass_ratio']/year/duration,'throughput_vs_installed_routine_factory':rate/p['routine_repaired_kg_s'],'factory_kg_if_repair_productivity_also_applies_to_construction':rate*p['factory_mass_per_repaired_kg_s'],'scope':'Construction-period throughput illustration; repair productivity is not a demonstrated new-build factory specification.'})
        solar=q and math.sqrt(a['nominal_solar_irradiance_W_m2_at_au']*a['comparison_solar_efficiency']*a['comparison_solar_delivery_efficiency']/a['comparison_solar_complete_kg_m2']/q)
        rows[name]={'scope':case['scope'],'specific_power_W_kg':q,'installed_kg_per_net_W':1/q,'net_J_per_consumed_H_kg':power/h,'net_GWh_per_consumed_H_kg':power/h/3.6e12,'mass_fraction':{k:v/m for k,v in p['mass_breakdown_kg'].items()},'receiver_C_kg':c,'resource_H_kg_accessed_for_receiver_at_chosen_C_H':c/a['resource_C_H_mass_ratio'],'methane_H_byproduct_kg':c*byproduct,'methane_H_byproduct_years_of_reactor_H_consumption':c*byproduct/h/year,'whole_plant_embodied_scenarios':[{'J_kg':e,'energy_J':e*m,'equivalent_export_years':e/q/year} for e in a['whole_plant_embodied_J_kg_scenarios']],'methane_chemical_step':chemicals,'lifetimes':lifetime,'refurbishment_scenarios':refurbishment,'refurbishment_J_kg_consuming_10pct_export':.1*power/p['routine_repaired_kg_s'],'H_processing_at_1GJ_kg_W':p['recovered_H_kg_s']*1e9,'H_processing_at_1GJ_kg_fraction_export':p['recovered_H_kg_s']*1e9/power,'source_electric_fraction_of_export':p['driver_W']/power,'construction_scenarios':build,'direct_solar_mass_parity_distance_au':solar,'solar_specific_power_at_5au_W_kg':a['nominal_solar_irradiance_W_m2_at_au']*a['comparison_solar_efficiency']*a['comparison_solar_delivery_efficiency']/a['comparison_solar_complete_kg_m2']/25,'export_fraction_of_nominal_solar_luminosity':power/a['nominal_solar_luminosity_W'],'validation':{'specific_power_identity':True,'receiver_areal_mass_identity':True,'mass_ledger_identity':True}}
        assert math.isclose(refurbishment[0]['W'],p['replacement_W'],rel_tol=1e-12)
        assert math.isclose(rows[name]['whole_plant_embodied_scenarios'][0]['energy_J']/power,rows[name]['whole_plant_embodied_scenarios'][0]['equivalent_export_years']*year,rel_tol=1e-12)
        rows[name]['validation']['existing_refurbishment_reproduced']=True
        rows[name]['validation']['payback_energy_identity']=True
    return {'scope':__doc__,'assumptions':a,'primary_sources':data['primary_sources'],'cases':rows}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--input',type=Path,default=Path(__file__).with_name('economics-inputs.json'));parser.add_argument('--out',type=Path,default=Path(__file__).with_suffix('.json'));args=parser.parse_args()
    data=json.loads(args.input.read_text());out=calculate(data)
    out['input_sha256']=sha(args.input);out['code_sha256']=sha(__file__)
    args.out.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:{'kg_W':v['installed_kg_per_net_W'],'payback_years_per_MJ_kg':v['whole_plant_embodied_scenarios'][0]['equivalent_export_years'],'receiver_C_kg':v['receiver_C_kg'],'methane_step_years_50pct':v['methane_chemical_step'][0]['receiver_chemical_step_equivalent_export_years'],'solar_parity_au':v['direct_solar_mass_parity_distance_au']} for k,v in out['cases'].items()},indent=2))
if __name__=='__main__':main()
