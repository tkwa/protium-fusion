"""Closed graphite processing/isotope-pool screen at an archived plant point.

Whole chosen skin treatment costs are scenario inputs, not established process
requirements. Molecular processing does not require atomic vaporization.
"""
import argparse,json,math,hashlib
from pathlib import Path
import material_cycle as mc
ap=argparse.ArgumentParser();ap.add_argument('--plant',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
d=json.loads(a.plant.read_text());d=d['top_detailed'][0] if 'top_detailed' in d else d;p=d['plant'];cy=d['cycle'];area=p['area_m2'];rate=p['helium_per_s'];year=365.25*86400
f13=cy['fresh_carbon_C13_number_fraction'];meanmass=(1-f13)*mc.atomic_kg('c12')+f13*mc.atomic_kg('c13');x13=f13*mc.atomic_kg('c13')/meanmass
isotopes={}
for n in ['c12','c13']:
 draw=sum(v['input_counts'].get(n,0)*cy['stage_attempts_per_He'][cat] for cat,v in cy['stage_pulse_data'].items())*rate
 returned=cy['returned_atoms_per_He'].get(n,0)*rate;fresh=cy['external_input_atoms_per_He'].get(n,0)*rate
 isotopes[n]={'fuel_draw_nuclei_s':draw,'collected_return_nuclei_s':returned,'fresh_supply_nuclei_s':fresh,'pool_balance_nuclei_s':returned+fresh-draw,
  'fuel_draw_kg_s':draw*mc.atomic_kg(n),'collected_return_kg_s':returned*mc.atomic_kg(n),'fresh_supply_kg_s':fresh*mc.atomic_kg(n)}
rows=[]
for thick in [50e-6,100e-6]:
 sigma=2000*thick;flow=area*sigma/year
 rows.append({'thickness_m':thick,'density_kg_m3':2000.,'skin_mass_kg':area*sigma,'annual_processing_kg_s':flow,
  'wall_C13_mass_fraction_if_natural':x13,'wall_C13_inventory_kg_if_natural':area*sigma*x13,
  'natural_C13_through_processing_kg_s':flow*x13,
  'scenarios':[{'electric_J_kg':c,'electric_W':c*flow,'electric_W_m2':c*sigma/year,'fraction_reference_net':c*flow/p['net_W']} for c in [2e6,20e6,100e6]]})
# NIST-JANAF standard-state values at 1000 K, kJ/mol.
Hco,Gco=-111.983,-200.275;Hw,Gw=-247.857,-192.590;Mc=.012011
thermo={'reaction':'C(graphite) + H2O(g) -> CO(g) + H2(g)','temperature_K':1000.,'standard_pressure_Pa':1e5,
 'deltaH_J_per_kg_C':(Hco-Hw)*1000/Mc,'deltaG_J_per_kg_C':(Gco-Gw)*1000/Mc,
 'reverse_standard_deltaG_electric_J_kg_at50pct':-(Gco-Gw)*1000/Mc/.5,
 'sources':['https://janaf.nist.gov/tables/C-093.html','https://janaf.nist.gov/tables/H-064.html'],
 'scope':'Reaction inventory only. Gas activities, kinetics, pressure work, isotope separation and deposition hardware set the actual route. A closed cycle returns H/O reagents; no extra nuclear energy is created, and heat/work are counted once. The tabulated endothermic heat is not all an electrical minimum.'}
checks=[{'name':n+' reservoir isotope balance','passed':abs(v['pool_balance_nuclei_s'])/v['fuel_draw_nuclei_s']<1e-11} for n,v in isotopes.items()]
out={'scope':__doc__,'reference_plant':{'selected_cases':p['selected_cases'],'helium_per_s':rate,'radius_m':p['radius_m'],'net_W':p['net_W']},'isotope_pool':isotopes,'skin_cases':rows,'molecular_route_thermodynamics':thermo,'checks':checks,'sources_sha256':{a.plant.name:mc.hashfile(a.plant)},'code_sha256':mc.hashfile(__file__),'module_sha256':mc.hashfile(mc.__file__)}
a.output.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps({'checks':checks,'isotope_pool':isotopes,'thermodynamics':thermo},indent=2))
