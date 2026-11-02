"""Closed material reservoirs: incomplete extraction creates stock, not loss.

Returns are eventual fractions. An extraction pass can leave retained residue;
only explicitly irreversible escape/disposal enters fresh-makeup accounting.
Nuclear side products are separately held for optional nuclear reprocessing.
"""
import argparse,json,math,hashlib
from pathlib import Path
import material_cycle as mc
ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();R=a.root.resolve()
d=json.loads((R/'component-design/material-cycle-validation.json').read_text());st=d['loaded_stages'];p=mc.solve_cycle(st,carbon_return=1.,hydrogen_return=1.);station=json.loads((R/'cold-sink-redesign/cases.json').read_text())['cases']['working_200K_60pct_IR11_bypass1'];H=station['helium_per_s'];area=station['area_m2'];DAY=86400.;YEAR=365.25*DAY
qgross=H*p['gross_fuel_kg_per_He'];coldout={}
for key,stage in st.items():mc.add(coldout,stage['cold_output'],p['stage_attempts_per_He'][key])
q=H*sum(v*mc.atomic_kg(n) for n,v in coldout.items() if n in mc.DEFAULT_STATES);checks=[]
def check(n,b):
 checks.append({'name':n,'passed':bool(b)})
 if not b:raise AssertionError(n)
def extraction(c,ell):
 """Irreversible loss first, extraction of c of surviving material, residue retries."""
 denom=1-(1-ell)*(1-c)
 return {'eventual_return':(1-ell)*c/denom,'eventual_irreversible_loss':ell/denom,'expected_processing_passes':1/denom}
rows=[]
for c in [.1,.5,.9,.99]:
 for ell in [0,1e-6,1e-5]:
  v=extraction(c,ell);check('eventualclosure '+str((c,ell)),abs(v['eventual_return']+v['eventual_irreversible_loss']-1)<1e-14)
  if ell==0:check('no invented loss '+str(c),abs(v['eventual_return']-1)<1e-14)
  rows.append({'extracted_fraction_per_pass':c,'irreversible_loss_per_pass':ell,**v,
   'monthly_cycle_residence_s_from_random_arrival':30*DAY*(v['expected_processing_passes']-.5),
   'monthly_cycle_mean_inventory_kg':q*30*DAY*(v['expected_processing_passes']-.5)})
skins=[]
for name,rho,thick in [('graphite10um',2000,1e-5),('graphite50um',2000,5e-5),('alumina50um',3900,5e-5),('Al100um',2300,1e-4),('Sn100um',7000,1e-4)]:
 sigma=rho*thick
 for days in [1,30,365.25]:
  period=days*DAY;flow=area*sigma/period;loading=q/area*period
  skins.append({'material_case':name,'rho_kg_m3':rho,'thickness_m':thick,'sigma_kg_m2':sigma,'batch_interval_days':days,
   'gross_skin_processed_kg_s':flow,'skin_over_fuel_throughput':flow/q,'end_of_batch_CNO_loading_kg_m2':loading,
   'CNO_mass_fraction_relative_skin_at_end':loading/sigma,'mean_uniform_age_trapped_CNO_kg':q*period/2,
   'cold_start_one_complete_period_CNO_stock_kg':q*period,
   'electric_W_by_J_kg':{str(cost):flow*cost for cost in [2e6,3e7,1e8]},
   'electric_W_m2_by_J_kg':{str(cost):sigma*cost/period for cost in [2e6,3e7,1e8]},
   'scope':'C/N-bearing retained fuel only; H/He require a separate prompt local gas loop. Whole chosen skin batch processed, no assertion that all skin must be chemically reduced or atomized. Costs are explicit scenarios; stock is recyclable fuel, not fresh ongoing loss.'})
# The active quartet excludes these nuclei, but the receiver still retains them.
side={k:v for k,v in p['ash_atoms_per_He'].items() if k not in (*mc.DEFAULT_STATES,'p','he4') and v>0}
sideflow={k:{'atoms_per_He':v,'kg_s':v*H*mc.atomic_kg(k),'one_year_inventory_kg':v*H*mc.atomic_kg(k)*YEAR} for k,v in side.items()}
O16=side.get('o16',0.);qNO=mc.mass_MeV({'o16':1,'p':2})-mc.mass_MeV({'n14':1,'he4':1})
closure={'stoichiometry':'O16 + 2p -> N14 + He4 + positron/neutrino (net atomic-mass convention)',
 'candidate_channels':['O16(p,gamma)F17','F17 beta+ O17','O17(p,alpha)N14'],
 'Q_atomic_MeV_per_restored_O16':qNO,'candidate_O16_rate_atoms_per_He':O16,
 'additional_protons_per_existing_He_if_all_O16_restored':2*O16,'additional_He_per_existing_He_if_all_O16_restored':O16,
 'scope':'Stoichiometric opportunity only. Requires its own rates, residence/preparation and apparatus. No energy or recovered-catalyst credit taken in the present station.'}
# Standard-state alumina reduction free-energy inventory, thermal cycle separate.
G298=1582.275e3/.1019613;G1000=1361.437e3/.1019613
out={'scope':__doc__,'perfect_eventual_return_cycle':p,'frozen_He_rate_s':H,'perfect_return_fresh_C_kg_s':p['fresh_carbon_kg_per_He']*H,
 'perfect_return_fresh_H_kg_s':p['fresh_hydrogen_kg_per_He']*H,'perfect_return_gross_fuel_kg_s':qgross,'retained_CNO_output_kg_s':q,'prompt_unburned_H_output_kg_s':coldout.get('p',0)*H*mc.atomic_kg('p'),'prompt_He_product_kg_s':coldout.get('he4',0)*H*mc.atomic_kg('he4'),
 'available_C_H_budget':.007,'carbon_budget_interpretation':'Available supply ceiling; no requirement to consume all of it. Ordinary retained fuel, rejects and catchers are inventory, not physical escape.',
 'extraction_cases':rows,'skin_batch_cases':skins,'held_nuclear_sideproducts':sideflow,'optional_O16_nuclear_recovery':closure,
 'alumina_thermodynamic_reference':{'source':'https://janaf.nist.gov/tables/Al-096.html','deltaG_reduction_J_kg_298K':G298,'deltaG_reduction_J_kg_1000K':G1000,
  '50pct_free_energy_conversion_electric_J_kg_1000K':2*G1000,'scope':'Full Al2O3->Al+O2 reduction at standard chemical activities, if chosen. It is not a mandatory regeneration step; sensible heat and thermal integration are separate.'},
 'checks':checks,'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'module_sha256':mc.hashfile(mc.__file__)}
a.output.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps({k:out[k] for k in ['perfect_return_fresh_C_kg_s','perfect_return_fresh_H_kg_s','perfect_return_gross_fuel_kg_s','held_nuclear_sideproducts','optional_O16_nuclear_recovery']},indent=2))
