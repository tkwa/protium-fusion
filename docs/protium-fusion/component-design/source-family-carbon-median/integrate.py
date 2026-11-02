"""Explicit source64 versus nominal64 material/plant comparison.
Keeps 5000 and 2000 kg/m2 receiver accounting as separate conditional cases.
No hydro; all physical records, profiles, compositions and durations stay fixed.
"""
from pathlib import Path
import sys,json,hashlib,math,argparse
P=Path(__file__).resolve().parent;ROOT=P.parent.parent
sys.path[:0]=[str(ROOT/'mixture-design'),str(ROOT/'component-design')]
import variant_integrate as vi
import material_cycle as mc
sha=vi.sha

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--work-dir',type=Path,default=Path('/home/tkwa/src/protium-fusion/output/source-family-integration'));ap.add_argument('--replay',type=Path);a=ap.parse_args();a.work_dir.mkdir(parents=True,exist_ok=True)
 receiver=vi.module('source_family_receiver',ROOT/'mixture-design/receiver.py');plant=vi.module('source_family_plant',ROOT/'mixture-design/plant_compare.py')
 codes=[Path(__file__),ROOT/'mixture-design/variant_integrate.py',ROOT/'mixture-design/plant_compare.py',ROOT/'mixture-design/receiver.py',ROOT/'component-design/material_cycle.py',ROOT/'cold-sink-redesign/station.py',ROOT/'open-loop-redesign/operating-cases.json',ROOT/'cold-sink-redesign/cases.json']
 hashes={str(p):sha(p) for p in codes};receiver.self_check();checks=[];sources={};families={}
 nom={}
 for cat,tag,h in [('n14','A','3'),('n15','B1','2'),('c12','B2','2'),('c13','C13','4')]:
  f=ROOT/f'mixture-design/rate-results/carbon64/{tag}-H{h}-fixed_energy-z64-k0p001-t1-rates-updated-median.json';d=vi.read(f);i=f.with_suffix('.inputs.json');ip=f.with_suffix('.initial.json');tr=Path(d['trace_path']);inp=vi.read(i)
  nom[cat]={'record':str(f.relative_to(ROOT)),'initial':str(ip.relative_to(ROOT)),'inputs':str(i.relative_to(ROOT)),'trace':str(tr),'duration_s':inp['duration_s'],'initial_Y':inp['initial_Y'],'eos':'FastEOS','result_sha256':sha(f),'initial_profile_sha256':sha(ip),'input_record_sha256':sha(i),'raw_trace_sha256':sha(tr)}
 nompath=P/'nominal64-stage-map.json';vi.write(nompath,{'family':'nominal64-A3-B1H2-B2H2-C13H4-carbon-median','stage_map':nom})
 if a.replay:
  old=vi.read(a.replay);archives={k:v['stage_inputs'] for k,v in old['families'].items()};sources[str(a.replay)]=sha(a.replay)
 for key,mapfile in [('nominal64',nompath),('source64',P/'stage-map.json')]:
  if a.replay:
   archive=archives[key];selected=[archive[k]['selected'] for k in vi.CAT];cs=[{'name':key+' explicit reduced-input replay','passed':True}];ss={}
  else:selected,archive,cs,ss=vi.extract_explicit(mapfile,ROOT,a.work_dir/key,receiver,mc)
  for x in selected:
   d=x['record'];r=d['rate_update'];assert d['zones']==64 and d['viscosity_mode']=='compression_cap' and d['opacity_m2_kg']==.001
   assert r['table_sha256']=='9166934c7f82868daf7ff95054c78f6078994ab1c0a8e53040a37c1d651f4b99' and r['mode']=='updated' and r['quantile']=='median'
  cycle=plant.cycle_for(selected,rc=1.,rh=1.)
  cases={}
  for sigma in [5000.,2000.]:
   pp=plant.best_rate(selected,cycle,radius=1.1e9,absorber_sigma=sigma);cs+=vi.result_checks(selected,cycle,pp);cases[str(int(sigma))]=pp
  families[key]={'stage_inputs':archive,'cycle':cycle,'capacity_plants':cases};checks.extend([dict(c,family=key) for c in cs]);sources.update(ss)
 nominal=families['nominal64'];changed=families['source64'];stagecompare={}
 for cat in vi.CAT:
  n=nominal['stage_inputs'][cat]['selected'];s=changed['stage_inputs'][cat]['selected'];nd=n['record'];sd=s['record']
  deltas={k:sd[k]/nd[k]-1 for k in ['mass_kg','initial_energy_J','elapsed_s','deposited_nuclear_over_initial_energy','desired_product_per_initial_catalyst']}
  for k in ['mass_kg','initial_energy_J','elapsed_s']:
   checks.append({'name':cat+' source/nominal '+k+' match','value':deltas[k],'passed':abs(deltas[k])<1e-11})
  checks.append({'name':cat+' source/nominal composition match','passed':s['stage_data']['initial_Y']==n['stage_data']['initial_Y']})
  stagecompare[cat]={'relative_changes':deltas,'nominal_gain':nd['deposited_nuclear_over_initial_energy'],'source_gain':sd['deposited_nuclear_over_initial_energy'],'nominal_desired_product':nd['desired_product_per_initial_catalyst'],'source_desired_product':sd['desired_product_per_initial_catalyst'],'nominal_source_rating_W':n['receiver']['source_rating']['peak_total_prepared_electric_W']/1.1,'source_profile_rating_W':s['receiver']['source_rating']['peak_total_prepared_electric_W']/1.1,'nominal_occupied_s':n['receiver']['timing_at_radius']['mass999']['occupied_s']*1.1,'source_occupied_s':s['receiver']['timing_at_radius']['mass999']['occupied_s']*1.1,'nominal_rate_kernel_sha256':nd['rate_update']['kernel_sha256'],'source_rate_kernel_sha256':sd['rate_update']['kernel_sha256']}
 comparisons={}
 selected=[changed['stage_inputs'][k]['selected'] for k in vi.CAT]
 for sig in ['5000','2000']:
  np=nominal['capacity_plants'][sig];sp=changed['capacity_plants'][sig];same=plant.station(selected,changed['cycle'],np['helium_per_s'],radius=1.1e9,absorber_sigma=float(sig));checks+=vi.result_checks(selected,changed['cycle'],same)
  comparisons[sig]={'source_at_nominal_He_rate':same,'capacity_relative_changes':{k:sp[k]/np[k]-1 for k in ['net_W','mass_kg','net_W_kg','helium_per_s','peak_source_electric_W','sequential_occupancy_mass999']},'same_He_rate_relative_changes':{k:same[k]/np[k]-1 for k in ['net_W','mass_kg','net_W_kg']}}
 cd={k:changed['cycle'][k]/nominal['cycle'][k]-1 for k in ['U_prepared_MeV_per_He','H_nonnu_MeV_per_He','fresh_C_H_mass_ratio','gross_fuel_kg_per_He','gross_H_output_atoms_per_He']}
 for f,h in hashes.items():assert sha(f)==h,('consumed code changed',f)
 for c in checks:c['passed']=bool(c['passed'])
 assert all(c['passed'] for c in checks)
 out={'schema':1,'scope':__doc__,'families':families,'stage_comparison':stagecompare,'cycle_relative_changes':cd,'receiver_scenario_comparisons':comparisons,'parameters':{'radius_m':1.1e9,'radiator_K':200.,'thermal_package_efficiency':.6,'driver_efficiency':.6,'assumed_energy_capture':.99,'eventual_ordinary_material_return':1.,'receiver_areal_mass_scenarios_kg_m2':[5000.,2000.],'maximum_main_window_occupancy':.8,'maximum_radiator_area_fraction':.9},'rate_source_note':'Nominal64 kernel56e07d differs from source649e3252 only by initialization and an assertion that the two carbon tables have identical temperature grids. Exact archived code diff and identical table establish numerical equivalence; hashes remain distinct. Other190 rates unchanged.','source_timing_scope':'The receiver adapter reconstructs source power and 99.9percent main-flight timing from each exact prescribed hot profile. It does not prove the actual cold acceleration/neutralization/collision waveform. B1 cold energy-delivery proxy is separately audited; controlled collision heat is converted stored kinetic energy, not an extra external flash heater. Slow seed delivery, overlap, finite outlets and 3D assembly remain explicit component tasks.','receiver_scope':'Both receiver masses use the same .99 capture assumption only for a controlled mass comparison. The2000 case does not inherit qualified capture, stopping, thermal cycling or mechanical support from5000; these require the separate receiver-thickness/photon/support packets.','checks':checks,'source_sha256':sources,'code_sha256':hashes}
 vi.write(P/'integration.json',out)
 print(json.dumps({'checks':len(checks),'stages':stagecompare,'cycle_relative_changes':cd,'plants':{k:{sig:{q:v[q] for q in ['feasible','net_W_kg','net_W','mass_kg','helium_per_s','fresh_C_H_mass_ratio','sequential_occupancy_mass999','radiator_area_fraction']} for sig,v in f['capacity_plants'].items()} for k,f in families.items()},'comparisons':{k:v['capacity_relative_changes'] for k,v in comparisons.items()}},indent=2))
if __name__=='__main__':main()
