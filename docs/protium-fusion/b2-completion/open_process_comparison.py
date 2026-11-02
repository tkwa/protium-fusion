"""Cross-check production open-loop report with independent 6x6 normalization."""
import argparse,json,sys,hashlib
from pathlib import Path
from independent_cycle_audit import STAGES,load_records,open_external_carbon
FIELDS={
 'fresh_carbon_atoms_per_helium':'fresh_carbon_atoms_per_He',
 'fresh_hydrogen_atoms_per_helium':'fresh_hydrogen_atoms_per_He',
 'fresh_C_to_H_mass_ratio':'fresh_C_H_mass_ratio',
 'fresh_carbon_kg_per_helium':'fresh_carbon_kg_per_He',
 'fresh_hydrogen_kg_per_helium':'fresh_hydrogen_kg_per_He',
 'closed_station_heat_MeV_per_helium':'H_nonnu_MeV_per_He',
 'preparation_MeV_per_helium':'U_prepared_MeV_per_He',
 'nuclear_Q_MeV_per_helium':'nuclear_Q_MeV_per_He',
 'neutrino_MeV_per_helium':'neutrino_MeV_per_He',
 'prompt_photon_MeV_per_helium':'escaped_prompt_MeV_per_He',
 'kinetic_MeV_per_helium':'bulk_kinetic_MeV_per_He',
 'hot_internal_MeV_per_helium':'hot_internal_MeV_per_He',
 'cold_heat_MeV_per_helium':'delayed_decay_heat_MeV_per_He',
 'selected_60driver_40heat_99capture_net_MeV_per_helium':'preaux_net_MeV_per_He',
 'fresh_C_kg_per_preaux_electric_J':'carbon_kg_per_preaux_J',
 'fresh_H_kg_per_preaux_electric_J':'hydrogen_kg_per_preaux_J'}

def compare_report(path):
 d=json.loads(Path(path).read_text())
 files={stage:d['stage_paths'][cat] for cat,(stage,n) in STAGES.items()}
 records=load_records(files);cases=[('selected',d['selected'])]
 for group in ['recovery_scan','H_recovery_sensitivity','C13feed_sensitivity']:
  cases.extend((group+'/'+str(i),v) for i,v in enumerate(d[group]))
 rows=[]
 for label,t in cases:
  if t.get('feasible',True) is False:
   rows.append({'case':label,'skipped':'Production solver reports infeasible; no operating-point comparison.'});continue
  assert (t['eta_driver'],t['eta_heat'],t['energy_capture_fraction'])==(.6,.4,.99),'Set independent efficiency formula to match deliberate new scenario.'
  ours=open_external_carbon(records,1-t['carbon_loss_per_cold_stage_output'],1-t['hydrogen_loss_per_cold_stage_output'],t['C13_number_fraction_in_fresh_carbon'])
  rel={}
  for key,other in FIELDS.items():
   if ours[key] is None or t[other] is None:
    assert ours[key] is None and t[other] is None;continue
   rel[key]=(ours[key]-t[other])/max(abs(ours[key]),abs(t[other]),1e-300)
  for cat,v in ours['attempts_per_helium'].items():rel['attempts_'+cat]=(v-t['stage_attempts_per_He'][cat])/max(abs(v),abs(t['stage_attempts_per_He'][cat]),1e-300)
  checks={'field_comparison_relative_below_1e_8':max(abs(x) for x in rel.values())<1e-8,
          'net_external_mass_Q_closes':abs(ours['mass_excess_Q_sum_minus_net_MeV'])<1e-8*max(abs(ours['nuclear_Q_MeV_per_helium']),1),
          'baryons_close':abs(ours['baryon_residual'])<1e-10*(ours['fresh_hydrogen_atoms_per_helium']+13*ours['fresh_carbon_atoms_per_helium']),
          'stationarity_closes':ours['matrix_residual']<1e-9}
  checks={k:bool(v) for k,v in checks.items()}
  rows.append({'case':label,'relative_differences':rel,'checks':checks,'independent':ours})
 return {'scope':'Independent six-variable fresh-supply/He-normalized solve versus four-state factory-free production process report. Shared pulse outputs and atomic data.','source_report':str(path),'sha256':hashlib.sha256(Path(path).read_bytes()).hexdigest(),
   'trajectory_status':d['trajectory_status'],'cases':rows,
   'checks_passed':sum(sum(x['checks'].values()) for x in rows if 'checks' in x),
   'checks_failed':sum(sum(not v for v in x['checks'].values()) for x in rows if 'checks' in x),
   'limits':['Matching arithmetic does not qualify pulse trajectories, chemical isotope recovery, collector nuclear survival or plant performance.',
    'Production constants use tabulated atomic masses in u; independent conversion uses mass-energy/c^2. Their small rounding difference is included in comparisons.']}
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--process-report',required=True);a.add_argument('--output',required=True);p=a.parse_args()
 r=compare_report(p.process_report);Path(p.output).write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
 print(json.dumps({k:r[k] for k in ['trajectory_status','checks_passed','checks_failed']},indent=2))
