"""Bounded arithmetic, source-hash, energy-accounting, syntax and local-link checks."""
from pathlib import Path
import hashlib,json,math,ast,re
BASE=Path(__file__).resolve().parent
ROOT=Path('/home/tkwa/.codex/worktrees/protium-fusion/website/docs/protium-fusion')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 checks=[]
 def ok(name,p):
  checks.append({'name':name,'passed':bool(p)})
  if not p:raise AssertionError(name)
 s=json.loads((BASE/'spectrum.json').read_text());h=json.loads((BASE/'hot-bypass.json').read_text());c=json.loads((BASE/'containment.json').read_text())
 for code,result in [('spectrum.py',s),('hot_bypass.py',h),('containment.py',c)]:
  ast.parse((BASE/code).read_text());ok('code hash '+code,sha(BASE/code)==result['code_sha256'])
 ok('hot input hash for containment',sha(BASE/'hot-bypass.json')==c['hot_bypass_sha256']);ok('spectrum input hash for hot',sha(BASE/'spectrum.json')==h['source_spectrum_sha256']);ok('spectrum input hash for containment',sha(BASE/'spectrum.json')==c['spectrum_sha256'])
 for rel,actual in [('component-design/nonideal_cycle.json',s['cycle_source_sha256']),('component-design/material_cycle.py',s['material_cycle_sha256']),('pulse_model.py',h['network_sha256']),('nuclear_screen.py',h['nuclear_screen_sha256'])]:ok('dependency hash '+rel,sha(ROOT/rel)==actual)
 for st,v in s['stages'].items():
  ok(st+' source hash',sha(v['source'])==v['sha256'])
  hv=h['stages'][st]
  for ext,digest in hv['sha256'].items():ok(st+' raw hash '+ext,sha(Path(v['source']).with_suffix(ext))==digest)
  ok(st+' sampled net Q quadrature within2%',abs(hv['net_rate_quadrature_relative_to_saved_integrated_Q'])<.02)
  a=hv['results']['all_positive_nonweak_reactions_as_photons'];b=hv['results']['radiative_captures'];ok(st+' loose source >= captures',a['gross_photon_source_J']>=b['gross_photon_source_J'])
  ok(st+' valid Doppler factor',1<=v['max_final_bulk_Doppler_factor']<1.10)
  ok(st+' thermal envelope margin',v['D1p10_endpoint_thermal_tail_screen']['speed_margin_over_sigma']>20)
  for r in a['rows']:
   ok(st+' edge envelope >= uniform '+str(r['kappa_first_interaction_m2_kg']),r['outer_cell_edge_all_outward_first_collision_bypass_J']>=r['uniform_radial_column_source_sensitivity_J'])
   ok(st+' bypass <= source '+str(r['kappa_first_interaction_m2_kg']),r['outer_cell_edge_all_outward_first_collision_bypass_J']<=a['gross_photon_source_J'])
 ok('primary energy sum',math.isclose(sum(x['MeV_per_He'] for x in s['primary_photon_lines']),s['primary_delayed_photon_MeV_per_He'],rel_tol=1e-14))
 ok('positron kinetic sum',math.isclose(sum(x['mean_positron_kinetic_MeV_per_He'] for x in s['parents'].values()),s['main_beta_kinetic_MeV_per_He'],rel_tol=1e-14))
 ok('annihilation count',math.isclose(next(x['photons_per_He'] for x in s['primary_photon_lines'] if x['energy_MeV']==.511),sum(2*s['decay_data'][p]['beta_plus_probability']*d['decays_per_He'] for p,d in s['parents'].items()),rel_tol=1e-14))
 for l in c['primary_lines']:
  ok('mu_en<=mu '+str(l['energy_MeV']),l['mu_en_m2_kg']<=l['mu_m2_kg'])
  ok('degraded floor<=incident '+str(l['energy_MeV']),l['minimum_degraded_energy_mu_en_m2_kg']<=l['mu_en_m2_kg'])
 for r in c['aggregate_localized_cascade_bounds']:
  prev=0
  for b in r['Doppler_boost_scenarios']:
   ok('boost bound monotonic '+str(r['sigma_kg_m2'])+' '+str(b['source_Doppler_factor']),b['plus_saved_prompt_at20MeV_lab_loss_fraction_thermal']>=prev);prev=b['plus_saved_prompt_at20MeV_lab_loss_fraction_thermal']
   ok('remaining allowance closure '+str(r['sigma_kg_m2'])+' '+str(b['source_Doppler_factor']),math.isclose(b['remaining_nonphoton_loss_fraction_for_99pct_capture']+b['plus_saved_prompt_at20MeV_lab_loss_fraction_thermal'],.01,abs_tol=1e-14))
 for r in c['KN_only_transport_sensitivity']:
  tag=str((r['sigma_kg_m2'],r['incident_energy_MeV']))
  ok('MC energy sum '+tag,math.isclose(r['deposited_fraction']+r['combined_leakage_fraction'],1,abs_tol=1e-12))
  ok('MC per-history conservation '+tag,r['max_energy_balance_error_MeV']<1e-12*r['incident_energy_MeV'])
  p=r['uncollided_fraction_exact'];N=r['histories'];ok('MC uncollided analytic '+tag,abs(r['uncollided_fraction_sample']-p)<6*math.sqrt(p*(1-p)/N)+2/N)
 for link in re.findall(r'\]\(([^)]+)\)',(BASE/'README.md').read_text()):
  if '://' not in link:ok('local link '+link,link=='integrity.json' or (BASE/link).exists())
 files=[p for p in BASE.iterdir() if p.is_file() and (p.suffix=='.py' or p.name in ['README.md','spectrum.json','hot-bypass.json','containment.json'])]
 return {'all_passed':all(x['passed'] for x in checks),'checks':checks,'number_of_checks':len(checks),'file_sha256':{p.name:sha(p) for p in sorted(files)},'scope':'Bounded source integrity and arithmetic checks, not validation of omitted coupled photon/electron physics.'}
if __name__=='__main__':
 d=main();(BASE/'integrity.json').write_text(json.dumps(d,indent=2)+'\n');print(json.dumps({'all_passed':d['all_passed'],'checks':d['number_of_checks']},indent=2))
