"""Bounded bookkeeping and Coulomb-barrier scales, not a stopping simulation."""
import sys,json,math,hashlib
from pathlib import Path
from independent_cycle_audit import BASE,STAGES,load_records,one_stage,MU
home=BASE/'open-loop-redesign';inputs=json.loads((home/'station-inputs.json').read_text());c=inputs['cycle']
files={stage:c['stage_pulse_data'][cat]['source'] for cat,(stage,p) in STAGES.items()};records=load_records(files)
nuclei={};total_baryons=0
for cat,(stage,p) in STAGES.items():
 s=one_stage(records[stage],{cat:1,'p':p});attempt=c['stage_attempts_per_He'][cat]
 for k,v in s['hot'].items():nuclei[k]=nuclei.get(k,0)+attempt*v
 total_baryons+=attempt*({'c12':12,'c13':13,'n14':14,'n15':15}[cat]+p)
N=sum(nuclei.values());NCNO=sum(v for k,v in nuclei.items() if k.startswith(('c','n','o')) and k!='n')
carbon_budget=(.007-c['fresh_C_H_mass_ratio'])*c['fresh_hydrogen_kg_per_He'];atoms=carbon_budget/(12*MU)
kin=c['bulk_kinetic_MeV_per_He'];rows=[]
for e in [1.,10.,100.]:
 residual=total_baryons*e/1000
 rows.append({'wall_entry_keV_per_u':e,'maximum_residual_bulk_MeV_per_He_if_all_baryons_at_ceiling':residual,
  'minimum_fraction_initial_bulk_kinetic_removed':max(0,1-residual/kin),
  'bare_Coulomb_2pi_eta_on_carbon':{s:31.29*6*z/math.sqrt(e) for s,z in [('p',1),('he4',2),('C',6),('N',7),('O',8)]},
  'log10_bare_nonresonant_Gamow_factor':{s:-31.29*6*z/math.sqrt(e)/math.log(10) for s,z in [('p',1),('he4',2),('C',6),('N',7),('O',8)]}})
out={'scope':'Quantitative wall-entry requirements and bare Coulomb-barrier scales. These are not integrated nuclear probabilities or measured reactor performance.',
 'station_inputs_sha256':hashlib.sha256((home/'station-inputs.json').read_bytes()).hexdigest(),
 'structural_C_kg_per_He':carbon_budget,'structural_C_atoms_per_He':atoms,
 'incident_nuclei_per_He':N,'incident_CNO_nuclei_per_He':NCNO,'incident_nuclei_by_species_per_He':nuclei,'baryons_per_He_all_passes':total_baryons,
 'initial_bulk_kinetic_MeV_per_He':kin,'allowable_mean_net_C_destructions_per_incident_nucleus':atoms/N,
 'allowable_mean_net_C_destructions_per_incident_CNO_nucleus_if_only_CNO_causes_loss':atoms/NCNO,
 'hot_tail_fraction_limits_if_destruction_yield_is_given':{str(y):atoms/(N*y) for y in [1,10,100]},
 'wall_entry_energy_scenarios':rows,
 'energy_interpretation':'Electromagnetic braking must transfer directed kinetic energy into captured fields/loads/coolant, not discard it. Retaining99% as heat/electric is not precluded by slowing; no hardware efficiency is calculated.',
 'sources':[
  {'url':'https://arxiv.org/abs/2308.16098','support':'LUNA measured proton captures on C12 and C13 down to60keV center-of-mass; radiative capture remains allowed.'},
  {'url':'https://arxiv.org/abs/2110.00303','support':'Direct C13(alpha,n)O16 data at230–300keV; this exothermic carbon-consuming neutron channel remains open below that measured window.'},
  {'url':'https://arxiv.org/abs/1909.07012','support':'Measured C12+C13 proton branch0.22(7)nb at2.323MeV center-of-mass; deep sub-barrier fusion is very small but nonzero.'},
  {'url':'https://archive.int.washington.edu/PROGRAMS/solar_papers/14N%2Bp/bemmerer.pdf','support':'Primary LUNA paper states sigma=S(E)/E exp(−2πη),2πη=31.29Z1Z2sqrt(mu/E); warns about resonance extrapolation.'},
  {'url':'https://www.nist.gov/pml/stopping-power-range-tables-electrons-protons-and-helium-ions','support':'PSTAR/ASTAR stopping and range data are separate inputs needed for probability integrals; elastic nuclear stopping does not itself mean transmutation.'}],
 'required_future_check':'For each projectile and wall isotope integrate reaction cross sections over the measured slowing spectrum and stopping trajectory, sum irreversible target-carbon changes and unrecycled ejecta, and compare with structural_C_atoms_per_He. Include neutral/high-energy tails, nuclear photons/neutrons, and damaged-material recycling.',
 'limits':['A bare Gamow factor is not a rigorous cross-section bound: S factors, resonances and electronic screening enter.',
 'No common positive kinetic threshold closes all CNO/alpha/proton reactions on natural carbon.',
 'A1keV/u wall-entry ceiling is a candidate deep-sub-barrier requirement, not a demonstrated decelerator or a certified nuclear-loss bound.',
 'The high-energy tail inequality is conditional on its actual destruction yield; values1,10,100 are examples, not estimated yields.',
 'Nuclear stopping in range tables means elastic recoil energy transfer; it must not be counted automatically as nuclear transmutation.']}
Path(sys.argv[1]).write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k] for k in ['structural_C_atoms_per_He','incident_nuclei_per_He','incident_CNO_nuclei_per_He','allowable_mean_net_C_destructions_per_incident_nucleus','hot_tail_fraction_limits_if_destruction_yield_is_given','wall_entry_energy_scenarios']},indent=2))
