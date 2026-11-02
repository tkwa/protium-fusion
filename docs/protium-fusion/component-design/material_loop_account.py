"""Finite material-return, separation and decay inventories; frozen pulse inputs.

No impurity hydrodynamics or component performance is implied. JSON output
separates exact ledger arithmetic, ideal thermodynamic floors, and choices.
"""
import argparse, json, hashlib, math, sys
from pathlib import Path
import numpy as np
ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
a=ap.parse_args();R=a.root.resolve();sys.path.insert(0,str(R));sys.path.insert(0,str(R/'open-loop-redesign'))
import process as p
paths=json.loads((R/'open-loop-redesign/corrected-stage-map.json').read_text());st=p.load_stages(R,paths)
c=json.loads((R/'open-loop-redesign/station-inputs.json').read_text())['cycle']
station=json.loads((R/'cold-sink-redesign/cases.json').read_text())['cases']['working_200K_60pct_IR11_bypass1']
old=json.loads((R/'open-loop-redesign/operating-cases.json').read_text())['cases']['high_source_rating_surface']
loading=json.loads((R/'membrane-receiver/pulse-loading.json').read_text())
H=station['helium_per_s'];MU=p.MU;KB=1.380649e-23;AREA=station['area_m2'];checks=[]
def check(name,truth):
 checks.append({'name':name,'passed':bool(truth)})
 if not truth:raise AssertionError(name)
def mass(n):return p.atomic_kg(n)
def element(n):return ''.join(x for x in n if x.isalpha()) if n not in ['p','n'] else ('h' if n=='p' else 'neutron')
def entropy(counts):
 total=sum(counts.values())
 return -sum(v*math.log(v/total) for v in counts.values() if v>0) if total>0 else 0.
flows={};combined={};bypassed={};mixed={};minsep=0.;minisotope=0.
for cat in p.STATES:
 s=st[cat];attempt=c['stage_attempts_per_He'][cat]*H;out={k:v*attempt for k,v in s['cold_output'].items()}
 inp={k:v*attempt for k,v in s['fuel_input'].items()};groups={}
 for n,v in out.items():
  combined[n]=combined.get(n,0)+v;groups.setdefault(element(n),{})[n]=v
 mixfloor=KB*300*entropy(out);isofloor=KB*300*sum(entropy(g) for g in groups.values())
 # Two principal isotope-cut feeds; pure stage-tagged streams remain separate.
 cutgroups={el:g for el,g in groups.items() if el in ['c','n'] and min([v/sum(g.values()) for v in g.values()] or [0])>0}
 # Report all groups, no arbitrary discard or purity classification in accounting.
 flows[cat]={'stage':s['pulse_data']['stage'],'pulse_interval_s':old['pulse_intervals_s'][cat],
 'input_kg_s':{k:v*mass(k) for k,v in inp.items()},'cold_output_kg_s':{k:v*mass(k) for k,v in out.items()},
 'cold_number_fraction_by_element':{el:{k:v/sum(g.values()) for k,v in g.items()} for el,g in groups.items()},
 'cold_elements_kg_s':{el:sum(v*mass(k) for k,v in g.items()) for el,g in groups.items()},
 'ideal_300K_atomic_all_species_separation_W':mixfloor,'ideal_300K_atomic_isotope_separation_W':isofloor,
 'single_pulse_input_kg':{k:v*mass(k)*old['pulse_intervals_s'][cat] for k,v in inp.items()}}
 minsep+=mixfloor;minisotope+=isofloor
CNO=sum(v*mass(n) for n,v in combined.items() if n in p.STATES)
catloss=c['carbon_loss_per_cold_stage_output'];struct=old['net_material_C_allowance_kg_s']
# Seven explicit interface allocations; sum of probabilities is only first order.
loss_alloc={'catcher_escape':150e-6,'chemical_transfer':50e-6,'isotope_system_leak':50e-6,'assay_purge':20e-6,'storage_transfer':20e-6,'launch_transfer':30e-6,'unallocated':1-(1-catloss)/np.prod([1-x for x in [150e-6,50e-6,50e-6,20e-6,20e-6,30e-6]])}
check('allocated serial survival exactly matches target',abs(np.prod([1-v for v in loss_alloc.values()])-(1-catloss))<1e-14)
# If a wrong-route fraction eps has extra probability b of removing a CNO
# nucleus from the reusable quartet, it consumes eps*b of each pass budget.
# This relation is a bound; b requires a tracer/burn calculation.
impurity=[{'extra_nonquartet_probability_per_impurity':b,'impurity_for_50ppm_extra_nonquartet_diversion':50e-6/b} for b in [1,.1,.01,.001]]
# Pure inert dilution: each intended successful attempt requires 1/(1-eps)
# preparation; full internal recirculation retains catalyst. Local burn may differ.
u=c['U_prepared_MeV_per_He'];hh=c['H_nonnu_MeV_per_He'];base=.6*.99*(hh+u)-u/.6
inert=[]
for eps in [.0001,.001,.01,.03,.05]:
 up=u/(1-eps);w=.6*.99*(hh+up)-up/.6
 inert.append({'wrong_inert_fraction':eps,'prepared_multiplier':1/(1-eps),'preaux_MeV_He':w,'preaux_relative_loss':(base-w)/base})
# Radiation half lives are preserved legacy network values; primary current
# TUNL recommendations shown separately to expose their negligible difference.
HL={'n13':597.9978026337243,'o14':70.59978360562287,'o15':121.99951623253882}
decay={}
for k,hl in HL.items():
 meanmass=0.;wallmassrate=0.;birthmassrate=0.;birthpower=0.
 for cat,s in st.items():
  tag=s['pulse_data']['stage'];v=loading['stages'][tag];sp=loading['species'].index(k)
  packet=v['spectral_packets'];mf=np.asarray(packet['mass_fluence_kg_m2']);ya=np.asarray(packet['Y_at_wall']);A=loading['A'][sp]
  wallmass=AREA*float(mf@ya[:,sp])*A
  birth= s['pulse_data']['mass_kg'] * json.loads(Path(s['source']).read_text())['output_mass_fractions'].get(k,0.)
  freq=1/old['pulse_intervals_s'][cat];wallmassrate+=wallmass*freq;birthmassrate+=birth*freq;meanmass+=wallmass*freq*hl/math.log(2)
 decay[k]={'network_half_life_s':hl,'birth_kg_s':birthmassrate,'wall_arrival_kg_s':wallmassrate,
 'fraction_surviving_to_wall':wallmassrate/birthmassrate,'mean_radioactive_inventory_after_wall_kg':meanmass,
 'additional_hold_s_for_remaining_fraction':{str(f):hl*math.log2(1/f) for f in [.01,.001,.0001,.00001]},
 'fraction_remaining_after_3h':2**(-10800/hl)}
# Isotope rotor equilibrium contrast is illustrative, no guessed SWU throughput.
rotors=[]
for T in [100,300]:
 for speed in [300,700,1000]:
  alpha=math.exp(MU*speed**2/(2*KB*T))
  rotors.append({'T_K':T,'tip_speed_m_s':speed,'one_u_radial_equilibrium_factor':alpha,
   'ideal_countercurrent_equilibrium_steps_50to99pct':math.log(99)/math.log(alpha),
   'scope':'Equilibrium radial factor; not a measured stage factor or real cascade length.'})
# Two-outlet isotope selector with Gaussian mass response; assignment error is
# numerical condition on separation/sigma, not an achieved hardware tolerance.
from scipy.special import ndtri
selector=[]
for e in [.01,.001,.0001,.00001]:
 sig=1/(2*ndtri(1-e));selector.append({'equal_abundance_wrong_cut_fraction':e,'sigma_mass_u_max':sig,
 'mass_resolution_m_over_FWHM_at14':14/(2.35482*sig),
 'sigma_r_over_radius_at14':sig/28})
cost=[]
for Jkg in [2e6,1e8,1e9,1e10,1e11]:
 W=old['all_stage_flow_kg_s']*Jkg
 cost.append({'specific_electric_J_kg_gross':Jkg,'W':W,'W_m2':W/AREA,'fraction_of_current_export':W/station['net_W']})
check('flow equals frozen baryonic stream within atomic mass offset',abs(sum(sum(s['input_kg_s'].values()) for s in flows.values())/old['all_stage_flow_kg_s']-1)<.001)
check('process return finite',0<catloss<.001)
check('all reported physical material outputs nonnegative',all(v>=0 for s in flows.values() for v in s['cold_output_kg_s'].values()))
result={'scope':'Historical maximum-loss budget-boundary arithmetic, ideal separation bounds and unused illustrative allocations. Current perfect-eventual-return and retained-reservoir accounting is in material-reservoir.json. No achieved or required physical loss/purity/hardware rate implied.',
 'station':{'helium_per_s':H,'area_m2':AREA,'net_W':station['net_W'],'net_W_kg':station['net_W_kg'],'gross_baryonic_fuel_kg_s':old['all_stage_flow_kg_s'],
 'fuel_fresh_C_H':.0056,'structural_C_H':.0014,'total_C_H':.007,'CNO_return_budget_boundary':1-catloss,'CNO_output_kg_s':CNO,'CNO_processing_loss_kg_s_at_budget_boundary':catloss*CNO,
 'fresh_C_kg_s':old['feed_C_kg_s'],'fresh_H_kg_s':old['feed_H_kg_s'],'net_structure_C_allowance_kg_s':struct},
 'stage_flows':flows,'combined_cold_output_kg_s':{k:v*mass(k) for k,v in combined.items()},
 'separation_floors':{'temperature_K':300,'atomic_species_fixed_stage_batch_W':minsep,'isotopes_within_elements_fixed_stage_batch_W':minisotope,
 'atomic_species_globally_mixed_W':KB*300*entropy(combined),'at50pct_minimum_free_energy_efficiency_electric_W':2*minsep,
 'meaning':'Ideal atomic mixture entropy change at the stated T; molecular chemistry, compression, irreversibility and hardware are separate.'},
 'illustrative_unused_serial_loss_allocation':loss_alloc,'impurity_extra_loss_bounds':impurity,'inert_impurity_energy_sensitivity':inert,'decay_inventory':decay,
 'rotor_equilibrium_examples':rotors,'mass_selector_gaussian_acceptance':selector,'gross_processing_cost_scan':cost,
 'frozen_one_day_working_inventory_kg':station['mass_breakdown_kg']['working_inventory'],
 'sources':{'TUNL_half_lives':'https://nucldata.tunl.duke.edu/nucldata/HalfLife.shtml','light_isotope_centrifuge_simulation':'https://doi.org/10.1016/S0168-9002(00)00266-7'},
 'source_sha256':{str(f.relative_to(R)):hashlib.sha256(f.read_bytes()).hexdigest() for f in [R/'open-loop-redesign/station-inputs.json',R/'open-loop-redesign/operating-cases.json',R/'cold-sink-redesign/cases.json',R/'membrane-receiver/pulse-loading.json']},
 'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'checks':checks}
a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
print(json.dumps({k:result[k] for k in ['station','separation_floors','decay_inventory','checks']},indent=2))
