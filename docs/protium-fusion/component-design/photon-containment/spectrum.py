"""Inventory-weighted primary delayed photons and conservative extra-radiation reserves."""
from pathlib import Path
import hashlib,json,sys,math
ROOT=Path('/home/tkwa/.codex/worktrees/protium-fusion/website/docs/protium-fusion');sys.path.insert(0,str(ROOT/'component-design'))
from material_cycle import load_stage,cold_decay,mass_MeV,DEFAULT_DECAYS
BASE=Path('/home/tkwa/src/protium-fusion/output/mixture-runs')
HSET={'A':3,'B1':2,'B2':2.5,'C13':4}
DATA={'n13':{'beta_plus_probability':.99803,'beta_mean_MeV':.49221,'Q_atomic_MeV':2.22047,'half_life_s':9.9589*60,'lines_MeV_per_decay':[],'source':'https://www.lnhb.fr/nuclides/N-13_tables.pdf','evaluation':'LNHB updated2026; ground-state transition'},'o15':{'beta_plus_probability':.999001,'beta_mean_MeV':.73347,'Q_atomic_MeV':2.75418,'half_life_s':122.266,'lines_MeV_per_decay':[],'source':'https://www.lnhb.fr/nuclides/O-15_tables.pdf','evaluation':'LNHB updated2026; ground-state transition'},'o14':{'beta_plus_probability':(.019+99.249+.61)/100,'beta_mean_weighted_MeV':(.019*.06437+99.249*.77055+.61*1.87595)/100,'Q_atomic_MeV':5.14304,'half_life_s':70.606,'lines_MeV_per_decay':[[1.63520,.00052],[2.312593,.99388],[3.94750,.0000211]],'source':'https://www.nndc.bnl.gov/ensnds/14/N/ec_decay.pdf','evaluation':'ENSDF July1990 evaluation; gamma intensities perparent, including cascades'}}
def main():
 cy=json.loads((ROOT/'component-design/nonideal_cycle.json').read_text());sel=next(r['cases']['control'] for r in cy['cases'] if r['mixture_H_per_catalyst']==HSET)
 stages={};counts={};source={}
 for s,h in HSET.items():
  tag=format(h,'.8g').replace('.','p');paths=list(BASE.rglob(f'{s}-H{tag}-fixed_energy-z64-k0p001-t1.json'));assert len(paths)==1
  st=load_stage(paths[0]);record=json.loads(paths[0].read_text());vmax=max(abs(v) for v in record['profile']['velocity_m_s']);bulk_beta=vmax/299792458.;w=sel['stage_attempts_per_He'][s];stages[s]={'attempts_per_He':w,'attempts_per_pulse':st['attempts_per_pulse'],'max_final_bulk_speed_m_s':vmax,'max_final_bulk_Doppler_factor':((1+bulk_beta)/(1-bulk_beta))**.5,'max_endpoint_temperature_K':max(record['profile']['T_K']),'hot_output':st['hot_output'],'source':str(paths[0]),'sha256':hashlib.sha256(paths[0].read_bytes()).hexdigest()}
  Denv=1.10;benv=(Denv**2-1)/(Denv**2+1);sig=math.sqrt(1.380649e-23*max(record['profile']['T_K'])/(13*1.66053906660e-27))/299792458.;z=(benv-bulk_beta)/sig
  stages[s]['D1p10_endpoint_thermal_tail_screen']={'minimum_parent_A':13,'one_dimension_thermal_speed_sigma_over_c':sig,'speed_margin_over_sigma':z,'Maxwell_speed_tail_probability_upper':math.erfc(z/math.sqrt(2))+math.sqrt(2/math.pi)*z*math.exp(-z*z/2),'scope':'Triangle-inequality speed bound using global maximum endpoint bulk speed and temperature for every parent. Maxwellian nonrelativistic thermal assumption; excludes subsequent pressure acceleration or external fields.'}
  for p,n in st['hot_output'].items():counts[p]=counts.get(p,0)+w*n
 photons={};parents={};other=[]
 for p,n in counts.items():
  if p in DATA:
   d=DATA[p];line511=2*d['beta_plus_probability']*n;photons[.511]=photons.get(.511,0)+line511
   for e,v in d['lines_MeV_per_decay']:photons[e]=photons.get(e,0)+n*v
   gamma=n*(1.022*d['beta_plus_probability']+sum(e*v for e,v in d['lines_MeV_per_decay']));beta=n*(d.get('beta_mean_weighted_MeV',d.get('beta_mean_MeV',0)*d['beta_plus_probability']))
   parents[p]={'decays_per_He':n,'primary_gamma_MeV_per_He':gamma,'mean_positron_kinetic_MeV_per_He':beta,'neutral_primary_nonnu_reference_MeV_per_He':gamma+beta}
  elif p in DEFAULT_DECAYS:
   fin,q,nu=cold_decay({p:1.});other.append({'parent':p,'decays_per_He':n,'whole_atomic_Q_MeV_per_He_upper':max(q,0)*n})
 rows=[{'energy_MeV':e,'photons_per_He':n,'MeV_per_He':e*n} for e,n in sorted(photons.items())];g=sum(x['MeV_per_He'] for x in rows);beta=sum(x['mean_positron_kinetic_MeV_per_He'] for x in parents.values());total=sel['H_nonnu_MeV_per_He']+sel['U_prepared_MeV_per_He']
 return {'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'material_cycle_sha256':hashlib.sha256((ROOT/'component-design/material_cycle.py').read_bytes()).hexdigest(),'cycle_source_sha256':hashlib.sha256((ROOT/'component-design/nonideal_cycle.json').read_bytes()).hexdigest(),'mixture':HSET,'cycle':{k:sel[k] for k in ['H_nonnu_MeV_per_He','U_prepared_MeV_per_He','delayed_decay_heat_MeV_per_He','escaped_prompt_MeV_per_He','preaux_net_MeV_per_He']},'stages':stages,'parents':parents,'primary_photon_lines':rows,'primary_delayed_photon_MeV_per_He':g,'primary_delayed_photon_fraction_nonnu_nuclear':g/sel['H_nonnu_MeV_per_He'],'primary_delayed_photon_fraction_total_thermal_boundary':g/total,'main_beta_kinetic_MeV_per_He':beta,'other_unstable_whole_Q_upper_MeV_per_He':sum(x['whole_atomic_Q_MeV_per_He_upper'] for x in other),'other_unstable':other,'decay_data':DATA,'scope':'Neutral-atom decay branching, positrons assumed to slow beforeannihilation. In-flight annihilation can redistribute both positron kinetic and annihilation rest energy into a harder photon component; bounds must move both together. The saved beta kinetic energy is reported separately for that purpose. Stripped EC suppression changes dominant annihilation counts by at most~0.2%. No change to original nuclear/neutrino ledger or capture efficiency.'}
if __name__=='__main__':
 o=main();Path(__file__).with_name('spectrum.json').write_text(json.dumps(o,indent=2)+'\n');print(json.dumps({k:v for k,v in o.items() if k not in ['stages','decay_data']},indent=2))
