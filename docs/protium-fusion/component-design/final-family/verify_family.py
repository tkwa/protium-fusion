"""Independent per-pulse energy and partition checks of a composed family."""
from pathlib import Path
import json,math,argparse
MEV=1.602176634e-13

def verify(path):
 d=json.loads(Path(path).read_text());h=d['hot_cycle'];p=d['precursor_cycle'];checks=[]
 def check(name,value,scale=1.,tol=1e-10):
  good=abs(value)<=tol*max(abs(scale),1e-30);checks.append({'name':name,'residual':value,'passed':good});assert good,(name,value,scale)
 cats={'A':'n14','B1':'n15','B2':'c12','C13':'c13'}
 def weighted(entries,key):return sum(h['stage_pulse_data'][cats[e['record']['stage']]]['pulses_per_He']*key(e)/MEV for e in entries)
 check('prepared hot from actual pulse records',weighted(d['hot_selected'],lambda e:e['record']['initial_energy_J'])-h['U_prepared_MeV_per_He'],h['U_prepared_MeV_per_He'])
 check('prepared precursor from explicit source energy',weighted(d['plant_selected'],lambda e:e['source_prepared_energy_J'])-p['U_prepared_MeV_per_He'],p['U_prepared_MeV_per_He'])
 prep=sum(h['stage_pulse_data'][k]['pulses_per_He']*v['preparation_escaped_radiation_J']/MEV for k,v in d['preparations'].items());check('preparation radiation per pulse',prep-d['preparation_escaped_MeV_per_He'],prep)
 hotout=sum(h[k] for k in ['bulk_kinetic_MeV_per_He','hot_internal_MeV_per_He','escaped_prompt_MeV_per_He','delayed_decay_heat_MeV_per_He']);check('hot channel partition',hotout-h['H_nonnu_MeV_per_He']-h['U_prepared_MeV_per_He'],hotout)
 check('precursor channel partition',hotout+prep-p['H_nonnu_MeV_per_He']-p['U_prepared_MeV_per_He'],hotout)
 check('preparation only driver debit',p['preaux_net_MeV_per_He']-(.6*.99*(h['H_nonnu_MeV_per_He']+h['U_prepared_MeV_per_He'])-p['U_prepared_MeV_per_He']/.6))
 check('spectrum retains hot boundary',d['spectrum']['thermal_boundary_MeV_per_He']-h['H_nonnu_MeV_per_He']-h['U_prepared_MeV_per_He'])
 for k in ['stage_attempts_per_He','external_input_atoms_per_He','ash_atoms_per_He','returned_atoms_per_He','nuclear_Q_MeV_per_He','neutrino_MeV_per_He','H_nonnu_MeV_per_He','bulk_kinetic_MeV_per_He','hot_internal_MeV_per_He','escaped_prompt_MeV_per_He','delayed_decay_heat_MeV_per_He']:
  assert h[k]==p[k];checks.append({'name':'unchanged nuclear/material/output '+k,'passed':True})
 for S,z in d['capacity_plants_preparation_external_no_recovery'].items():
  rate=z['helium_per_s'];Qh=(h['H_nonnu_MeV_per_He']+h['U_prepared_MeV_per_He'])*MEV*rate;Up=p['U_prepared_MeV_per_He']*MEV*rate
  check(S+' hot capture once',z['captured_pulse_W']-.99*Qh,Qh);check(S+' driver precursor once',z['driver_W']-Up/.6,Up)
  check(S+' external escape includes prep increment once',z['escaping_source_W']-(.01*Qh+d['preparation_energy_increment_MeV_per_He']*MEV*rate),Qh)
  check(S+' total mass',sum(z['mass_breakdown_kg'].values())-z['mass_kg'],z['mass_kg']);check(S+' nuclear global power',z['nuclear_non_neutrino_W']-z['escaping_source_W']-z['net_W']-z['cold_rejection_W'],z['nuclear_non_neutrino_W'])
 return {'family':d['family'],'checks_passed':len(checks),'checks':checks}
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('integration',type=Path);a.add_argument('--out',type=Path,required=True);v=a.parse_args();d=verify(v.integration);v.out.write_text(json.dumps(d,indent=2)+'\n');print(json.dumps({'family':d['family'],'checks':d['checks_passed']}))
