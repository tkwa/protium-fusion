"""Receiver mass versus explicitly routed missed heat, on one frozen pulse family.

99% hot capture is a comparison choice, not a mandatory design requirement.
Uses a separate generalized copy of the original plant ledger. All missing hot
heat is either allowed to escape all equipment or charged to the cold system.
"""
from pathlib import Path
import argparse,hashlib,importlib.util,json,math,sys
import numpy as np
from scipy.optimize import brentq
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import plant_compare as pc
import station_generalized as generalized

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def load_module(name,path):
 spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def optimized(selected,cycle,sigma,capture,hot_fraction):
 radius=1.1e9
 occ=sum(cycle['stage_pulse_data'][pc.CAT[s['record']['stage']]]['pulses_per_He']*s['receiver']['timing_at_radius']['mass999']['occupied_s'] for s in selected)*radius/pc.RREF
 peak=max(s['receiver']['source_rating']['peak_total_prepared_electric_W']*pc.RREF/radius for s in selected)
 ratecap=min(.8/occ,peak*.6/(cycle['U_prepared_MeV_per_He']*pc.cold.MEV))
 opts={'energy_capture_fraction':capture,'hot_eligible_fraction':hot_fraction}
 def at(f):return pc.station(selected,cycle,ratecap*f,radius,absorber_sigma=sigma,cold_parameters=opts)
 ans=at(1.)
 if ans.get('radiator_area_fraction',0)>.9:
  lo=None
  for f in [.5,.25,.125,.0625,.03125,.015625]:
   z=at(f)
   if 'radiator_area_fraction' in z and z['radiator_area_fraction']<.9:lo=f;break
  if lo is None:raise ValueError('No physical rate bracket')
  f=brentq(lambda x:at(x)['radiator_area_fraction']-.9,lo,1.,xtol=1e-13)*(1-1e-10)
  ans=at(f)
 return ans

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--inputs',type=Path,required=True);ap.add_argument('--plant',type=Path,required=True);ap.add_argument('--photon-code',type=Path,required=True);ap.add_argument('--spectrum',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--columns',default='500,750,1000,1250,1500,2000,3000,5000');args=ap.parse_args()
 inputs=read(args.inputs);base=read(args.plant)['top_detailed'][0];selected=[inputs['cases'][label] for label in base['summary']['cases'].values()];cycle=pc.cycle_for(selected)
 assert all(s['record']['zones']==128 for s in selected)
 assert abs(cycle['U_prepared_MeV_per_He']/base['cycle']['U_prepared_MeV_per_He']-1)<1e-12
 photon=load_module('current_photons',args.photon_code);spectrum=read(args.spectrum)
 assert abs(spectrum['thermal_boundary_MeV_per_He']/(cycle['H_nonnu_MeV_per_He']+cycle['U_prepared_MeV_per_He'])-1)<1e-12
 original=pc.best_rate(selected,cycle,absorber_sigma=5000.)
 pc.cold=generalized
 default=optimized(selected,cycle,5000.,.99,1.)
 checkfields=['net_W','mass_kg','helium_per_s','cold_rejection_W','radiator_area_fraction','driver_W','support_W','replacement_W']
 checks={k:bool(abs(default[k]/original[k]-1)<1e-11) for k in checkfields};assert all(checks.values())
 rows=[]
 for sigma in [float(v) for v in args.columns.split(',')]:
  bound=photon.evaluate(spectrum,sigma,mode='broad_D1p10',hot_bypass_ppm=6.4)
  for extra in [0.,.005,.01]:
   loss=bound['leakage_upper_fraction_thermal']+extra;assert 0<=loss<1
   for routing in ['external_escape','cold_equipment']:
    capture,hot=(1-loss,1.) if routing=='external_escape' else (1.,1-loss)
    p=optimized(selected,cycle,sigma,capture,hot)
    assert p['feasible'],p.get('reasons')
    assert max(abs(v) for v in p['closure_relative_residuals'].values())<1e-10
    rows.append({'receiver_sigma_kg_m2':sigma,'photon_bound_fraction_HplusU':bound['leakage_upper_fraction_thermal'],'additional_missed_hot_fraction_HplusU':extra,'missed_hot_routing':routing,'total_hot_fraction':1-loss,'mW_kg':1000*p['net_W_kg'],'net_W':p['net_W'],'mass_kg':p['mass_kg'],'He_s':p['helium_per_s'],'occupancy':p['sequential_occupancy_mass999'],'radiator_fraction':p['radiator_area_fraction'],'heat_in_cold_equipment_W':p['cold_heat_components_W']['captured_source_ineligible_for_hot_engine'],'escape_W':p['escaping_source_W'],'plant':p})
 output={'scope':__doc__,'family':spectrum['family'],'selected_cases':base['summary']['cases'],'cycle':cycle,'default_reproduction_checks':checks,'default_reproduction':default,'rows':rows,'input_sha256':{str(p):sha(p) for p in [args.inputs,args.plant,args.photon_code,args.spectrum,Path(pc.__file__),Path(generalized.__file__)]},'code_sha256':sha(__file__),'interpretation':['All rows are conditional accounts, subject to the photon-bound assumptions, finite-depth surface heating, backing heat capacity, guide support and actual component qualification.','External escape and complete deposition in cold equipment are limiting routings. Intermediate cases require a specified hot/cold geometry.','The0,0.5%,1% extra missed-hot fractions are sensitivity allocations, not calculated losses or requirements. The separate1W/m² hot-to-cold bypass remains charged in every row.','A smaller receiver column does not create a loss of intercepted ordinary matter. Skin retention and processing remain a separate material account.','Cold radiator, source peak, service mass and capacity equations are rebuilt consistently. Source hardware efficiencies and installed ratings retain their explicit current assumptions.']}
 args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
 print(json.dumps({'checks':checks,'rows':[{k:v for k,v in r.items() if k!='plant'} for r in rows if r['additional_missed_hot_fraction_HplusU']==.005]},indent=2))
if __name__=='__main__':main()
