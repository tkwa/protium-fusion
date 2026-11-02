"""Current median-carbon128 photon spectrum and localized-cascade containment.

The callable evaluate(spectrum, sigma_kg_m2, ...) returns a conditional bound
for heat missed by the hot carbon receiver. It does not set mandatory capture.
No Monte Carlo or new hydrodynamic evolution is performed.
"""
from pathlib import Path
import argparse,bisect,hashlib,json,math,sys
C=299792458.;KB=1.380649e-23;MU=1.66053906660e-27;SIGMAS=[500,750,1000,1250,1500,2000,3000,5000]
CAT={'A':'n14','B1':'n15','B2':'c12','C13':'c13'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def write(p,d):Path(p).write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def coefficient(table,energy_MeV,column=2):
 if not .001<=energy_MeV<=20:raise ValueError('Coefficient requires1keV–20MeV explicit table range')
 energies=[v[0] for v in table];j=min(max(bisect.bisect_right(energies,energy_MeV)-1,0),len(table)-2);x0,y0=table[j][0],.1*table[j][column];x1,y1=table[j+1][0],.1*table[j+1][column]
 return math.exp(math.log(y0)+math.log(energy_MeV/x0)/math.log(x1/x0)*math.log(y1/y0))
def absorption_floor(table,max_energy_MeV):
 v=min(coefficient(table,max_energy_MeV),min(.1*row[2] for row in table if row[0]<=max_energy_MeV));return math.floor(v*100000)/100000

def evaluate(spectrum,sigma_kg_m2,*,mode='broad_D1p10',hot_bypass_ppm=0.,cold_intercept_fraction=0.):
 if sigma_kg_m2<0 or not 0<=cold_intercept_fraction<=1 or hot_bypass_ppm<0:raise ValueError('Invalid column or partition')
 if mode not in ['nominal_lines','broad_rest','broad_actual_bulk','broad_D1p10']:raise ValueError('Unknown spectral scenario')
 table=spectrum['NIST_carbon_table'];thermal=spectrum['thermal_boundary_MeV_per_He'];lines=spectrum['primary_photon_lines'];beta=spectrum['main_beta_kinetic_MeV_per_He'];other=spectrum['other_unstable_whole_Q_upper_MeV_per_He'];prompt=spectrum['cycle']['escaped_prompt_MeV_per_He'];D={'nominal_lines':1.,'broad_rest':1.,'broad_actual_bulk':spectrum['max_endpoint_bulk_D'],'broad_D1p10':1.10}[mode]
 atten=lambda E:math.exp(-absorption_floor(table,E)*sigma_kg_m2)
 nuclear=sum(D*l['MeV_per_He']*atten(D*l['energy_MeV']) for l in lines if l['energy_MeV']!=.511)
 anni=sum(l['MeV_per_He'] for l in lines if l['energy_MeV']==.511)
 if mode=='nominal_lines':annileak=anni*atten(.511);betaleak=beta*atten(5.14304)
 else:annileak=D*anni*atten(D*5.14304);betaleak=D*beta*atten(D*5.14304)
 hardprompt=prompt*atten(20.);otherreserve=D*other;hotreserve=hot_bypass_ppm*1e-6*thermal
 total=nuclear+annileak+betaleak+hardprompt+otherreserve+hotreserve;fraction=total/thermal
 if fraction>1:raise ValueError('Screen exceeds total available thermal boundary')
 return {'sigma_kg_m2':sigma_kg_m2,'mode':mode,'D':D,'thermal_boundary_MeV_per_He':thermal,'components_leakage_MeV_per_He':{'nuclear_lines':nuclear,'annihilation':annileak,'beta_kinetic_reserved_as_photons':betaleak,'saved_grey_prompt_at20MeV_lab':hardprompt,'other_unstable_whole_Q_full_loss':otherreserve,'carried_hot_bypass_full_loss':hotreserve},'hot_bypass_sensitivity_ppm':hot_bypass_ppm,'leakage_upper_MeV_per_He':total,'leakage_upper_fraction_thermal':fraction,'hot_capture_lower_fraction_for_included_components':1-fraction,'partition_at_this_bound':{'hot_fraction':1-fraction,'cold_equipment_fraction':cold_intercept_fraction*fraction,'external_escape_fraction':(1-cold_intercept_fraction)*fraction},'cold_intercept_fraction':cold_intercept_fraction,'scope':'Conditional energy-absorption floor with local charged-secondary deposition, no upscatter, closed-cavity returns and full carbon column after the photon source. This is a spectral scenario, not a predicted spectrum or complete apparatus capture certificate.'}

def extract(research,inputs,plant,sources):
 sys.path.insert(0,str(research/'component-design'));import material_cycle as mc
 selected=plant['top_detailed'][0];cycle=selected['cycle'];stages={};counts={};data=sources['decay_data'];hashes={};thermal=cycle['H_nonnu_MeV_per_He']+cycle['U_prepared_MeV_per_He']
 for stage,label in selected['summary']['cases'].items():
  entry=inputs['cases'][label];record=entry['record'];st=entry['stage_data'];assert record['zones']==128 and record['mixture_inputs']['rate_update']['quantile']=='median';w=cycle['stage_attempts_per_He'][CAT[stage]];vmax=max(abs(v) for v in record['profile']['velocity_m_s']);temp=max(record['profile']['T_K']);bulk=vmax/C;D=math.sqrt((1+bulk)/(1-bulk));sigma=math.sqrt(KB*temp/(13*MU))/C;bmax=(1.1**2-1)/(1.1**2+1);z=(bmax-bulk)/sigma
  tail=math.erfc(z/math.sqrt(2))+math.sqrt(2/math.pi)*z*math.exp(-z*z/2) if z>=0 else 1.
  stages[stage]={'label':label,'source':entry['path'],'record_sha256':sha(entry['path']),'attempts_per_He':w,'attempts_per_pulse':st['attempts_per_pulse'],'hot_output':st['hot_output'],'max_endpoint_speed_m_s':vmax,'max_endpoint_speed_c':bulk,'max_endpoint_bulk_D':D,'max_endpoint_temperature_K':temp,'D1p10_endpoint_thermal_tail_screen':{'minimum_parent_A':13,'one_dimension_thermal_speed_sigma_over_c':sigma,'margin_sigma':z,'Maxwell_speed_tail_probability_upper':tail,'scope':'Triangle-inequality envelope with global endpoint maximum bulk speed and temperature assigned to every parent, Maxwellian thermal speeds. Does not constrain later acceleration or an energetic non-Maxwellian population.'}}
  hashes[entry['path']]=sha(entry['path'])
  for p,n in st['hot_output'].items():counts[p]=counts.get(p,0)+w*n
 photons={};parents={};other=[]
 for p,n in counts.items():
  if p in data:
   d=data[p];photons[.511]=photons.get(.511,0)+2*d['beta_plus_probability']*n
   for E,number in d['lines_MeV_per_decay']:photons[E]=photons.get(E,0)+n*number
   gamma=n*(1.022*d['beta_plus_probability']+sum(E*v for E,v in d['lines_MeV_per_decay']));beta=n*d.get('beta_mean_weighted_MeV',d.get('beta_mean_MeV',0)*d['beta_plus_probability']);parents[p]={'decays_per_He':n,'gamma_MeV_per_He':gamma,'mean_positron_kinetic_MeV_per_He':beta,'neutral_primary_nonnu_MeV_per_He':gamma+beta}
  elif p in mc.DEFAULT_DECAYS:
   _,q,_=mc.cold_decay({p:1.});other.append({'parent':p,'decays_per_He':n,'whole_atomic_Q_MeV_per_He_upper':max(q,0)*n})
 lines=[{'energy_MeV':E,'photons_per_He':v,'MeV_per_He':E*v,'floor_mu_en_m2_kg':absorption_floor(sources['NIST_carbon_table'],E)} for E,v in sorted(photons.items())];gamma=sum(x['MeV_per_He'] for x in lines);beta=sum(x['mean_positron_kinetic_MeV_per_He'] for x in parents.values());fullybeta511=sum(1.022*(1-data[p]['beta_plus_probability'])*v['decays_per_He'] for p,v in parents.items())
 return {'schema':1,'family':'median-carbon128 A3/B1H2/B2H2/C13H4 prescribed-flow','stages':stages,'cycle':{k:cycle[k] for k in ['H_nonnu_MeV_per_He','U_prepared_MeV_per_He','delayed_decay_heat_MeV_per_He','escaped_prompt_MeV_per_He','bulk_kinetic_MeV_per_He','preaux_net_MeV_per_He']},'thermal_boundary_MeV_per_He':thermal,'parents':parents,'primary_photon_lines':lines,'primary_delayed_photon_MeV_per_He':gamma,'main_beta_kinetic_MeV_per_He':beta,'other_unstable_whole_Q_upper_MeV_per_He':sum(v['whole_atomic_Q_MeV_per_He_upper'] for v in other),'other_unstable':other,'primary_gamma_plus_beta_minus_existing_delayed_heat_MeV_per_He':gamma+beta-cycle['delayed_decay_heat_MeV_per_He'],'max_endpoint_bulk_D':max(v['max_endpoint_bulk_D'] for v in stages.values()),'fully_suppressed_EC_additional_annihilation_MeV_per_He':fullybeta511,'D1p10_extra_source_energy_drawn_from_bulk_kinetic_MeV_per_He':.1*(gamma+beta),'NIST_carbon_table':sources['NIST_carbon_table'],'decay_data':data,'source_sha256':hashes,'code_sha256':sha(__file__),'scope':'Inventory-weighted neutral-atom delayed photons and mean beta energy. D1.10 source uplift is drawn from the already counted bulk kinetic budget; it is not new nuclear energy. All grey prompt is kept as its saved energy account, with a separate20MeV laboratory spectral scenario.'}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--research-root',type=Path,required=True);ap.add_argument('--inputs',type=Path,required=True);ap.add_argument('--plant',type=Path,required=True);ap.add_argument('--sources',type=Path,required=True);ap.add_argument('--out-dir',type=Path,required=True);a=ap.parse_args();a.out_dir.mkdir(parents=True,exist_ok=True);sources=read(a.sources);s=extract(a.research_root,read(a.inputs),read(a.plant),sources);s['input_file_sha256']={str(p):sha(p) for p in [a.inputs,a.plant,a.sources,a.research_root/'component-design/material_cycle.py']};write(a.out_dir/'spectrum.json',s)
 rows=[{'sigma_kg_m2':S,'scenarios':{mode:evaluate(s,S,mode=mode) for mode in ['nominal_lines','broad_rest','broad_actual_bulk','broad_D1p10']},'broad_D1p10_plus_carried6p4ppm':evaluate(s,S,hot_bypass_ppm=6.4)} for S in SIGMAS]
 result={'spectrum_sha256':sha(a.out_dir/'spectrum.json'),'code_sha256':sha(__file__),'columns_kg_m2':SIGMAS,'rows':rows,'carried_hot_bypass_scope':'Optional6.4ppm of current H+U as full unrecovered energy; inherited conservative sensitivity from old selected64 positive-flow screen, not recomputed for current128 or updated rates; may double-count saved grey prompt.','geometry_scope':'Any escaping photon must cross the stated credited hot carbon column. Parent/material routing behind that shield or interception by cold hardware needs explicit partition. No MC reused as a current128 result.'};write(a.out_dir/'containment.json',result)
 print(json.dumps({'boundary_MeV_He':s['thermal_boundary_MeV_per_He'],'parents':s['parents'],'max_bulk_D':s['max_endpoint_bulk_D'],'minimum_thermal_margin_sigma':min(v['D1p10_endpoint_thermal_tail_screen']['margin_sigma'] for v in s['stages'].values()),'rows':[{'sigma':x['sigma_kg_m2'],'nominal_pct':100*x['scenarios']['nominal_lines']['leakage_upper_fraction_thermal'],'broad_pct':100*x['scenarios']['broad_D1p10']['leakage_upper_fraction_thermal'],'broad_pluscarry_pct':100*x['broad_D1p10_plus_carried6p4ppm']['leakage_upper_fraction_thermal']} for x in rows]},indent=2))
if __name__=='__main__':main()
