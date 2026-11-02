"""Independent current-rate128 identities, ranking, replay and rate-boundary checks."""
import argparse,pathlib,json,hashlib,importlib.util,math,sys
ap=argparse.ArgumentParser();ap.add_argument('--research-root',required=True,type=pathlib.Path);ap.add_argument('--inputs',required=True,type=pathlib.Path);ap.add_argument('--results-dir',required=True,type=pathlib.Path);ap.add_argument('--output',required=True,type=pathlib.Path);a=ap.parse_args();r=a.research_root
sp=importlib.util.spec_from_file_location('verify128plant',r/'mixture-design/plant_compare.py');pc=importlib.util.module_from_spec(sp);sp.loader.exec_module(pc)
x=json.loads(a.inputs.read_text());cases=x['cases'];checks=[]
def ck(name,error,limit=1e-10):
 ok=bool(abs(error)<=limit);checks.append({'name':name,'absolute_error':float(abs(error)),'limit':limit,'passed':ok});assert ok,(name,error)
for sigma in [5000,2000]:
 d=json.loads((a.results_dir/f'current-rate128-sigma{sigma}.json').read_text());rank=d['ranked'];ck(str(sigma)+' all36combinations',len(rank)-36,0);ck(str(sigma)+' no rejection',len(d['rejected']),0);ck(str(sigma)+' sorted score',max([0]+[rank[i+1]['mW_kg']-rank[i]['mW_kg'] for i in range(len(rank)-1)]),0)
 for j,row in enumerate(d['top_detailed']):
  p=row['plant'];c=row['cycle'];s=row['summary'];selected=[cases[k] for k in s['cases'].values()];cy=pc.cycle_for(selected);replay=pc.best_rate(selected,cy,absorber_sigma=sigma)
  ck(f'{sigma} row{j} score',p['net_W']/p['mass_kg']*1000/s['mW_kg']-1)
  ck(f'{sigma} row{j} mass',sum(p['mass_breakdown_kg'].values())/p['mass_kg']-1)
  ck(f'{sigma} row{j} receiver mass',p['mass_breakdown_kg']['receiver']/(4*math.pi*p['radius_m']**2*sigma)-1)
  electrical=p['gross_electric_W']-p['driver_W']-p['ordinary_auxiliary_W']-p['pump_W']-p['replacement_W']-p['support_W'];ck(f'{sigma} row{j} net electrical',electrical/p['net_W']-1)
  nuclear=p['nuclear_non_neutrino_W'];ck(f'{sigma} row{j} globalenergy',(nuclear-p['escaping_source_W']-p['net_W']-p['cold_rejection_W'])/nuclear)
  ck(f'{sigma} row{j} fresh carbon H',p['fresh_C_kg_s']/p['fresh_H_kg_s']/c['fresh_C_H_mass_ratio']-1)
  for k,v in c['stage_pulse_data'].items():ck(f'{sigma} row{j} cadence{k}',p['pulse_interval_s'][k]*v['pulses_per_He']*p['helium_per_s']-1)
  occ=sum(c['stage_pulse_data'][pc.CAT[sx['record']['stage']]]['pulses_per_He']*sx['receiver']['timing_at_radius']['mass999']['occupied_s'] for sx in selected)*p['radius_m']/1e9*p['helium_per_s'];ck(f'{sigma} row{j} occupancy',occ/p['sequential_occupancy_mass999']-1)
  ck(f'{sigma} row{j} material replay',cy['U_prepared_MeV_per_He']/c['U_prepared_MeV_per_He']-1)
  ck(f'{sigma} row{j} plant replay',replay['net_W_kg']/p['net_W_kg']-1)
 best=d['top_detailed'][0];selected=[cases[k] for k in best['summary']['cases'].values()];p=best['plant'];cy=best['cycle'];lo=pc.station(selected,cy,p['helium_per_s']*(1-1e-5),absorber_sigma=sigma);hi=pc.station(selected,cy,p['helium_per_s']*(1+1e-5),absorber_sigma=sigma)
 ck(str(sigma)+' radiator boundary crossed',0 if lo['radiator_area_fraction']<.9<hi['radiator_area_fraction'] else 1,0)
 ck(str(sigma)+' occupancy below criterion',max(0,p['sequential_occupancy_mass999']-.8),0)
 ck(str(sigma)+' C/H below availablebudget',max(0,p['fresh_C_H_mass_ratio']-.007),0)
 ck(str(sigma)+' actual-composition64notmixed',0 if all(v['record']['zones']==128 and v['record']['mixture_inputs']['rate_update']['quantile']=='median' for v in selected) else 1,0)
out={'scope':__doc__,'passed':sum(c['passed'] for c in checks),'total':len(checks),'checks':checks,'code_sha256':hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()};a.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'passed':out['passed'],'total':out['total']}))
