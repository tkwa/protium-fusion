"""Independent arithmetic, envelope, coefficient and energy-partition checks."""
from pathlib import Path
import json,math,hashlib,importlib.util
P=Path(__file__).resolve().parent
sp=importlib.util.spec_from_file_location('current_photon128',P/'current_photon128.py');m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)
s=json.loads((P/'spectrum.json').read_text());d=json.loads((P/'containment.json').read_text());checks=[]
def ck(name,truth):
 checks.append({'name':name,'passed':bool(truth)});assert truth,name
ck('spectrumcodehash',s['code_sha256']==m.sha(P/'current_photon128.py'));ck('containmentcodehash',d['code_sha256']==m.sha(P/'current_photon128.py'));ck('spectruminputhash',d['spectrum_sha256']==m.sha(P/'spectrum.json'))
for p,h in s['source_sha256'].items():ck('recordhash '+p,m.sha(p)==h)
for p,h in s['input_file_sha256'].items():ck('inputhash '+p,m.sha(p)==h)
counts={}
for v in s['stages'].values():
 for p,n in v['hot_output'].items():counts[p]=counts.get(p,0)+v['attempts_per_He']*n
for p,v in s['parents'].items():ck('decayweight '+p,math.isclose(counts[p],v['decays_per_He'],rel_tol=1e-14))
ck('annihilationphotoncount',math.isclose(next(v['photons_per_He'] for v in s['primary_photon_lines'] if v['energy_MeV']==.511),2*sum(counts[p]*v['beta_plus_probability'] for p,v in s['decay_data'].items()),rel_tol=1e-14))
ck('primaryphotonenergy',math.isclose(sum(v['MeV_per_He'] for v in s['primary_photon_lines']),s['primary_delayed_photon_MeV_per_He'],rel_tol=1e-14))
ck('betakineticenergy',math.isclose(sum(v['mean_positron_kinetic_MeV_per_He'] for v in s['parents'].values()),s['main_beta_kinetic_MeV_per_He'],rel_tol=1e-14))
ck('DenergyavailablefromK',s['D1p10_extra_source_energy_drawn_from_bulk_kinetic_MeV_per_He']<s['cycle']['bulk_kinetic_MeV_per_He'])
for stage,v in s['stages'].items():
 b=v['max_endpoint_speed_c'];ck('Ddefinition '+stage,math.isclose(math.sqrt((1+b)/(1-b)),v['max_endpoint_bulk_D'],rel_tol=1e-14));ck('D1.10coversbulk '+stage,v['max_endpoint_bulk_D']<1.10);ck('thermalmargin '+stage,v['D1p10_endpoint_thermal_tail_screen']['margin_sigma']>9)
for row in s['NIST_carbon_table']:ck('NISTknot '+str(row[0]),math.isclose(m.coefficient(s['NIST_carbon_table'],row[0]),.1*row[2],rel_tol=1e-12))
for E in [.511,1.6352,2.312593,3.9475,5.14304,1.1*5.14304,20.]:
 k=m.absorption_floor(s['NIST_carbon_table'],E);v=[m.coefficient(s['NIST_carbon_table'],math.exp(math.log(.001)+(math.log(E)-math.log(.001))*i/1000)) for i in range(1001)];ck('degradedenergyfloor '+str(E),k<=min(v)*(1+1e-12))
previous={}
for row in d['rows']:
 S=row['sigma_kg_m2'];modes=['nominal_lines','broad_rest','broad_actual_bulk','broad_D1p10'];last=0
 for mode in modes:
  v=row['scenarios'][mode];result=m.evaluate(s,S,mode=mode);ck(f'callableexact {S} {mode}',v==result);loss=v['leakage_upper_fraction_thermal'];ck(f'spectralordering {S} {mode}',loss>=last);last=loss;ck(f'columnordering {S} {mode}',mode not in previous or loss<previous[mode]);previous[mode]=loss
  ck(f'componentclosure {S} {mode}',math.isclose(sum(v['components_leakage_MeV_per_He'].values())/s['thermal_boundary_MeV_per_He'],loss,rel_tol=1e-14))
 for f in [0,.5,1]:
  p=m.evaluate(s,S,cold_intercept_fraction=f)['partition_at_this_bound'];ck(f'coldescapepartition {S} {f}',math.isclose(sum(p.values()),1,abs_tol=1e-14));ck(f'coldescapesplit {S} {f}',math.isclose(p['cold_equipment_fraction'],f*last,abs_tol=1e-14))
 carry=row['broad_D1p10_plus_carried6p4ppm'];ck(f'carry6p4ppm {S}',math.isclose(carry['leakage_upper_fraction_thermal']-last,6.4e-6,abs_tol=1e-14))
out={'passed':sum(x['passed'] for x in checks),'total':len(checks),'checks':checks,'scope':__doc__,'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};(P/'validation.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'passed':out['passed'],'total':out['total']}))
