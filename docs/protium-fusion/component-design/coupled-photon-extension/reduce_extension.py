from pathlib import Path
import copy,json,hashlib,gzip,statistics,math,sys
P=Path(__file__).resolve().parent;R=Path('/home/tkwa/.codex/worktrees/protium-fusion/website/docs/protium-fusion');B=R/'component-design/coupled-photon-transport';sys.path.insert(0,str(B));from reweight import evaluate
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,d):p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
base=json.loads((B/'basis.json').read_text());basis={'schema':1,'columns':{},'source_sha256':base['source_sha256'],'scope':base['scope'],'base_basis_reference':{'path':str(B/'basis.json'),'sha256':sha(B/'basis.json')},'raw_inputs':{}};(P/'events').mkdir(exist_ok=True);checks=[]
for S in [200,350]:
 node={'sigma_kg_m2':S,'responses':[]};basis['columns'][str(S)]=node
 for kind in ['broad','prompt20']:
  name=f'S{S}-{kind}';p=P/(name+'.json');d=json.loads(p.read_text());assert d['code_sha256']==base['source_sha256'];ev={'parameters':d['parameters'],'columns':['deposited_MeV','escaped_calorimetric_MeV','vacuum_deposit_MeV','residual_inside_kinetic_MeV'],'records':[]}
  for r in d['records']:
   vals=[[e[k] for k in ev['columns']] for e in r['events']];E=r['incident_calorimetric_energy_MeV'];a=[v[1]/E for v in vals];N=len(a);mean=math.fsum(a)/N;se=math.sqrt(math.fsum((v-mean)**2 for v in a)/(N*(N-1)));errs=[math.fsum(v)/E-1 for v in vals]
   for n,v in [('N',N==r['histories']),('mean',abs(mean-r['leakage_fraction'])<1e-13),('SE',abs(se-r['sampling_standard_error'])<1e-13),('energy',max(abs(x) for x in errs)<1e-7),('finished',all(v[3]==0 for v in vals))]:checks.append({'name':f'{name}:{E}:{n}','passed':v});assert v
   ev['records'].append({'particle':'gamma','kinetic_energy_MeV':E,'calorimetric_energy_MeV':E,'values':vals});node['responses'].append({'particle':'gamma','mode':d['parameters']['mode'],'energy_MeV':E,'calorimetric_energy_MeV':E,'histories':N,'leakage_fraction':mean,'deposited_fraction':r['deposited_fraction'],'sampling_standard_error':se,'sampling_80pct_half_width':1.2815515655446004*se,'max_relative_event_energy_residual':max(abs(x) for x in errs),'mean_absolute_relative_energy_residual':statistics.mean(abs(x) for x in errs),'event_archive':f'events/{name}.json.gz','event_archive_record':len(ev['records'])-1,'raw_name':name})
  ep=P/'events'/f'{name}.json.gz'
  with ep.open('wb') as f:
   with gzip.GzipFile(filename='',mode='wb',fileobj=f,mtime=0) as g:g.write(json.dumps(ev,separators=(',',':')).encode())
  basis['raw_inputs'][name]={'path':str(p),'sha256':sha(p),'log_sha256':sha(P/(name+'.log')),'source_sha256':d['code_sha256'],'event_archive_sha256':sha(ep)}
write(P/'basis-extension.json',basis)
combined=copy.deepcopy(base)
for v in combined['columns'].values():
 for q in v['responses']:q['event_archive']='../coupled-photon-transport/'+q['event_archive']
combined['columns'].update(copy.deepcopy(basis['columns']));combined['raw_inputs'].update(basis['raw_inputs']);combined['extension_provenance']={'base_basis_sha256':sha(B/'basis.json'),'extension_basis_sha256':sha(P/'basis-extension.json'),'scope':'Original columns unchanged; original event archive paths adjusted to the neighboring frozen packet. New columns have only broad delayed and20MeV prompt kernels.'};write(P/'combined-basis.json',combined)
sp=Path('/home/tkwa/src/protium-fusion/output/final-family/delivered-B1-C13H2/spectrum.json');s=json.loads(sp.read_text());write(P/'final-family-spectrum.json',s);rows=[]
for S in [200,250,350,500,750,1000,2000]:
 a=evaluate(combined,s,S,hot_bypass_ppm=6.4);b=evaluate(combined,s,S,prompt_policy='all_to_cold',hot_bypass_ppm=6.4);assert a['sampling_80pct_half_width']<.001;rows.append({'sigma_kg_m2':S,'broad_EM20':a,'broad_all_prompt_to_cold':b})
write(P/'weighted-results.json',{'family':s['family'],'spectrum_source_path':str(sp),'spectrum_sha256':sha(sp),'rows':rows});write(P/'validation.json',{'checks_passed':len(checks),'checks':checks,'new_histories':sum(r['histories'] for c in basis['columns'].values() for r in c['responses']),'max_relative_event_energy_error':max(r['max_relative_event_energy_residual'] for c in basis['columns'].values() for r in c['responses'])})
print(json.dumps([{'sigma':r['sigma_kg_m2'],'loss':r['broad_EM20']['leakage_fraction'],'fullprompt':r['broad_all_prompt_to_cold']['leakage_fraction'],'half80':r['broad_EM20']['sampling_80pct_half_width']} for r in rows],indent=2))
