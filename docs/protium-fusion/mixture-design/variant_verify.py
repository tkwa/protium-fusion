"""Explicit variant selection, baseline and replay regression checks."""
import argparse,json,importlib.util,tempfile,copy,hashlib
from pathlib import Path
ap=argparse.ArgumentParser();ap.add_argument('--utility',type=Path,required=True);ap.add_argument('--research-root',type=Path,required=True);ap.add_argument('--reference-result',type=Path,required=True);ap.add_argument('--source-result',type=Path,required=True);ap.add_argument('--replay-result',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
spec=importlib.util.spec_from_file_location('variant_under_test',a.utility);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
checks=[]
def ck(n,b):
 checks.append({'name':n,'passed':bool(b)})
 if not b:raise AssertionError(n)
ref=m.read(a.reference_result);src=m.read(a.source_result);rep=m.read(a.replay_result)
ck('reference exact capacity reproduction',all(v==0 for v in ref['reference_comparison']['capacity_relative_changes'].values()))
ck('reference exact cycle reproduction',all(v==0 for v in ref['reference_comparison']['cycle_relative_changes'].values()))
ck('source/replay cycle physical fields identical',all(v==rep['cycle'][k] for k,v in src['cycle'].items() if k!='mass_defect_Q_residual_MeV_per_He'))
ck('source/replay subtraction residual remains negligible',max(abs(x['cycle']['mass_defect_Q_residual_MeV_per_He']) for x in [src,rep])<1e-10)
ck('source/replay full plant equality',src['capacity_plant']==rep['capacity_plant'])
ck('all integration checks pass',all(c['passed'] for d in [ref,src,rep] for c in d['checks']))
root=a.research_root.resolve();pro=root/'component-design/selected-source-stage-map.json';mapping=m.read(root/'component-design/selected-source-paths.json')
class NeverReduce:
 def reduce_case(self,*args,**kw):raise AssertionError('Invalid source was allowed into receiver')
with tempfile.TemporaryDirectory(prefix='protium-variant-guards-') as td:
 td=Path(td)
 cases=[]
 bad=copy.deepcopy(mapping);bad.pop('c13');cases.append(('missing stage rejected',bad))
 bad=copy.deepcopy(mapping);bad['n14']=mapping['c12'];cases.append(('wrong stage label rejected',bad))
 bad=copy.deepcopy(mapping);bad['n14']={'record':mapping['n14'],'result_sha256':'0'*64};cases.append(('mismatched source hash rejected',bad))
 bad=copy.deepcopy(mapping);bad['n14']={'record':mapping['n14'],'eos':'NonidealEOS'};cases.append(('unadapted EOS rejected',bad))
 bad=copy.deepcopy(mapping);bad['n14']={'record':mapping['n14'],'duration_s':1.};cases.append(('unfinished physical duration rejected',bad))
 for title,case in cases:
  path=td/'map.json';m.write(path,case)
  try:m.extract_explicit(path,root,td/'work',NeverReduce(),None,provenance_path=pro)
  except ValueError:ck(title,True)
  else:ck(title,False)
out={'scope':__doc__,'checks':checks,'passed':sum(c['passed'] for c in checks),'source_sha256':{str(p):m.sha(p) for p in [a.utility,a.reference_result,a.source_result,a.replay_result]},'code_sha256':m.sha(__file__)};m.write(a.output,out);print(json.dumps({'passed':out['passed'],'total':len(checks)},indent=2))
