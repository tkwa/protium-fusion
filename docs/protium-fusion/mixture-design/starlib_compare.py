"""Compare selected archived rates with the official STARLIB2022 table."""
import json,zipfile,pathlib,sys,hashlib,argparse
import numpy as np
sys.path.insert(0,'/home/tkwa/.codex/worktrees/protium-fusion/website/docs/protium-fusion');import pulse_model as pm
here=pathlib.Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--archive',type=pathlib.Path,required=True);ap.add_argument('--output',type=pathlib.Path,default=here/'starlib-key-rate-comparison.json');a=ap.parse_args()
wanted={('p','n13','o14'):'N13_p_to_O14_reaclib',('p','n14','o15'):'N14_p_to_O15_reaclib',('p','n15','he4','c12'):'N15_p_to_He4_C12_reaclib'}
with zipfile.ZipFile(a.archive) as f:lines=f.read('starlib.dat').decode().splitlines()
rows=[]
for i,l in enumerate(lines):
 if not l.startswith(('4 ','5 ')):continue
 nuclei=tuple(l[5:40].split())
 if nuclei not in wanted:continue
 name=wanted[nuclei];tab=np.array([list(map(float,s.split())) for s in lines[i+1:i+61]])
 r=next(r for r in pm.RATES if r.fname==name);data=[]
 for t,v,factor in tab:
  if t<.1:continue
  old=r.eval(t*1e9);data.append({'T9':t,'STARLIB':v,'factor_uncertainty':factor,'REACLIB':old,'ratio':v/old})
 rows.append({'rate':name,'header':l,'rows':data,'raw_record':'\n'.join(lines[i:i+61])+'\n'})
assert len(rows)==3,len(rows)
a.output.write_text(json.dumps({'official_archive_sha256':hashlib.sha256((a.archive).read_bytes()).hexdigest(),'rates':rows},indent=2)+'\n')
for r in rows:
 print(r['rate'],r['header']);print([(v['T9'],v['ratio'],v['factor_uncertainty']) for v in r['rows'] if v['T9'] in [.3,.5,1.,1.5,2.]])
