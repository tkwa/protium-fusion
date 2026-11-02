from pathlib import Path
import json,hashlib,sys,re,argparse
import numpy as np
ap=argparse.ArgumentParser();ap.add_argument('--text',type=Path,required=True);ap.add_argument('--pdf',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
sys.path.insert(0,'/home/tkwa/.codex/worktrees/protium-fusion/website/docs/protium-fusion');import pulse_model as pm
s=a.text.read_text();start=s.index('Table 35:');section=s[start:s.index('REV =',start)];rows=[]
for line in section.replace('−','-').splitlines():
 try:v=list(map(float,line.split()))
 except ValueError:continue
 if len(v)==8:rows.extend([v[:4],v[4:]])
rows.sort();assert len(rows)==56,len(rows)
r=next(r for r in pm.RATES if r.fname=='N15_p_to_He4_C12_reaclib')
data=[]
for t,med,lo,hi in rows:
 old=r.eval(t*1e9);data.append(dict(T9=t,adopted=med,low=lo,high=hi,saved_REACLIB=old,adopted_over_saved=med/old,low_over_saved=lo/old,high_over_saved=hi/old))
a.out.write_text(json.dumps({'source':'Xu et al. 2013, NACRE II, Table35, PDFpage96','doi':'10.1016/j.nuclphysa.2013.09.007','pdf_url':'https://arxiv.org/pdf/1310.7099','pdf_sha256':hashlib.sha256(a.pdf.read_bytes()).hexdigest(),'confidence':'Adopted/low/high model evaluation limits; not assigned a Gaussian confidence level here.','rate_units':'cm3 mol-1 s-1','table_text':section,'rows':data},indent=2)+'\n')
for d in data:
 if d['T9'] in [.2,.25,.3,.35,.4,.45,.5,.7,1.,1.5,2.]:print({k:round(v,6) for k,v in d.items() if k in ['T9','adopted_over_saved','low_over_saved','high_over_saved']})
