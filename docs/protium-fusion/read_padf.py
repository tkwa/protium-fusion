"""Read TAB1 production cross sections in PADF-2 MF10/MT5 (ENDF-6)."""
import re, json, hashlib
from pathlib import Path

def number(s):
 s=s.strip()
 if not s: return 0.
 if 'e' not in s.lower(): s=re.sub(r'(?<=\d)([+-]\d+)$',r'e\1',s)
 return float(s)

def read(path):
 raw=Path(path).read_bytes()
 lines=[s for s in raw.decode('latin1').splitlines() if len(s)>=75 and s[70:72].strip()=='10' and s[72:75].strip()=='5']
 def fields(i): return [number(lines[i][j:j+11]) for j in range(0,66,11)]
 count=int(fields(0)[4]); i=1; out={}
 for k in range(count):
  qm,qi,za,state,nr,np=fields(i); nr,np,za,state=map(int,[nr,np,za,state]); i+=1
  v=[]
  for _ in range((nr*2+5)//6): v+=fields(i); i+=1
  interp=[[int(v[j]),int(v[j+1])] for j in range(0,2*nr,2)]
  v=[]
  for _ in range((np*2+5)//6): v+=fields(i); i+=1
  out[f'{za}_{state}']={'product_ZA':za,'state':state,'QM_raw':qm,'QI_raw':qi,'interpolation':interp,'energy_MeV':[v[j]*1e-6 for j in range(0,2*np,2)],'sigma_barn':[v[j+1] for j in range(0,2*np,2)]}
 assert i==len(lines),(i,len(lines))
 return {'source_sha256':hashlib.sha256(raw).hexdigest(),'evaluation':'PADF-2 C13 EVAL-MAR26, distributed 2026-04-23','source_url':'https://bwsyncandshare.kit.edu/s/YwCwKH5qSjBPjfx?dir=/Data/Z__6_41','units_note':'Cross sections barn, incident laboratory energy MeV. Q raw fields appear in MeV despite ENDF convention; do not use them in energy accounting.','channels':out}

if __name__=='__main__':
 import numpy as np
 root=Path(__file__).parent
 d=read(root/'beam-data'/'padf2-c13-2026.ndf'); (root/'padf2-c13.json').write_text(json.dumps(d,indent=2)+'\n')
 for k in ['1_0','1002_0','5010_0','6012_0','6013_0','7013_0','7014_0']:
  v=d['channels'][k]; print(k,'interp',v['interpolation'],'sigma at4,5,6,7,8,10MeV',np.interp([4,5,6,7,8,10],v['energy_MeV'],v['sigma_barn']).round(6))
 print('parsed',len(d['channels']),'channels')
