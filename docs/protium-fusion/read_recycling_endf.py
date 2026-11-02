"""Read MF3 TAB1 tables; resonance reconstruction is deliberately not supplied."""
from pathlib import Path
import sys,hashlib,json
from read_padf import number

def read(path):
    raw=Path(path).read_bytes(); sections={}
    for s in raw.decode('latin1').splitlines():
        if len(s)<75:continue
        mf=int(s[70:72] or 0);mt=int(s[72:75] or 0)
        if mf==3 and mt>0:sections.setdefault(mt,[]).append(s)
    out={}
    for mt,lines in sections.items():
        def f(i):return [number(lines[i][j:j+11]) for j in range(0,66,11)]
        qm,qi,_,lr,nr,np=f(1);nr=int(nr);np=int(np);i=2;v=[]
        for _ in range((2*nr+5)//6):v+=f(i);i+=1
        interp=[[int(v[j]),int(v[j+1])] for j in range(0,2*nr,2)];v=[]
        for _ in range((2*np+5)//6):v+=f(i);i+=1
        assert i==len(lines),(mt,i,len(lines))
        out[mt]={'QM_MeV':qm*1e-6,'QI_MeV':qi*1e-6,'LR':int(lr),'interpolation':interp,
                 'energy_MeV':[v[j]*1e-6 for j in range(0,2*np,2)],
                 'sigma_barn':[v[j+1] for j in range(0,2*np,2)]}
    return {'source_sha256':hashlib.sha256(raw).hexdigest(),'source_file':Path(path).name,'MF3':out}

if __name__=='__main__':
    import numpy as np
    for path in sys.argv[1:]:
        d=read(path);print(Path(path).name,'MTs',list(d['MF3']))
        for mt,t in d['MF3'].items():
            if mt<200:print(mt,t['QM_MeV'],t['LR'],t['interpolation'],np.interp([2.53e-8,2,4,6,10],t['energy_MeV'],t['sigma_barn']))
