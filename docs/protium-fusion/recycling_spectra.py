"""MF6 spectrum moments for the exclusive recycling reactions.

Read ENDF LAW1 spectra with separate discrete lines and continuous histograms.
LCT3 light products are in the incident pair CM frame. Heavy residual spectra
are kept only as diagnostics; their lab-frame moments are never added to CM
moments. No interpolation of an already multiplied yield/energy product.
"""
from pathlib import Path
import argparse,gzip,json,zipfile,tempfile,hashlib
import numpy as np
from endf_parserpy import EndfParserFactory
from catalyst_network import HOT,SYMBOL,EMIT,thermal_yield
from catalyst_catalog import thermal,nucmass
from nonneutral_screen import MEV,C
from math import sqrt,pi
from numpy.polynomial.legendre import leggauss
GX,GW=leggauss(16)
ROOT=Path(__file__).parent


def moments(sub):
    assert sub['LAW']==1,sub['LAW']
    assert sub['LEP'] in (1,2),sub['LEP']
    assert all(sub['INT']==2),sub['INT']
    assert all(sub['yields']['INT']==2)
    m0=[];m1=[]
    for nd,ep,b in zip(sub['ND'],sub['Ep'],sub['b']):
        p=b[:,0]; nd=int(nd)
        assert np.all(ep>=0),'Primary photon flag needs explicit treatment'
        assert np.all(np.diff(ep[nd:])>=0)
        norm=float(sum(p[:nd]));mean=float(sum(p[:nd]*ep[:nd]))
        x=ep[nd:];y=p[nd:]
        if sub['LEP']==1:
            norm+=float(sum(y[:-1]*np.diff(x)))
            mean+=float(sum(y[:-1]*np.diff(x*x)/2))
        else:
            dx=np.diff(x)
            norm+=float(sum(dx*(y[:-1]+y[1:])/2))
            mean+=float(sum(dx*(y[:-1]*(2*x[:-1]+x[1:])+y[1:]*(x[:-1]+2*x[1:]))/6))
        m0.append(norm);m1.append(mean*1e-6)
    return {'ZA':int(sub['ZAP']),'state':int(sub['LIP']),
      'energy_MeV':(sub['E']*1e-6).tolist(),'norm':m0,'mean_MeV':m1,
      'yield_energy_MeV':(sub['yields']['Eint']*1e-6).tolist(),'yield':sub['yields']['yi'].tolist(),
      'frame':'CM' if sub['ZAP']<=2004 else 'stationary_target_lab'}


def extract(archive):
    parser=EndfParserFactory.create(array_type='numpy');out={}
    with zipfile.ZipFile(archive) as zf,tempfile.TemporaryDirectory() as tmp:
        for z,a in HOT:
            name=f'p-{SYMBOL[z]}{a:03d}.tendl';path=Path(tmp)/name;raw=zf.read(name);path.write_bytes(raw)
            mts=sorted(set(int(x[72:75]) for x in raw.decode().splitlines() if len(x)>=75 and x[70:72].strip()=='6')-set([0,2,5]))
            assert set(mts)<=set(EMIT)
            sections=parser.parsefile(str(path),include=[(6,mt) for mt in mts])[6]
            result={}
            for mt in mts:
                sec=sections[mt];assert sec['LCT']==3
                result[str(mt)]=[moments(s) for s in sec['subsection']]
            out[str(z*1000+a)]={'source_sha256':hashlib.sha256(raw).hexdigest(),'channels':result}
    return {'scope':__doc__,'archive_sha256':hashlib.sha256(Path(archive).read_bytes()).hexdigest(),'targets':out}


def weighted_moment(tab,sub,T,z,a,which='mean_MeV'):
    """<sigma v yield(E) mean_energy(E)>; units MeV m3/s for mean_MeV."""
    m1=nucmass(1,1);m2=nucmass(z,a);alpha=m2/(m1+m2);mu=m1*m2/(m1+m2)*MEV/C**2
    grid=np.unique(np.r_[tab['energy_MeV'],sub['energy_MeV'],sub['yield_energy_MeV']])
    lo=grid[:-1,None];hi=grid[1:,None];half=(hi-lo)/2;e=(hi+lo)/2+half*GX
    ss=np.interp(e,tab['energy_MeV'],tab['sigma_barn'],left=0,right=0)
    yy=np.interp(e,sub['yield_energy_MeV'],sub['yield'],left=0,right=0)
    mom=np.interp(e,sub['energy_MeV'],sub[which],left=0,right=0)
    ec=alpha*e
    integ=np.sum(half*alpha*GW*ss*yy*mom*ec*np.exp(-ec/T))
    return float(sqrt(8*MEV/(pi*mu))*T**-1.5*integ*1e-28)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--archive',type=Path,default=ROOT/'beam-data/tendl2025-proton-recycling.zip');ap.add_argument('--out',type=Path,default=ROOT/'recycling-spectrum-moments.json.gz');args=ap.parse_args()
    d=extract(args.archive);args.out.write_bytes(gzip.compress(json.dumps(d,separators=(',',':'),allow_nan=False).encode(),mtime=0))
    worst=[]
    for za,t in d['targets'].items():
        for mt,subs in t['channels'].items():
            for sub in subs:
                for E,n in zip(sub['energy_MeV'],sub['norm']):
                    y=np.interp(E,sub['yield_energy_MeV'],sub['yield'])
                    if y>1e-10:worst.append((abs(n-1),za,mt,sub['ZA'],E,n))
    print('targets',len(d['targets']),'bytes',args.out.stat().st_size,'worst norms',sorted(worst,reverse=True)[:8])
