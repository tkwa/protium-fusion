"""Independent numerical checks of resolved source quadrature and Bessel streaming.

This validates evaluation of the specified weak static-screened operator. It
makes no assertion about a stationary population, physical screening, curved
ion orbits, electron supply, collection or reactor power.
"""
from pathlib import Path
import argparse,json,time
import numpy as np
from scipy.special import jv
from fast_bessel_grid import integer_bessel_grid,integer_bessel_pairs
from finite_gyro_spectrum import helical_coefficients as original
from moment_helical_spectrum import helical_coefficients as moments
from aligned_helical_spectrum import helical_coefficients as aligned
from streaming_helical_spectrum import helical_coefficients as streamed


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--data',type=Path,default=Path(__file__).parent)
    args=ap.parse_args();start=time.monotonic();checks=[];comparisons=[]
    def check(name,value,tol):
        checks.append({'name':name,'passed':bool(value<tol),'metric':float(value),'upper_limit':tol})
    x=np.unique(np.r_[0.,np.geomspace(1e-10,1000,500),np.arange(1,501)-.01,np.arange(1,501),np.arange(1,501)+.01])
    nmax=1100;n=np.arange(nmax+1)
    table=integer_bessel_grid(x,nmax);reference=jv(n[:,None],x[None,:])
    check('Stable Bessel recurrence against independent scipy.jv',np.max(np.abs(table-reference)),1e-12)
    check('Bessel squared completeness over 0<=x<=1000',np.max(np.abs(table[0]**2+2*np.sum(table[1:]**2,axis=0)-1)),1e-12)
    check('Bessel weighted squared identity',np.max(np.abs(2*np.sum(n[1:,None]**2*table[1:]**2,axis=0)-x*x/2)/np.maximum(1,x*x/2)),1e-12)
    streamed_error=0.
    for ni,current,following in integer_bessel_pairs(x,nmax-1):
        streamed_error=max(streamed_error,float(np.max(np.abs(current-table[ni]))),float(np.max(np.abs(following-table[ni+1]))))
    check('Streaming and full Bessel tables agree',streamed_error,1e-14)
    test_states=[(1.,.8,-.7,1.),(0.,.8,-.7,1.),(1.,0.,0.,1.),(1.,.8,-.7,1e4),(.6,-1.7,.9,.7)]
    for i,(vp,vz,U,Omega) in enumerate(test_states):
        parameters=(vp,vz,U,1.,Omega,.001,1.,4.)
        target=original(*parameters,nk=80,nmu=48,nphi=48)
        for name,function,grid,tolerance in [('moments',moments,(80,48,24),1e-12),('streamed',streamed,(80,48,24),1e-12),('aligned',aligned,(96,64,48),1e-7)]:
            got=function(*parameters,nk=grid[0],nmu=grid[1],nphi=grid[2])
            for key in ['drift','diffusion','divergence','polarization']:
                err=np.linalg.norm(got[key]-target[key])/max(1e-15,np.linalg.norm(target[key]))
                check(f'{name} versus original, state {i}, {key}',err,tolerance)
    def read(name):return json.loads((args.data/name).read_text())
    def compare(label,left,right,tol):
        a=read(left);b=read(right);row={'name':label,'left':left,'right':right,'relative_changes':{}}
        for key in ['perpendicular_heating_eV_s','parallel_heating_eV_s','total_heating_eV_s','electron_radial_drift_m_s']:
            err=abs(a['matched'][key]-b['matched'][key])/abs(b['matched'][key]);row['relative_changes'][key]=err
            if tol is not None:check(label+': '+key,err,tol)
        v=b['species'][1]['vth_m_s'];S=np.diag([1/v,1/v,1/v**2])
        for key in ['drift','diffusion']:
            aa=np.array(a['matched'][key]);bb=np.array(b['matched'][key])
            if key=='drift':aa=S@aa;bb=S@bb
            else:aa=S@aa@S;bb=S@bb@S
            err=float(np.linalg.norm(aa-bb)/np.linalg.norm(bb));row['relative_changes']['normalized_'+key]=err
            if tol is not None:check(label+': normalized '+key,err,tol)
        comparisons.append(row)
    compare('Independent aligned versus magnetic-axis integration','matched-aligned-warm-cut1-128-64-64.json','matched-moment-warm-cut1-384-192-32.json',1e-7)
    compare('Polar refinement at factor 1','matched-moment-warm-cut1-256-64-32.json','matched-moment-warm-cut1-256-192-32.json',1e-12)
    compare('Radial refinement at factor 1','matched-moment-warm-cut1-256-192-32.json','matched-moment-warm-cut1-384-192-32.json',2e-5)
    compare('Radial refinement at factor 2','matched-moment-warm-cut2-512-192-32.json','matched-moment-warm-cut2-768-192-32.json',1e-11)
    compare('Streaming versus full special-function table at factor 2','matched-moment-warm-cut2-512-192-32.json','matched-streaming-warm-cut2-512-192-32.json',1e-12)
    compare('Independent recurrence and radial refinement at factor 4','matched-moment-warm-cut4-1024-256-32.json','matched-streaming-warm-cut4-1536-256-32.json',1e-11)
    compare('Three-grid refinement at factor 8','matched-streaming-warm-cut8-2048-384-32.json','matched-streaming-warm-cut8-3072-512-48.json',1e-8)
    for c,lower,upper in [(1,'matched-moment-warm-cut1-384-192-32.json','matched-moment-warm-cut2-768-192-32.json'),(2,'matched-moment-warm-cut2-768-192-32.json','matched-streaming-warm-cut4-1536-256-32.json'),(4,'matched-streaming-warm-cut4-1536-256-32.json','matched-streaming-warm-cut8-3072-512-48.json')]:
        compare(f'Filter scale {c} to {2*c}',lower,upper,None)
    for name in sorted(set(c[k] for c in comparisons for k in ['left','right'])):
        d=read(name);v=d['species'][1]['vth_m_s'];S=np.diag([1/v,1/v,1/v**2]);ev=np.linalg.eigvalsh(S@np.array(d['matched']['diffusion'])@S)
        check('Positive covariance '+name,-float(ev.min()/ev.max()),1e-13)
        check('Positive complementary logarithm '+name,-min(s['short_logarithm'] for s in d['species']),0.)
    out={'scope':__doc__,'passed':sum(c['passed'] for c in checks),'total':len(checks),'checks':checks,'source_comparisons':comparisons,'seconds':time.monotonic()-start}
    (args.data/'helical-resolution-validation.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:out[k] for k in ['passed','total','seconds']}),flush=True)
    for c in checks:
        if not c['passed']:print('FAILED',c,flush=True)
    if out['passed']!=out['total']:raise SystemExit(1)


if __name__=='__main__':main()
