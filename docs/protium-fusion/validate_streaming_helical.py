"""Regression of the helical operator's 28 mathematical checks.

Runs the streamed recurrence/moment backend through limiting operators,
finite-difference divergence, Maxwell bath polarization, the cold boundary,
positive covariance and independent quadrature refinement. Source-state
resolution and matching are checked separately in validate_helical_resolution.py.
"""
from pathlib import Path
import argparse,json,time
import numpy as np
from finite_gyro_spectrum import unmagnetized_gyro_average
from streaming_helical_spectrum import helical_coefficients
from guiding_center_collisions import coefficients


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--outdir',type=Path,default=Path(__file__).parent);args=ap.parse_args()
    start=time.monotonic();checks=[];cases=[]
    def check(name,ok,metric=None):checks.append({'name':name,'passed':bool(ok),'metric':metric})
    def evaluate(vp=1.,vz=.8,U=-.7,omega=1.,**kw):
        config={'nk':80,'nmu':48,'nphi':24};config.update(kw)
        c=helical_coefficients(vp,vz,U,1.,omega,.001,1.,4.,**config)
        cases.append({'vperp':vp,'vz':vz,'U':U,'Omega':omega,'quadrature':kw,
                      **{k:(v.tolist() if isinstance(v,np.ndarray) else v) for k,v in c.items()}})
        return c
    for vp in [0.,1.]:
        c=evaluate(vp=vp,omega=1e4)
        raw=coefficients(.8,-.7,np.array([[c['spectral_log'],.001,1.]]))
        expectedD=np.array([[raw[3],raw[5],0.],[raw[5],raw[4],0.],[0.,0.,0.]])
        expectedA=np.array([raw[0],raw[1],0.])
        for key,want in [('diffusion',expectedD),('drift',expectedA)]:
            err=float(np.linalg.norm(c[key]-want)/np.linalg.norm(want))
            check(f'Strong-field distant limit {key}, vperp={vp}',err<1e-7,err)
    base=evaluate()
    h=1e-4
    plus=evaluate(vz=.8+h);minus=evaluate(vz=.8-h)
    ep=evaluate(vp=np.sqrt(1+2*h));em=evaluate(vp=np.sqrt(1-2*h))
    fd=(plus['diffusion'][:,1]-minus['diffusion'][:,1])/(2*h)+(ep['diffusion'][:,2]-em['diffusion'][:,2])/(2*h)
    err=float(np.linalg.norm(fd-base['divergence'])/np.linalg.norm(base['divergence']))
    check('Divergence matches independent state finite differences',err<5e-8,err)
    pol=-.001*base['diffusion']@np.array([.7,.8,1.])
    err=float(np.linalg.norm(pol-base['polarization'])/np.linalg.norm(pol))
    check('Maxwell bath polarization balance',err<1e-12,err)
    fine=evaluate(nk=112,nmu=64,nphi=64)
    for key in ['drift','diffusion']:
        err=float(np.linalg.norm(base[key]-fine[key])/np.linalg.norm(fine[key]))
        check(f'Wavevector quadrature refinement {key}',err<1e-7,err)
    cold=evaluate(vp=0.)
    almost=evaluate(vp=1e-3)
    check('Cold perpendicular boundary has zero energy diffusion',abs(cold['diffusion'][2,2])<1e-20)
    check('Cold perpendicular boundary heats inward',cold['drift'][2]>0,float(cold['drift'][2]))
    err=float(np.linalg.norm(cold['drift']-almost['drift'])/np.linalg.norm(cold['drift']))
    check('Continuous cold-boundary drift',err<3e-6,err)
    errors=[]
    for omega in [.1,.05]:
        c=evaluate(omega=omega)
        zero=unmagnetized_gyro_average(1.,.8,-.7,1.,.001,c['spectral_log'])
        row=[]
        for key in ['drift','diffusion']:
            err=float(np.linalg.norm(c[key]-zero[key])/np.linalg.norm(zero[key]))
            row.append(err);check(f'Approach straight-orbit limit {key}, Omega={omega}',err<.003,err)
        errors.append(row)
    ratio=np.array(errors[1])/errors[0]
    check('Weak-field errors shrink on halving gyrofrequency',bool(np.all(ratio<.4)),ratio.tolist())
    for i,c in enumerate(cases):
        D=np.array(c['diffusion']);ev=np.linalg.eigvalsh(D)
        check(f'Positive covariance case {i}',ev.min()>-1e-12*ev.max(),float(ev.min()/ev.max()))
    out={'scope':__doc__,'passed':sum(c['passed'] for c in checks),'total':len(checks),'checks':checks,'cases':cases,'seconds':time.monotonic()-start}
    (args.outdir/'streaming-helical-validation.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(out['passed'],'/',out['total'],'checks pass in',out['seconds'],'s')
    for c in checks:
        if not c['passed']:print('FAILED',c)
    if out['passed']!=out['total']:raise SystemExit(1)


if __name__=='__main__':main()
