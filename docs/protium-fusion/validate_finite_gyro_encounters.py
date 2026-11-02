"""Independent real-space versus spectral finite-gyroradius checks."""
from pathlib import Path
import argparse,json
from math import pi
import numpy as np
from scipy.integrate import quad
from scipy.special import spherical_jn
from finite_gyro_encounters import radial_field


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--outdir',type=Path,default=Path(__file__).parent);args=ap.parse_args();root=args.outdir
    checks=[]
    def check(name,ok,metric=None):checks.append({'name':name,'passed':bool(ok),'metric':metric})
    for r in [0.,1e-6,.01,.0199,.0201,.1,1.,10.]:
        V,h,hp=radial_field(r)
        def bessel_ratio(k,order):
            x=k*r
            if x<1e-3:
                return (1/3-x*x/30+x**4/840) if order==1 else (1/15-x*x/210+x**4/7560)
            return spherical_jn(order,x)/x**order
        Vq=-2/pi*quad(lambda k:k*k/(k*k+1)*np.exp(-k*k/32)*np.sinc(k*r/pi),0,40,epsabs=1e-12)[0]
        hq=2/pi*quad(lambda k:k**4/(k*k+1)*np.exp(-k*k/32)*bessel_ratio(k,1),0,40,epsabs=1e-12)[0]
        hpq=-2/pi*quad(lambda k:k**6/(k*k+1)*np.exp(-k*k/32)*bessel_ratio(k,2),0,40,epsabs=1e-12)[0]
        for name,a,b in [('potential',V,Vq),('radial_force',h,hq),('force_gradient',hp,hpq)]:
            err=float(abs(a-b)/max(abs(b),1e-10))
            check(f'Independent inverse Fourier {name} r={r}',err<2e-8,err)
    rows=[]
    for eta,omega in [(0.,-1.),(.1,1.)]:
        path=root/f'finite-gyro-encounters-96-24-24-eta{eta:g}-omega{omega:g}.json'
        real=json.loads(path.read_text())
        name=f'finite-gyro-beam-spectrum-eta{eta:g}.json'
        spectrum=json.loads((root/name).read_text());w=np.array(real['w']);b=np.r_[w,1.]
        comparisons={}
        for key in ['diffusion','drift']:
            got=np.array(real[key]);want=np.array(spectrum[key]);err=float(np.linalg.norm(got-want)/np.linalg.norm(want))
            check(f'Real-space versus spectral {key}, eta={eta}, omega={omega}',err<1e-5,err);comparisons[key]=err
        for order in ['first','second']:
            err=real[f'{order}_order_energy_error']
            check(f'Individual orbit {order}-order energy identity, eta={eta}',err<1e-7,err)
        D=np.array(spectrum['diffusion']);A=np.array(spectrum['drift'])
        err=float(np.linalg.norm(D@b)/np.linalg.norm(D))
        check(f'Spectral covariance work identity, eta={eta}',err<1e-12,err)
        residual=b@A+.5*(D[2,2]+eta*np.trace(D[:3,:3]))
        err=float(abs(residual)/max(np.linalg.norm(A),1.))
        check(f'Spectral mean work identity including recoil, eta={eta}',err<2e-7,err)
        check(f'Encounter solver terminal success, eta={eta}',real['solver_success'])
        rows.append({'encounter_file':path.name,'spectral_file':name,'relative_errors':comparisons})
    reference=json.loads((root/'finite-gyro-beam-spectrum-eta0.json').read_text())
    errors=[]
    for name in ['finite-gyro-encounters-40-8-8-eta0-omega-1.json','finite-gyro-encounters-64-16-16-eta0-omega-1.json','finite-gyro-encounters-96-24-24-eta0-omega-1.json']:
        result=json.loads((root/name).read_text())
        errors.append([float(np.linalg.norm(np.array(result[key])-reference[key])/np.linalg.norm(reference[key])) for key in ['diffusion','drift']])
    check('Impact and gyroangle refinement reduces both errors',bool(np.all(np.diff(np.array(errors),axis=0)<0)),errors)
    out={'scope':__doc__,'passed':sum(c['passed'] for c in checks),'total':len(checks),'checks':checks,'comparisons':rows}
    (root/'finite-gyro-encounter-validation.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(out['passed'],'/',out['total'],'checks pass')
    for c in checks:
        if not c['passed']:print('FAILED',c)
    if out['passed']!=out['total']:raise SystemExit(1)


if __name__=='__main__':main()
