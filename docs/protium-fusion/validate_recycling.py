"""Independent integral and conservation checks for the recycling audit."""
from pathlib import Path
import json,zipfile,tempfile
from math import sqrt,pi
import numpy as np
from scipy.special import gammainc,gammaincc,gamma
from endf_parserpy import EndfParserFactory
from catalyst_network import load_data,audit,cold_routes,HOT,EMIT,ARCHIVE
from catalyst_catalog import thermal,nucmass,atomicmass
from nonneutral_screen import MEV,C
ROOT=Path(__file__).parent

def analytic_thermal(tab,T,z,a,power=0):
    alpha=nucmass(z,a)/(nucmass(z,a)+nucmass(1,1));m1=nucmass(1,1);m2=nucmass(z,a)
    mu=m1*m2/(m1+m2)*MEV/C**2
    E=np.array(tab['energy_MeV'])*alpha;s=np.array(tab['sigma_barn']);total=0.
    def moment(lo,hi,n):
        if lo/T>5:delta=gammaincc(n+1,lo/T)-gammaincc(n+1,hi/T)
        else:delta=gammainc(n+1,hi/T)-gammainc(n+1,lo/T)
        return T**(n+1)*gamma(n+1)*delta
    for lo,hi,sl,sh in zip(E[:-1],E[1:],s[:-1],s[1:]):
        if hi==lo:continue
        b=(sh-sl)/(hi-lo);a0=sl-b*lo
        total+=a0*moment(lo,hi,1+power)+b*moment(lo,hi,2+power)
    return sqrt(8*MEV/(pi*mu))*T**-1.5*total*1e-28

def main():
    data=load_data();checks=[]
    def check(name,passed,metric=None):
        checks.append({'name':name,'passed':bool(passed),'metric':metric})
    parser=EndfParserFactory.create(array_type='numpy')
    with zipfile.ZipFile(ARCHIVE) as zf,tempfile.TemporaryDirectory() as temp:
        for za,d in data.items():
            f=Path(temp)/d['source_file'];f.write_bytes(zf.read(f.name))
            alt=parser.parsefile(str(f),include=[3])[3]
            equal=all(np.array_equal(np.array(t['energy_MeV']),alt[mt]['xstable']['E']*1e-6)
                      and np.array_equal(t['sigma_barn'],alt[mt]['xstable']['xs']) for mt,t in d['MF3'].items())
            check(f'{za}: independent MF3 reader',equal)
            for mt in [4,102,107,111]:
                if mt not in d['MF3']:continue
                for power in [0,1]:
                    t=d['MF3'][mt];v=thermal(t,1.6,*za,power);a=analytic_thermal(t,1.6,*za,power)
                    err=abs(v-a)/max(abs(a),1e-100)
                    check(f'{za} MT{mt} moment{power}: incomplete-gamma integral',err<2e-8,float(err))
    for T in [1.,1.3,1.6,1.9,2.2]:
        d=audit(data,T);b=abs(d['baryon_ledger_error_per_primary_pn']);q=abs(d['net_atomic_mass_energy_ledger_error_MeV'])
        check(f'T{T}: baryon conservation including feed and waste',b<2e-12,b)
        check(f'T{T}: independent atomic-mass energy ledger',q<5e-9,q)
        check(f'T{T}: stationary auxiliary isotope balance',d['matrix_nonprimary_balance_max']<1e-12,d['matrix_nonprimary_balance_max'])
        stock=d['Ti49_feed_for_boundary_per_primary_pn']-d['boundary_residual_sum']-d['high_energy_discard_events_per_primary_pn']
        check(f'T{T}: fresh catalyst equals all heavy losses',abs(stock)<2e-12,float(stock))
        check(f'T{T}: nonnegative processor rates',min(d['processor_pn_rates'].values())>=-1e-14)
    routes=cold_routes((25,54))
    check('Mn54 weak probabilities sum to one',abs(sum(r['fraction'] for r in routes)-1)<1e-15)
    check('Mn53 promptly enters proton processor',cold_routes((25,53))[0]['target']==(25,53))
    check('V50 radiative-capture path reaches V51',cold_routes((23,50))[0]['target']==(23,51))
    out={'scope':__doc__,'passed':all(c['passed'] for c in checks),'checks':checks,'count':len(checks)}
    (ROOT/'recycling-validation.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print('passed',out['passed'],'checks',len(checks))
    for c in checks:
        if not c['passed']:print(c)
    if not out['passed']:raise SystemExit(1)

if __name__=='__main__':main()
