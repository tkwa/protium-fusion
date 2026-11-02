"""Independent moments, equilibrium convergence, optimization and allocation checks."""
from pathlib import Path
import json,gzip,itertools
from math import pi,sqrt
import numpy as np
from scipy.special import gamma,gammainc,gammaincc
from catalyst_network import load_data,HOT,label
from catalyst_catalog import nucmass
from catalyst_penning import equilibrium
from recycling_spectra import moments,weighted_moment
from network_ion_budget import geometry
from allocate_network_channels import choose
from network_plant import prepare
from nonneutral_screen import MEV,C
ROOT=Path(__file__).parent

def exact_product(tab,sub,T,z,a):
    m1=nucmass(1,1);m2=nucmass(z,a);alpha=m2/(m1+m2);mu=m1*m2/(m1+m2)*MEV/C**2
    grid=np.unique(np.r_[tab['energy_MeV'],sub['energy_MeV'],sub['yield_energy_MeV']])*alpha
    curves=[(np.asarray(xx)*alpha,np.asarray(yy)) for xx,yy in [(tab['energy_MeV'],tab['sigma_barn']),(sub['yield_energy_MeV'],sub['yield']),(sub['energy_MeV'],sub['mean_MeV'])]]
    result=0.
    for i,(lo,hi) in enumerate(zip(grid[:-1],grid[1:])):
        # Interpolation outside a table is identically zero. Avoid interpreting
        # a nonzero endpoint as a ramp over the preceding unsupported interval.
        mid=(lo+hi)/(2*alpha)
        if any(mid<min(xx) or mid>max(xx) for xx in [tab['energy_MeV'],sub['yield_energy_MeV'],sub['energy_MeV']]):continue
        poly=np.array([1.])
        for xx,yy in curves:
            # Select the native interval by its interior. Duplicate ENDF
            # energies encode a discontinuity, with distinct left/right limits.
            j=np.searchsorted(xx,(lo+hi)/2,side='right')-1
            slope=(yy[j+1]-yy[j])/(xx[j+1]-xx[j])
            poly=np.polynomial.polynomial.polymul(poly,[yy[j]-slope*xx[j],slope])
        for k,c in enumerate(poly):
            n=k+2
            delta=gammaincc(n,lo/T)-gammaincc(n,hi/T) if lo/T>5 else gammainc(n,hi/T)-gammainc(n,lo/T)
            result+=c*T**n*gamma(n)*delta
    return sqrt(8*MEV/(pi*mu))*T**-1.5*result*1e-28


def main():
    checks=[]
    def check(name,passed,metric=None):checks.append({'name':name,'passed':bool(passed),'metric':metric})
    spec=json.loads(gzip.decompress((ROOT/'recycling-spectrum-moments.json.gz').read_bytes()));data=load_data()
    energy=json.loads((ROOT/'network-ion-budget.json').read_text());fam=json.loads(gzip.decompress((ROOT/'network-channel-family.json.gz').read_bytes()))
    plants=json.loads((ROOT/'network-plant.json').read_text())['cases']+json.loads((ROOT/'network-plant-sensitivity.json').read_text())['cases']
    sub={'LAW':1,'LEP':1,'INT':np.array([2]),'yields':{'INT':np.array([2]),'Eint':np.array([1.,2.]),'yi':np.array([1.,1.])},'ZAP':0,'LIP':0,'E':np.array([1.]),'ND':[2],'Ep':[np.array([3.,1.,0.,2.])],'b':[np.array([[.2],[.3],[.25],[0.]])]}
    v=moments(sub);check('discrete lines and histogram normalization',abs(v['norm'][0]-1)<1e-15);check('discrete lines and histogram first moment',abs(v['mean_MeV'][0]-1.4e-6)<1e-18)
    for za in HOT:
        key=label(za);t=spec['targets'][key];worst=0.
        for mt,subs in t['channels'].items():
            for s in subs:
                for E,n in zip(s['energy_MeV'],s['norm']):
                    if np.interp(E,s['yield_energy_MeV'],s['yield'])>1e-10:worst=max(worst,abs(n-1))
        check(f'{key}: all populated MF6 spectra normalize',worst<2e-5,worst)
        for mt in [4,102,103]:
            for s in t['channels'][str(mt)]:
                if s['ZA'] not in [0,1]:continue
                tab=data[za]['MF3'][mt];numerical=weighted_moment(tab,s,1.6,*za);exact=exact_product(tab,s,1.6,*za)
                err=abs(numerical-exact)/max(abs(exact),1e-99)
                check(f'{key} MT{mt} product{s["ZA"]}: analytic polynomial-exponential moment',err<3e-7,float(err))
    check('thermal/direct ledger matches conserving network heat',abs(energy['energy_ledger_error_MeV'])<1e-12,energy['energy_ledger_error_MeV'])
    negative=sum(s['pn_per_primary']*min(0,x['diagnostic_heavy_CM_remainder_MeV_per_pn']) for s in energy['stages'] for x in s['events'])
    check('integrated negative CM remainder is negligible',abs(negative)<1e-8,negative)
    for plant in plants:
        a=plant['assumptions'];name=plant.get('variant',f'B{a["B"]}-h{a["heat"]}-end{a["end_in"]}')
        check(name+': installed mass ledger',abs(plant['mass_ledger_error_kg'])/plant['total_mass_kg']<1e-14,plant['mass_ledger_error_kg'])
        check(name+': energy ledger',abs(plant['energy_ledger_error_W'])/plant['nuclear_after_capture_W']<1e-14,plant['energy_ledger_error_W'])
        check(name+': nonnegative mass allocation',min(plant['mass_breakdown_kg'].values())>=0)
        vol=sum(r['channels']*pi*r['radius_m']**2*plant['active_length_m'] for r in plant['channel_allocations'])
        expected=a['packing']*pi*a['core_radius']**2*plant['active_length_m']
        check(name+': channel packing volume',abs(vol/expected-1)<1e-13,float(vol/expected-1))
        err=max(abs(r['pn_rate_s']/(plant['primary_pn_rate_s']*energy['network']['processor_pn_rates'][label(tuple(r['ZA']))])-1) for r in plant['channel_allocations'])
        check(name+': all processor flows allocated',err<1e-13,float(err))
    baseline=next(r for r in plants if r['assumptions']['B']==140 and r['assumptions']['heat']==.7 and r['assumptions']['end_in']==.999)
    e,choice=prepare(energy,fam,baseline['source_field_T']);ratio=choice['history'][-1]['ratio_MeV_m3_s']
    check('fractional allocation fixed point',abs(choice['history'][-1]['fixed_point_residual_MeV'])<1e-12)
    # The saved selected equilibria are recomputed with tighter ODE tolerance
    # and triple the radial sampling, directly at the finite source field.
    for key,r in baseline['selected_source_profiles'].items():
        s=next(s for s in energy['stages'] if label(tuple(s['ZA']))==key)
        eq=equilibrium(*r['ZA'],r['h0'],r['c0'],r['xi'],nout=12001,rtol=1e-11)
        alt=geometry(eq,r['xwall'],s,T=1.6,B=baseline['source_field_T'])
        for field,tol in [('reaction_rate_per_m_s',1e-5),('total_loss_MeV_per_pn',2e-5),('wall_E_V_m',1e-5)]:
            err=abs(alt[field]/r[field]-1);check(key+': direct field scaling and tighter '+field,err<tol,float(err))
    raw_choice=choose(energy,fam);rr=raw_choice['history'][-1]['ratio_MeV_m3_s']
    for key,selected in raw_choice['selected'].items():
        objective=selected['total_loss_MeV_per_pn']+rr*selected['area_per_pn_rate_m3_s']
        best=min(r['total_loss_MeV_per_pn']+rr*r['area_per_pn_rate_m3_s'] for r in fam['configurations'][key])
        check(key+': finite-family fractional optimality certificate',abs(objective-best)<1e-12,float(objective-best))
    toy={'bare_electricity_MeV_per_primary':20.,'network':{'processor_pn_rates':{'a':1.,'b':.25}}}
    groups={'configurations':{'a':[{'total_loss_MeV_per_pn':l,'area_per_pn_rate_m3_s':v} for l,v in [(1.,4.),(4.,2.),(10.,1.)]],'b':[{'total_loss_MeV_per_pn':l,'area_per_pn_rate_m3_s':v} for l,v in [(2.,3.),(3.,1.)]]}}
    result=choose(toy,groups)
    brute=max((20-x['total_loss_MeV_per_pn']-.25*y['total_loss_MeV_per_pn'])/(x['area_per_pn_rate_m3_s']+.25*y['area_per_pn_rate_m3_s']) for x,y in itertools.product(groups['configurations']['a'],groups['configurations']['b']))
    check('fractional optimization agrees with exhaustive toy enumeration',abs(result['history'][-1]['ratio_MeV_m3_s']-brute)<1e-13)
    ee,cc=prepare(energy,fam,140,heat=1,heater=1,direct=1,potential_in=1,potential_out=1,end_in=1,end_out=1,extraction=0,atomic=0)
    maxloss=max(abs(r['total_loss_MeV_per_pn']) for r in cc['selected'].values())
    check('unit-efficiency source limit has no processing loss',maxloss<2e-14,float(maxloss))
    out={'scope':__doc__,'passed':all(c['passed'] for c in checks),'count':len(checks),'checks':checks}
    (ROOT/'network-plant-validation.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print('passed',out['passed'],'count',len(checks))
    for c in checks:
        if not c['passed']:print(c)
    if not out['passed']:raise SystemExit(1)

if __name__=='__main__':main()
