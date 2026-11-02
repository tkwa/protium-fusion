"""Cold-return isotope-flow audit with TENDL exclusive proton channels.

This is not a complete plant model. Capture selection, extraction and converter
design remain unproved. Boundary products and evaluated-data inconsistencies
are reported, never silently returned to Ti49.
"""
from pathlib import Path
import sys,json,gzip,re,zipfile,tempfile,hashlib,argparse
from collections import defaultdict
from functools import lru_cache
from math import sqrt,pi,log
import numpy as np
from numpy.polynomial.legendre import leggauss
from endf_parserpy import EndfParserFactory
from read_recycling_endf import read
from catalyst_catalog import thermal,nucmass,atomicmass,nubase_records,nuclear
from nonneutral_screen import MEV,C,Q,DD_CHARGED_PER_NEW_N,DD_THERMAL_PER_NEW_N
ROOT=Path(__file__).parent
OUT=ROOT
ARCHIVE=ROOT/'beam-data'/'tendl2025-proton-recycling.zip'
HOT=[(22,49),(23,51),(24,52),(24,53),(24,54),(25,53),(25,55),(26,56),(26,57),(26,58),(27,59),(28,59),(28,60),(28,61),(28,62),(28,64)]
SYMBOL={22:'Ti',23:'V',24:'Cr',25:'Mn',26:'Fe',27:'Co',28:'Ni'}
NB=nubase_records();GX,GW=leggauss(16)
# Emitted light species, defined by ENDF-6 MT numbers. MT103 is inelastic
# proton scattering; incident-proton elastic scattering uses MT2.
EMIT={4:{(0,1):1},16:{(0,1):2},17:{(0,1):3},22:{(0,1):1,(2,4):1},
      28:{(0,1):1,(1,1):1},37:{(0,1):4},102:{},103:{(1,1):1},
      104:{(1,2):1},105:{(1,3):1},106:{(2,3):1},107:{(2,4):1}}
EMIT.update({11:{(0,1):2,(1,2):1},24:{(0,1):2,(2,4):1},25:{(0,1):3,(2,4):1},
 29:{(0,1):1,(2,4):2},30:{(0,1):2,(2,4):2},32:{(0,1):1,(1,2):1},
 33:{(0,1):1,(1,3):1},34:{(0,1):1,(2,3):1},35:{(0,1):1,(1,2):1,(2,4):2},
 36:{(0,1):1,(1,3):1,(2,4):2},41:{(0,1):2,(1,1):1},42:{(0,1):3,(1,1):1},
 44:{(0,1):1,(1,1):2},45:{(0,1):1,(1,1):1,(2,4):1},108:{(2,4):2},
 109:{(2,4):3},111:{(1,1):2},112:{(1,1):1,(2,4):1},113:{(1,3):1,(2,4):2},
 114:{(1,2):1,(2,4):2},115:{(1,1):1,(1,2):1},116:{(1,1):1,(1,3):1},117:{(1,2):1,(2,4):1}})

def label(za):
    if za is None:return 'unknown'
    return str(za[0]*1000+za[1])

@lru_cache(maxsize=None)
def weak_isomer_margin(za):
    """Upper excitation margin if every produced atom occupies a weak isomer.

    NUBASE states 8/9 are excluded (isobaric analog / special records).
    IT-only isomers return their excitation electromagnetically.
    """
    values=[]
    for (z,a,state),r in NB.items():
        if (z,a)!=za or not 0<state<8:continue
        if not any(s in r['decay_modes'] for s in ['B+','B-','EC']):continue
        match=re.search(r'[0-9]+(?:\.[0-9]*)?',r['record'][42:54])
        if match:values.append(float(match.group())*1e-3)
    return max(values,default=0.)

def load_data(archive=ARCHIVE):
    p=EndfParserFactory.create(array_type='numpy');out={}
    with zipfile.ZipFile(archive) as zf,tempfile.TemporaryDirectory(prefix='protium-tendl-') as temp:
        manifest=json.loads(zf.read('download-manifest.json'))
        expected={r['name']:r['sha256'] for r in manifest}
        for za in HOT:
            z,a=za;name=f'p-{SYMBOL[z]}{a:03d}.tendl';raw=zf.read(name)
            assert hashlib.sha256(raw).hexdigest()==expected[name]
            path=Path(temp)/name;path.write_bytes(raw);d=read(path)
            d6=p.parsefile(str(path),include=[(6,5)])[6][5]
            yields=[]
            for sub in d6['subsection']:
                if sub['ZAP']==0:continue
                y=sub['yields'];assert all(y['INT']==2)
                yields.append({'ZA':int(sub['ZAP']),'state':int(sub['LIP']),
                               'energy_MeV':(y['Eint']*1e-6).tolist(),'yield':y['yi'].tolist()})
            d['MT5_yields']=yields;out[za]=d
    return out


def thermal_yield(tab,y,T,z,a,power=0):
    """Integrate the product of separately interpolated sigma and yield."""
    alpha=nucmass(z,a)/(nucmass(z,a)+nucmass(1,1))
    mu=nucmass(z,a)*nucmass(1,1)/(nucmass(z,a)+nucmass(1,1))*MEV/C**2
    grid=np.unique(np.r_[tab['energy_MeV'],y['energy_MeV']])
    lo=grid[:-1,None];hi=grid[1:,None];half=(hi-lo)/2
    e=(hi+lo)/2+half*GX
    xs=np.interp(e,tab['energy_MeV'],tab['sigma_barn'],left=0,right=0)
    yy=np.interp(e,y['energy_MeV'],y['yield'],left=0,right=0)
    cm=e*alpha
    integral=np.sum(half*alpha*GW*xs*yy*cm**(1+power)*np.exp(-cm/T))
    return float(sqrt(8*MEV/(pi*mu))*T**-1.5*integral*1e-28)

def cold_path(za):
    """Prescribed ideal isotope-separated neutron captures / cold decays."""
    cur=za;steps=[];seen=set()
    for _ in range(40):
        if cur in HOT:return cur,steps,None
        if cur in seen:return None,steps,'cold cycle'
        seen.add(cur);z,a=cur;r=NB.get((z,a,0))
        if not r:return None,steps,'missing decay data'
        mode=r['decay_modes']
        if r['half_life_raw']!='stbl' and cur!=(23,50):
            # A genuinely mixed return is retained as a boundary for now.
            if ('B-' in mode) and ('B+' in mode or 'EC' in mode):
                return None,steps,'mixed weak branches: '+mode
            if mode.startswith(('B+=100','EC=100')):nxt=(z-1,a);kind='weak'
            elif mode.startswith('B-=100'):nxt=(z+1,a);kind='weak'
            else:return None,steps,'unhandled decay: '+mode
        elif (z<=22 and a<=50) or cur in [(23,50),(24,50),(26,54),(28,58)]:
            nxt=(z,a+1);kind='capture'
        else:return None,steps,'stable outside selected processor set'
        q=atomicmass(*cur)+(atomicmass(0,1) if kind=='capture' else 0)-atomicmass(*nxt)
        assert q>=0,(cur,nxt,kind,q)
        steps.append({'from':cur,'to':nxt,'kind':kind,'Q_MeV':q})
        cur=nxt
    return None,steps,'path too long'

def cold_routes(za,depth=0):
    assert depth<10
    target,steps,err=cold_path(za)
    terminal=tuple(steps[-1]['to']) if steps else za
    # NUBASE percentages: the Mn54 beta-minus branch is 0.93e-4 percent.
    # Its positron branch and EC both return Cr54. Preserve the weak branch
    # to Fe54 rather than treating the entire Mn54 flow as a poison.
    if err and terminal==(25,54):
        out=[]
        for nextza,frac in [((24,54),1-9.3e-7),((26,54),9.3e-7)]:
            q=atomicmass(*terminal)-atomicmass(*nextza)
            for r in cold_routes(nextza,depth+1):
                out.append({**r,'fraction':frac*r['fraction'],
                    'steps':steps+[{'from':terminal,'to':nextza,'kind':'weak','Q_MeV':q}]+r['steps']})
        return out
    return [{'target':target,'steps':steps,'boundary':err,'terminal':terminal,'fraction':1.}]

def stage(za,d,T):
    z,a=za;tabs=d['MF3'];sv=thermal(tabs[4],T,z,a);events=[]
    assert set(tabs)-{2,3,5} <= set(EMIT)
    assert all(law==2 for t in tabs.values() for _,law in t['interpolation'])
    for mt,light in EMIT.items():
        if mt not in tabs:continue
        tab=tabs[mt];rate=thermal(tab,T,z,a)/sv
        if not rate:continue
        ecm=thermal(tab,T,z,a,1)/(sv*rate)
        residual=(z+1-sum(k[0]*v for k,v in light.items()),a+1-sum(k[1]*v for k,v in light.items()))
        routes=cold_routes(residual)
        qp=atomicmass(z,a)+atomicmass(1,1)-atomicmass(*residual)-sum(atomicmass(*k)*v for k,v in light.items())
        events.append({'mt':mt,'rate':rate,'Ecm_MeV':ecm,'residual':residual,
                       'light':light,'routes':routes,'Q_prompt_MeV':qp})
    # Independent mean product yields are sufficient for linear isotope
    # flow. Their mass/charge/one-heavy-residual sum rules are audited below.
    mt5=thermal(tabs[5],T,z,a)/sv
    prods=[]
    for y in d['MT5_yields']:
        rate=thermal_yield(tabs[5],y,T,z,a)/sv
        if not rate:continue
        res=divmod(y['ZA'],1000)
        if res[1]>4:
            routes=cold_routes(res)
        else:routes=[]
        prods.append({'rate':rate,'residual':res,'state':y['state'],'routes':routes})
    return {'ZA':za,'sv':sv,'events':events,'mt5_rate':mt5,'mt5_products':prods,
            'exclusive_to_nonelastic_ratio':sum(thermal(t,T,z,a) for mt,t in tabs.items() if mt not in [2,3])/thermal(tabs[3],T,z,a),
            'mt5_heavy_per_event':sum(x['rate'] for x in prods if x['residual'][1]>4)/mt5 if mt5 else 0,
            'mt5_A_error_per_pn':sum(x['residual'][1]*x['rate'] for x in prods)-(a+1)*mt5,
            'mt5_Z_error_per_pn':sum(x['residual'][0]*x['rate'] for x in prods)-(z+1)*mt5}

def audit(data,T=1.6,mt5_policy='discard'):
    stages=[stage(za,data[za],T) for za in HOT];ix={za:i for i,za in enumerate(HOT)}
    matrix=np.zeros((len(HOT),len(HOT)))
    for j,s in enumerate(stages):
        for e in s['events']:
            matrix[j,j]-=e['rate']
            for r in e['routes']:
                if r['target'] in ix:matrix[ix[r['target']],j]+=e['rate']*r['fraction']
        matrix[j,j]-=s['mt5_rate']
        if mt5_policy=='evaluated_yields':
            for e in s['mt5_products']:
                for r in e['routes']:
                    if r['target'] in ix:matrix[ix[r['target']],j]+=e['rate']*r['fraction']
    exposure=np.r_[1.,np.linalg.solve(matrix[1:,1:],-matrix[1:,0])]
    assert min(exposure)>-1e-12
    balance=matrix@exposure
    light=defaultdict(float);captures=defaultdict(float);decays=defaultdict(float);boundaries=defaultdict(float)
    nu=0.;capture_q=0.;prompt_q=0.;mt5_discard=0.;mt5_waste_mass=0.;mt5_waste_A=0.;isomer_margin=0.
    def process(e,w):
        nonlocal nu,capture_q,isomer_margin
        isomer_margin+=w*weak_isomer_margin(tuple(e['residual']))
        for r in e['routes']:
            wr=w*r['fraction']
            for st in r['steps']:
                isomer_margin+=wr*weak_isomer_margin(tuple(st['to']))
                if st['kind']=='capture':captures[tuple(st['from'])]+=wr;capture_q+=wr*st['Q_MeV'];light[(0,1)]-=wr
                else:decays[tuple(st['from'])]+=wr;nu+=wr*st['Q_MeV']
            if r['boundary']:boundaries[(tuple(r['terminal']),r['boundary'])]+=wr
    for s,f in zip(stages,exposure):
        for e in s['events']:
            w=f*e['rate'];light[(1,1)]-=w
            for p,v in e['light'].items():light[p]+=w*v
            prompt_q+=w*e['Q_prompt_MeV'];process(e,w)
        light[(1,1)]-=f*s['mt5_rate']
        if mt5_policy=='discard':
            mt5_discard+=f*s['mt5_rate']
            mt5_waste_mass+=f*s['mt5_rate']*(atomicmass(*s['ZA'])+atomicmass(1,1))
            mt5_waste_A+=f*s['mt5_rate']*(s['ZA'][1]+1)
        else:
            prompt_q+=f*s['mt5_rate']*(atomicmass(*s['ZA'])+atomicmass(1,1))
            for e in s['mt5_products']:
                w=f*e['rate'];r=e['residual'];prompt_q-=w*atomicmass(*r)
                if r[1]>4:process(e,w)
                else:light[r]+=w
    # Hypothetical downstream light-fuel completion, separate from boundaries:
    # n+p -> D; net 2D -> alpha; T+p -> alpha through DT plus recapture;
    # He3+n -> alpha. No external D/T/He3 fuel is supplied.
    n=light[(0,1)];d=light[(1,2)];t=light[(1,3)];he3=light[(2,3)];he4=light[(2,4)]
    assert n>=he3
    alpha=he4+t+he3+(n-he3+d)/2
    fuel=-light[(1,1)]+n-he3+t
    qlight=(n-he3)*Q['np_capture']+(n-he3+d)*.5*Q['net_D']
    qlight+=t*(atomicmass(1,3)+atomicmass(1,1)-atomicmass(2,4))
    qlight+=he3*(atomicmass(2,3)+atomicmass(0,1)-atomicmass(2,4))
    # Ground-state atomic weak Q plus a worst-population margin for known
    # weakly decaying isomers is assigned to neutrinos. The actual branching
    # and gamma/charged energy partition remain unresolved.
    ground_nu=nu;nu+=isomer_margin
    nuclear_heat=prompt_q+capture_q+qlight-isomer_margin
    feed=-balance[0]
    waste_mass=mt5_waste_mass+sum(atomicmass(*z)*v for (z,err),v in boundaries.items())
    waste_A=mt5_waste_A+sum(z[1]*v for (z,err),v in boundaries.items())
    expected_heat=fuel*atomicmass(1,1)+feed*atomicmass(22,49)-alpha*atomicmass(2,4)-waste_mass-nu
    pure_H_heat=alpha*Q['net_H']-nu
    inventories=[]
    for za,r in decays.items():
        half=nuclear(*za).tau
        if half is not None:
            inventories.append({'nuclide':label(za),'decays_per_primary_pn':r,'mean_life_s':float(half/log(2)),
                                'kg_per_primary_pn_s':float(r*half/log(2)*atomicmass(*za)*MEV/C**2)})
    return {'T_MeV':T,'mt5_policy':mt5_policy,'processor_pn_rates':{label(z):float(f) for z,f in zip(HOT,exposure)},
      'matrix_nonprimary_balance_max':float(max(abs(balance[1:]))),'Ti49_feed_for_boundary_per_primary_pn':float(-balance[0]),
      'light_net_before_D_completion':{label(z):v for z,v in light.items()},
      'capture_steps_per_primary_pn':{label(z):v for z,v in captures.items()},
      'cold_weak_steps_per_primary_pn':{label(z):v for z,v in decays.items()},
      'boundary_products':[{'nuclide':label(z),'reason':err,'rate':v} for (z,err),v in sorted(boundaries.items(),key=lambda x:-x[1])],
      'boundary_residual_sum':sum(boundaries.values()),'all_neutrons_captured_on_catalysts':sum(captures.values()),
      'high_energy_discard_events_per_primary_pn':mt5_discard,
      'discard_electricity_reserve_1000MeV_per_event_MeV':1000*mt5_discard,
      'net_atomic_mass_energy_ledger_error_MeV':nuclear_heat-expected_heat,
      'baryon_ledger_error_per_primary_pn':fuel+49*feed-4*alpha-waste_A,
      'trace_makeup_transmutation_heat_MeV':nuclear_heat-pure_H_heat,
      'heat_after_excluding_positive_makeup_transmutation_MeV':min(nuclear_heat,pure_H_heat),
      'net_neutrons_available_for_light_fuels':n,'net_alpha_per_primary_pn':alpha,'net_H_per_primary_pn':fuel,
      'H_minus_4alpha_per_primary_pn':fuel-4*alpha,'neutrino_upper_MeV':nu,
      'ground_state_neutrino_bound_MeV':ground_nu,'known_weak_isomer_excitation_margin_MeV':isomer_margin,
      'prompt_reaction_Q_MeV':prompt_q,'cold_capture_Q_MeV':capture_q,'downstream_light_Q_MeV':qlight,
      'nuclear_heat_excluding_boundary_restoration_MeV':nuclear_heat,'decay_inventories':sorted(inventories,key=lambda x:-x['kg_per_primary_pn_s']),
      'data_checks':[{'target':label(s['ZA']),'pn_sv':s['sv'],'exclusive_to_nonelastic_ratio':s['exclusive_to_nonelastic_ratio'],
                     'mt5_rate':s['mt5_rate'],'mt5_heavy_per_event':s['mt5_heavy_per_event'],
                     'mt5_A_error_per_pn':s['mt5_A_error_per_pn'],'mt5_Z_error_per_pn':s['mt5_Z_error_per_pn']} for s in stages]}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--archive',type=Path,default=ARCHIVE);args=ap.parse_args()
    data=load_data(args.archive);out={'scope':__doc__,'archive_sha256':hashlib.sha256(args.archive.read_bytes()).hexdigest(),'cases':[audit(data,T) for T in [1.,1.3,1.6,1.9,2.2]],
                       'raw_yield_diagnostic_1p6':audit(data,1.6,mt5_policy='evaluated_yields')}
    (OUT/'catalyst-network-audit.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    for d in out['cases']:
        print('T',d['T_MeV'],'processors',d['processor_pn_rates'])
        print('capture',d['all_neutrons_captured_on_catalysts'],'netn',d['net_neutrons_available_for_light_fuels'],'H/4He error',d['H_minus_4alpha_per_primary_pn'])
        print('boundary',d['boundary_residual_sum'],d['boundary_products'][:8])
        print('inventories',d['decay_inventories'][:3])
