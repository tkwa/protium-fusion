"""Photon-associated internal-pair diagnostic for the recycled catalyst network.

The ENDF spectra do not identify each transition multipolarity or emitting
level. Each pure E1/M1/E2/M2/E3/M3 assignment is a conditional scenario, not an
actual branching reconstruction or an uncertainty bound. E0, atomic cascade
changes on stripping, and missing/high-energy photon contributions are separate.
"""
from pathlib import Path
import argparse, os, subprocess, json, gzip, hashlib, zipfile, tempfile
import xml.etree.ElementTree as ET
import numpy as np
from numpy.polynomial.legendre import leggauss
from endf_parserpy import EndfParserFactory
from catalyst_network import HOT,SYMBOL,EMIT,load_data,stage,label
from recycling_spectra import weighted_moment
ROOT=Path(__file__).parent
MULT=['E1','M1','E2','M2','E3','M3']
THRESHOLD=1.0219979

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def bricc_query(home,z,kev):
    env={**os.environ,'BrIccHome':str(home)}
    r=subprocess.run([str(home/'briccs'),'-Z',str(z),'-g',f'{kev:.6f}'],env=env,cwd=home,
        capture_output=True,text=True,check=True,timeout=5)
    tree=ET.fromstring(r.stdout)
    pair=tree.find("AllCC[@Shell='IPF']")
    if pair is None:raise ValueError((z,kev,r.stdout))
    vals={e.attrib['Mult']:float(e.text) for e in pair.findall('MultCC') if 'Mult' in e.attrib}
    assert set(MULT)<=set(vals),(z,kev,vals)
    return [vals[m] for m in MULT]

def generate(home,out,step=.025):
    # Finer threshold spacing; grid aligned so halving --step nests the grids.
    grid=np.unique(np.r_[np.linspace(1.1,1.3,int(round(.2/(step/5)))+1),
                         np.linspace(1.3,2.,int(round(.7/(step/2.5)))+1),
                         np.linspace(2.,8.,int(round(6/step))+1)])
    tables={}
    for z in range(18,30):
        tables[str(z)]=[bricc_query(home,z,e*1000) for e in grid]
        print('coefficient table Z',z,'nodes',len(grid),flush=True)
    d={'scope':__doc__,'mult':MULT,'energy_MeV':grid.tolist(),'alpha_pair_over_gamma':tables,
       'source_manifest':json.loads((home/'download-manifest.json').read_text()),
       'executable_version':'BrIccS v2.3d (13-Sep-2022)',
       'pair_dataset':'ScPcc (Schlueter and Soff 1979), 1.1–8 MeV',
       'step_MeV':step,'interpolation':'linear between generated coefficients'}
    out.write_bytes(gzip.compress(json.dumps(d,separators=(',',':'),allow_nan=False).encode(),mtime=0))
    return d

class Coefficients:
    def __init__(self,d):
        self.grid=np.array(d['energy_MeV']);self.data={int(z):np.array(v) for z,v in d['alpha_pair_over_gamma'].items()}
    def values(self,z,e):
        e=np.asarray(e);rows=self.data[int(z)]
        return np.stack([np.interp(e,self.grid,rows[:,i],left=0.,right=0.) for i in range(len(MULT))],axis=-1)

# Every function below integrates the normalized conditional spectrum, then
# multiplies separately interpolated yield and cross section at incident energy.
BASE=['photon_count','photon_energy_MeV','below_pair_count','near_threshold_count','tabulated_count','above_table_count',
      'below_pair_energy_MeV','near_threshold_energy_MeV','tabulated_energy_MeV','above_table_energy_MeV']
KEYS=BASE+[f'{kind}_{m}' for kind in ['pair_ratio','pair_branch','pair_transition_energy_MeV','pair_kinetic_MeV'] for m in MULT]

def functions(e,z,c):
    e=np.asarray(e);alpha=c.values(z,e);branch=alpha/(1+alpha)
    masks=np.stack([e<THRESHOLD,(e>=THRESHOLD)&(e<1.1),(e>=1.1)&(e<=8.),e>8.],axis=-1).astype(float)
    return np.concatenate([np.ones(e.shape+(1,)),e[...,None],masks,masks*e[...,None],alpha,branch,
                           branch*e[...,None],branch*np.maximum(e-THRESHOLD,0)[...,None]],axis=-1)

def fold_spectrum(sub,z,c,order=6):
    assert sub['LAW']==1 and sub['LEP'] in (1,2)
    assert all(sub['INT']==2) and all(sub['yields']['INT']==2)
    gx,gw=leggauss(order);out=[]
    for nd,ep,b in zip(sub['ND'],sub['Ep'],sub['b']):
        nd=int(nd);x=np.asarray(ep)*1e-6;y=np.asarray(b[:,0]);assert np.all(x>=0)
        val=np.sum(functions(x[:nd],z,c)*y[:nd,None],axis=0)
        # Break at every spectrum node, threshold, and coefficient knot. This
        # makes spectrum density and alpha piecewise linear on each interval.
        if len(x)>nd+1:
            xx=x[nd:];yy=y[nd:]*1e6 # ENDF density per eV -> per MeV
            grid=np.unique(np.r_[xx,c.grid[(c.grid>xx[0])&(c.grid<xx[-1])],
                         np.array([THRESHOLD,1.1,8.])[(np.array([THRESHOLD,1.1,8.])>xx[0])&(np.array([THRESHOLD,1.1,8.])<xx[-1])]])
            lo=grid[:-1,None];hi=grid[1:,None];half=(hi-lo)/2;e=(lo+hi)/2+half*gx
            if sub['LEP']==1:
                iy=np.minimum(np.searchsorted(xx,e,side='right')-1,len(xx)-2);density=yy[iy]
            else:density=np.interp(e,xx,yy)
            val+=np.sum(functions(e,z,c)*(density*half*gw)[...,None],axis=(0,1))
        out.append(val)
    out=np.asarray(out)
    return {'ZA':0,'state':int(sub['LIP']),'energy_MeV':(sub['E']*1e-6).tolist(),
       'yield_energy_MeV':(sub['yields']['Eint']*1e-6).tolist(),'yield':sub['yields']['yi'].tolist(),
       **{k:out[:,i].tolist() for i,k in enumerate(KEYS)}}

def ledger(root,coeff,order=6):
    archive=root/'beam-data/tendl2025-proton-recycling.zip';c=Coefficients(coeff)
    data=load_data(archive);base=json.loads((root/'network-ion-budget.json').read_text())
    parser=EndfParserFactory.create(array_type='numpy');stages=[]
    with zipfile.ZipFile(archive) as zf,tempfile.TemporaryDirectory() as tmp:
        for za in HOT:
            flow=base['network']['processor_pn_rates'][label(za)]
            if flow==0:continue
            z,a=za;s=stage(za,data[za],base['T_MeV']);name=f'p-{SYMBOL[z]}{a:03d}.tendl'
            path=Path(tmp)/name;path.write_bytes(zf.read(name));mts=[e['mt'] for e in s['events']]
            sections=parser.parsefile(str(path),include=[(6,mt) for mt in mts+[5]])[6];events=[]
            for ev in s['events']:
                mt=ev['mt'];sec=sections[mt];assert sec['LCT']==3
                values={k:0. for k in KEYS};resZ=ev['residual'][0]
                for sub in sec['subsection']:
                    if sub['ZAP']!=0:continue
                    folded=fold_spectrum(sub,resZ,c,order)
                    for k in KEYS:values[k]+=weighted_moment(data[za]['MF3'][mt],folded,base['T_MeV'],z,a,k)/s['sv']
                events.append({'mt':mt,'residual':list(ev['residual']),'events_per_stage_pn':ev['rate'],
                    'coefficient_Z_proxy':resZ,'per_stage_pn':values})
            total={k:sum(e['per_stage_pn'][k] for e in events) for k in KEYS}
            discard={k:0. for k in KEYS}
            assert sections[5]['LCT']==3
            for sub in sections[5]['subsection']:
                if sub['ZAP']!=0:continue
                folded=fold_spectrum(sub,z+1,c,order)
                for k in KEYS:discard[k]+=weighted_moment(data[za]['MF3'][5],folded,base['T_MeV'],z,a,k)/s['sv']
            stages.append({'ZA':list(za),'pn_per_primary':flow,'per_stage_pn':total,'events':events,
                'discarded_MT5_photons_per_stage_pn':discard,'discarded_MT5_coefficient_Z_proxy':z+1})
            print('folded',za,'photons',total['photon_count'],'E1 pairs',total['pair_branch_E1'],flush=True)
    total={k:sum(s['pn_per_primary']*s['per_stage_pn'][k] for s in stages) for k in KEYS}
    discarded={k:sum(s['pn_per_primary']*s['discarded_MT5_photons_per_stage_pn'][k] for s in stages) for k in KEYS}
    return {'scope':__doc__,'T_MeV':base['T_MeV'],'quadrature_order':order,'threshold_MeV':THRESHOLD,
      'coefficients_scope':'Final-residual charge proxy; pure multipoles assigned uniformly to existing photon spectra. No E0.',
      'photon_associated_ratio_definition':'Integral of alpha_pair/gamma times the evaluated gamma yield.',
      'branch_scenario_definition':'Integral of alpha/(1+alpha) treating the evaluated gamma cascade as a hypothetical gamma-only transition cascade.',
      'no_revised_plant_allocation':True,'per_primary':total,'stages':stages,
      'discarded_MT5_photons_per_primary':discarded,
      'discarded_MT5_note':'Photons from the existing high-energy discard branch; kept separate from the exclusive-network energy identity. Compound-target charge is a proxy.',
      'inputs':{'archive_sha256':sha(archive),'network_ion_budget_sha256':sha(root/'network-ion-budget.json')},
      'assumed_unrecovered_40MV_positron_work_MeV_per_primary':{m:40*total['pair_branch_'+m] for m in MULT}}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=ROOT);ap.add_argument('--bricc-home',type=Path,default=Path('/tmp/protium-sources/bricc'))
    ap.add_argument('--coefficients',type=Path,default=ROOT/'bricc-pair-coefficients.json.gz')
    ap.add_argument('--generate',action='store_true');ap.add_argument('--step',type=float,default=.025);ap.add_argument('--order',type=int,default=6)
    ap.add_argument('--out',type=Path,default=ROOT/'prompt-pair-supply.json');args=ap.parse_args()
    coeff=generate(args.bricc_home,args.coefficients,args.step) if args.generate else json.loads(gzip.decompress(args.coefficients.read_bytes()))
    result=ledger(args.root,coeff,args.order);result['inputs']['coefficient_sha256']=sha(args.coefficients)
    args.out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(result['per_primary'],indent=2))
