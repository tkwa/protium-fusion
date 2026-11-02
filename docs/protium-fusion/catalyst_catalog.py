"""Screen proton/neutron/weak-decay catalysts using PADF-2 and NUBASE2020.

This is a ranked nuclear-data diagnostic, not a reactor calculation. It does
not credit independently evaluated production channels as exclusive branches.
Nuclide ground states with a single reported beta-plus/EC return are retained;
isomers and mixed/particle-emission branches are flagged separately.
"""
from pathlib import Path
import re,json,zipfile,hashlib,argparse,gzip,tempfile
from math import pi,sqrt
import numpy as np
from numpy.polynomial.legendre import leggauss
import pynucastro as pyn
from read_padf import read
from nonneutral_screen import MEV,EPS0,MU0,E_CHARGE,C
ROOT=Path(__file__).parent
NUBASE=Path(pyn.__file__).parent/'nucdata/AtomicMassEvaluation/nubase_4.mas20.txt'
ME=.51099895
GX,GW=leggauss(16)


def nuclear(Z,A):return pyn.Nucleus.from_Z_A(Z,A)
def atomicmass(Z,A):return nuclear(Z,A).mass
def nucmass(Z,A):return atomicmass(Z,A)-Z*ME

def nubase_records():
 out={}
 for line in NUBASE.read_text().splitlines():
  if line.startswith('#') or not line.strip():continue
  A=int(line[:3]);Z=int(line[4:7]);state=int(line[7])
  out[(Z,A,state)]={'A':A,'Z':Z,'state':state,'name':line[11:17].strip(),
   'half_life_raw':line[69:80].strip(),'decay_modes':line[119:209].strip(),
   'record':line}
 return out


def thermal(table,T,Z,A,power=0):
 alpha=nucmass(Z,A)/(nucmass(Z,A)+nucmass(1,1));mu=nucmass(Z,A)*nucmass(1,1)/(nucmass(Z,A)+nucmass(1,1))*MEV/C**2
 en=alpha*np.array(table['energy_MeV']);xs=np.array(table['sigma_barn'])
 lo=en[:-1,None];hi=en[1:,None];half=(hi-lo)/2
 ee=(lo+hi)/2+half*GX;ss=(xs[:-1,None]+xs[1:,None])/2+(xs[1:,None]-xs[:-1,None])/2*GX
 integ=float(np.sum(half*GW*ss*ee**(1+power)*np.exp(-ee/T)))
 return sqrt(8*MEV/(pi*mu))*T**-1.5*integ*1e-28


def field_upper(sv,Z,A,heat,B=100,fill=.2,strength=1e6,mult=5):
 m=nucmass(Z,A)*MEV/C**2;S=Z*EPS0*B*B/(2*m)
 rate=S*S/(4*Z)*sv;mf=mult*B*B/(2*MU0)/strength
 return fill*rate*heat*MEV/mf


def run(archive):
 nb=nubase_records();meta=[];allcurves={};mp=atomicmass(1,1);mn=atomicmass(0,1)
 qH=(4*mp-atomicmass(2,4))/2
 with zipfile.ZipFile(archive) as z, tempfile.TemporaryDirectory(prefix='protium-padf-') as temp:
  tmp=Path(temp)
  for name in z.namelist():
   mt=re.fullmatch(r'p(\d{2})(\d{3})\.ndf',name)
   if not mt:continue
   Z,A=map(int,mt.groups());target=nb.get((Z,A,0));daughter=nb.get((Z+1,A,0))
   if not target or not daughter:continue
   raw=z.read(name);path=tmp/name;path.write_bytes(raw);d=read(path)
   header=next(s[:66].strip() for s in raw.decode('latin1').splitlines() if 'EVAL-' in s)
   d['evaluation']=header;d['source_file']=name;d['target_Z']=Z;d['target_A']=A
   chans=d['channels'];key=f'{(Z+1)*1000+A}_0';pn=chans.get(key)
   qpn=atomicmass(Z,A)+mp-atomicmass(Z+1,A)-mn
   modes=daughter['decay_modes'];unambiguous=modes in ['B+=100','EC=100']
   stable=target['half_life_raw']=='stbl'
   row={'target':target['name'],'daughter':daughter['name'],'Z':Z,'A':A,
        'target_stable':stable,'daughter_half_life_raw':daughter['half_life_raw'],
        'daughter_decay_modes':modes,'single_ground_state_return':unambiguous,
        'Q_pn_MeV':qpn,'lab_pn_threshold_MeV':max(0,-qpn*(1+nucmass(1,1)/nucmass(Z,A))),
        'decay_Q_atomic_MeV':atomicmass(Z+1,A)-atomicmass(Z,A),
        'source_file':name,'source_sha256':d['source_sha256'],'evaluation':header,
        'pn_ground_table_present':pn is not None,
        'pn_additional_isomers':[k for k in chans if k.startswith(f'{(Z+1)*1000+A}_') and k!=key]}
   qpa=atomicmass(Z,A)+mp-atomicmass(Z-1,A-3)-atomicmass(2,4)
   qpd=atomicmass(Z,A)+mp-atomicmass(Z,A-1)-atomicmass(1,2)
   row['Q_palpha_MeV']=qpa;row['Q_pd_MeV']=qpd
   if pn and unambiguous and stable:
    allcurves[name]=d
    selected={'pn':pn,'p_gamma':chans.get(f'{(Z+1)*1000+A+1}_0'),
              'p_alpha_residual':chans.get(f'{(Z-1)*1000+A-3}_0'),
              'p_d_residual':chans.get(f'{Z*1000+A-1}_0'),
              'inelastic_residual':chans.get(f'{Z*1000+A}_0'),
              'inclusive_neutron':chans.get('1_0'),'inclusive_D':chans.get('1002_0')}
    # Native linear data intervals, zero where no evaluation is present.
    assert all(all(law==2 for _,law in t['interpolation']) for t in selected.values() if t)
    rows=[]
    for T in [.2,.3,.4,.5,.6,.8,1.,1.2,1.5,2.,2.5,3.]:
     sv=thermal(pn,T,Z,A);er=thermal(pn,T,Z,A,1)/sv if sv else None
     ratios={k:thermal(t,T,Z,A)/sv if t is not None and sv else None for k,t in selected.items() if k!='pn'}
     ratios={k:v for k,v in ratios.items()}
     # Zero neutrino energy provides an intentionally loose upper bound;
     # no practical energy conversion or plant loads are counted here.
     upper=field_upper(sv,Z,A,qH)
     rows.append({'ion_T_MeV':T,'pn_ground_sigma_v_m3_s':sv,'reaction_Ecm_MeV':er,
                  'production_ratios_to_pn_ground':ratios,
                  'ideal_nuclear_W_kg_field_mass_only_zero_neutrinos':upper})
    row['thermal_cases']=rows
    row['sigma_pn_barn_at_3_5_7_10_MeV']=np.interp([3,5,7,10],pn['energy_MeV'],pn['sigma_barn']).tolist()
   meta.append(row)
 out={'scope':'Ground-state return catalog and loose field-mass nuclear ceiling; no net energy, no recycling, no full equilibrium or complete plant.',
      'archive':Path(archive).name,'archive_sha256':hashlib.sha256(Path(archive).read_bytes()).hexdigest(),
      'nubase_sha256':hashlib.sha256(NUBASE.read_bytes()).hexdigest(),
      'nubase_citation':'Kondev et al., Chinese Physics C45,030001(2021), NUBASE2020 bundled with pynucastro3.1.0.',
      'field_assumptions':{'B_T':100,'plasma_fill_fraction':.2,'specific_strength_J_kg':1e6,'mass_factor':5,'thermal_efficiency':1,'neutrino_energy':0},
      'candidates':meta}
 (ROOT/'catalyst-catalog.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
 (ROOT/'catalyst-pn-curves.json.gz').write_bytes(gzip.compress((json.dumps({'source':out['archive_sha256'],'evaluations':allcurves},separators=(',',':'),allow_nan=False)+'\n').encode(),mtime=0))
 for row in meta:
  if 'thermal_cases' in row:
   x=next(v for v in row['thermal_cases'] if v['ion_T_MeV']==1.)
   rat=x['production_ratios_to_pn_ground']
   print(row['target'],'->',row['daughter'],'half',row['daughter_half_life_raw'],'Q',round(row['Q_pn_MeV'],3),
         'sv@1',f"{x['pn_ground_sigma_v_m3_s']:.3e}",'alpha',None if rat['p_alpha_residual'] is None else round(rat['p_alpha_residual'],4),
         'pg',None if rat['p_gamma'] is None else round(rat['p_gamma'],4),'inel',None if rat['inelastic_residual'] is None else round(rat['inelastic_residual'],4),flush=True)
 print('targets',len(meta),'single-return candidates',len(allcurves))

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--archive',default=str(ROOT/'beam-data'/'padf2-z6-z41-2026-04-23.zip'));args=ap.parse_args();run(args.archive)
