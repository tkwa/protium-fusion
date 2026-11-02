"""Selected H-rich recipe FEOS merger entropy; direct methane FPEOS check.
Uses the library built by launcher_feos.py; no frozen hydro/result modified.
"""
from pathlib import Path
import argparse,ctypes as c,os,json,math,hashlib,shutil
import numpy as np
from scipy.optimize import root,brentq
from launcher_screen import Table
P=Path(__file__).resolve().parent;ROOT=P.parent
ap=argparse.ArgumentParser();ap.add_argument('--cache',type=Path,required=True);ap.add_argument('--library-cache',type=Path,required=True);a=ap.parse_args();a.cache=a.cache.resolve();a.cache.mkdir(parents=True,exist_ok=True);lib=a.library_cache.resolve()/'FEOS/Code/libfeos.so'
shutil.copyfile(a.library_cache/'FEOS/EOS-Data/FEOS_TF-Table_1197.dat',a.cache/'FEOS_TF-Table_1197.dat')
recipe={'C13H4':{'stage':'C13','h':4,'A':13,'Z':6,'rho0':439*17/16,'T0':102.4,'file':Path('/home/tkwa/src/protium-fusion/output/mixture-runs/grid/C13-H4-fixed_energy-z64-k0p001-t1')},'AH3':{'stage':'A','h':3,'A':14,'Z':7,'rho0':850.,'T0':150.,'file':Path('/home/tkwa/src/protium-fusion/output/mixture-runs/refine/A-H3-fixed_energy-z64-k0p001-t1')}}
db='';
for j,z in enumerate(recipe.values(),1):
 n=9100+j;z['entity']=j
 db+=f'\nMaterial-{n}:\n[{n}]_SESAME-Number = 1000\n[{n}]_Treference = {z["T0"]/11604.51812}\n[{n}]_Rhoreference = {z["rho0"]/1000}\n[{n}]_Bulk-Modulus = 2.0e10\n[{n}]_Number-of-Elements = 2\n[{n}]_A[1] = {z["A"]}\n[{n}]_Z[1] = {z["Z"]}\n[{n}]_X[1] = 1\n[{n}]_A[2] = 1\n[{n}]_Z[2] = 1\n[{n}]_X[2] = {z["h"]}\n'
(a.cache/'FEOS_Material-DB.dat').write_text(db);os.chdir(a.cache)
l=c.CDLL(str(lib));D=c.c_double;I=c.c_int;DP=c.POINTER(D);IP=c.POINTER(I)
l.FEOS_Initialize.argtypes=[I,I];l.FEOS_Init_Mat.argtypes=[I,I,I,I,DP,I,IP,IP,DP];l.FEOS_Get_EOS.argtypes=[I,I,I,D,D,IP,DP,DP,DP,DP,DP,DP];l.FEOS_Initialize(2,0)
for j in [1,2]:
 ne=I();ni=I();tt=D();l.FEOS_Init_Mat(j,9100+j,0,0,None,0,c.byref(ne),c.byref(ni),c.byref(tt))
AR=7.5657e-16

def st(j,rho,T):
 b=I();p=D();e=D();s=D();f=D();q=(D*2)();qt=D();l.FEOS_Get_EOS(j,0,0,rho/1000,T/11604.51812,c.byref(b),c.byref(p),c.byref(e),c.byref(s),c.byref(f),q,c.byref(qt))
 return dict(p=p.value*.1+AR*T**4/3,e=e.value*1e-4+AR*T**4/rho,s=s.value*1e-4/11604.51812+4*AR*T**3/(3*rho),q=list(q))
def shock(z,q):
 j=z['entity'];rho0=z['rho0'];s0=st(j,rho0,z['T0'])
 def f(x):
  rho,T=np.exp(x);s=st(j,rho,T);return [(s['e']-s0['e']-q)/q,(.5*(s['p']+s0['p'])*(1/rho0-1/rho)-q)/q]
 rr=root(f,np.log([rho0*4,max(1e4,q/1e4)]));rho,T=np.exp(rr.x);err=float(max(abs(np.array(f(rr.x)))))
 if err>1e-5:raise ValueError('shockroot failed')
 return dict(rho_kg_m3=rho,T_K=T,**st(j,rho,T),q_J_kg=q,error=err)
def quant(x,w):
 ii=np.argsort(x);return np.interp([.1,.5,.9],np.cumsum(w[ii])-.5*w[ii],np.array(x)[ii]).tolist()
checks=[];results={}
for name,z in recipe.items():
 j=z['entity'];src=z['file'].with_suffix('.initial.json');inp=z['file'].with_suffix('.inputs.json');p=json.loads(src.read_text());inputs=json.loads(inp.read_text());m=np.array(p['mass_fraction_of_total']);rr=np.r_[0,p['outer_radius_m']];rho=(4*np.pi/3*1e8*30**3)*m/(4*np.pi/3*np.diff(rr**3));T=np.array(p['T_K']);v=np.array(p['velocity_m_s']);K=.25*(v*v+np.r_[0.,v[:-1]**2]);S=np.array([st(j,r,t)['s'] for r,t in zip(rho,T)]);active=(np.cumsum(m)-.5*m)>=.001;w=m[active]/sum(m[active]);N=2000;rows=[]
 for ii in np.where(active)[0]:
  sh=shock(z,K[ii]*(2/N-4/(3*N*N)));temp=math.exp(brentq(lambda x:st(j,rho[ii],math.exp(x))['s']-sh['s'],math.log(5),math.log(3e9)));rows.append(dict(index=int(ii),mapped_T_K=temp,old_T_K=T[ii],ratio=temp/T[ii],shock=sh));checks.append(sh['error']<1e-5)
 meanq=float(m@K)*(2/N-4/(3*N*N));meansh=shock(z,meanq);meanS=float(m@S);qm=math.exp(brentq(lambda x:shock(z,math.exp(x))['s']-meanS,math.log(1e7),math.log(1e11)));km=float(m@K);nm=brentq(lambda n:km*(2/n-4/(3*n*n))-qm,100,1e6)
 comparison=None
 if name=='C13H4':
  tab=Table('C10H40_EOS_09-18-20.txt',17);native_rho=tab.m/tab.V0;tab.V0=tab.m/z['rho0'];fp=tab.shock(meanq)
  comparison=dict(scope='DirectCH4firstprinciplesP/E;C13inertialmass17u/formula. InitialV resettochosen102Kliquiddensity,referencecoldE unchanged;omits smallcoldcohesive/pressureenergy difference. Originaltableinitialdensityretainedbelow.',original_FPEOS_reference_rho_kg_m3=native_rho,cold_rho_kg_m3=z['rho0'],FPEOS=fp,FEOS_over_FPEOS_T=meansh['T_K']/fp['temperature_K'] if fp['success'] else None,FEOS_over_FPEOS_compression=(meansh['rho_kg_m3']/z['rho0'])/fp['compression'] if fp['success'] else None)
 results[name]=dict(stage=z['stage'],H_per_catalyst=z['h'],effective_N=N,cold_reference={'rho_kg_m3':z['rho0'],'T_K':z['T0'],'bulk_modulus_Pa':2e9,'pressure_Pa':0.,'density_scope':'C13methane naturalCH4density439×17/16 isotopevolumeassumption' if name=='C13H4' else 'ChosenNH3solidmatrix850kg/m3;notmeasuredforactualcapsule'},initial_Y=inputs['initial_Y'],physical_duration_s=inputs['duration_s'],profile_source=str(src),profile_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),inputs_source=str(inp),inputs_sha256=hashlib.sha256(inp.read_bytes()).hexdigest(),initial_entropy_massmean_J_kg_K=meanS,mean_shock=meansh,mean_entropy_matching_N=nm,active_mass_fraction=float(sum(m[active])),temperature_ratio_massmean=float(w@np.array([r['ratio'] for r in rows])),temperature_ratio_mass_q10_q50_q90=quant([r['ratio'] for r in rows],w),cells=rows,primary_comparison=comparison)
 checks += [len(rows)>0,np.all(np.isfinite([r['mapped_T_K'] for r in rows]))]
out={'scope':'Selectedrecipe local FEOS entropy mapping; C13H4directFPEOScontactcomparison;AH3hasnoNH3firstprinciplesmixturecheckyet. PrescribedhotstateKproxy,inner0.1%seedunchanged.','code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'material_database_sha256':hashlib.sha256(db.encode()).hexdigest(),'FEOS_library_sha256':hashlib.sha256(lib.read_bytes()).hexdigest(),'cases':results,'checks':{'count':len(checks),'passed':int(sum(checks)),'all_passed':bool(all(checks))}}
(P/'launcher-mixture-entropy.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out['checks']));print(json.dumps({k:{'Nmatch':v['mean_entropy_matching_N'],'Tratio':v['temperature_ratio_mass_q10_q50_q90'],'primary':v['primary_comparison']} for k,v in results.items()},indent=2))
