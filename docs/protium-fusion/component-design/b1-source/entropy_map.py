"""B1-specific nonideal contact entropy map; no 3D assembly claim.
FEOS16.7 mixture, chosen condensed hydrazine reference, equilibrium radiation.
The actual B1 mass/profile/duration are used, rather than the A/C13 density.
"""
from pathlib import Path
import argparse,ctypes as c,os,json,math,hashlib,shutil
import numpy as np
from scipy.optimize import root,brentq
P=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--cache',type=Path,required=True);ap.add_argument('--library-cache',type=Path,required=True);ap.add_argument('--output',type=Path,default=P/'entropy-map.json');a=ap.parse_args();a.cache=a.cache.resolve();a.cache.mkdir(parents=True,exist_ok=True);a.output=a.output.resolve();lib=a.library_cache.resolve()/'FEOS/Code/libfeos.so'
shutil.copyfile(a.library_cache/'FEOS/EOS-Data/FEOS_TF-Table_1197.dat',a.cache/'FEOS_TF-Table_1197.dat')
recipes=[dict(name='H2-rho1062p5',h=2,rho=1062.5,T=250.,folder='grid',tag='2'),dict(name='H2-rho850',h=2,rho=850.,T=250.,folder='grid',tag='2'),dict(name='H2-rho1200',h=2,rho=1200.,T=250.,folder='grid',tag='2'),dict(name='H1p5-rho1000',h=1.5,rho=1000.,T=250.,folder='refine',tag='1p5')]
db=''
for j,z in enumerate(recipes,1):
 z['entity']=j;n=9200+j
 db+=f'\nMaterial-{n}:\n[{n}]_SESAME-Number = 1000\n[{n}]_Treference = {z["T"]/11604.51812}\n[{n}]_Rhoreference = {z["rho"]/1000}\n[{n}]_Bulk-Modulus = 2.0e10\n[{n}]_Number-of-Elements = 2\n[{n}]_A[1] = 15\n[{n}]_Z[1] = 7\n[{n}]_X[1] = 1\n[{n}]_A[2] = 1\n[{n}]_Z[2] = 1\n[{n}]_X[2] = {z["h"]}\n'
(a.cache/'FEOS_Material-DB.dat').write_text(db);os.chdir(a.cache)
l=c.CDLL(str(lib));D=c.c_double;I=c.c_int;DP=c.POINTER(D);IP=c.POINTER(I)
l.FEOS_Initialize.argtypes=[I,I];l.FEOS_Init_Mat.argtypes=[I,I,I,I,DP,I,IP,IP,DP];l.FEOS_Get_EOS.argtypes=[I,I,I,D,D,IP,DP,DP,DP,DP,DP,DP];l.FEOS_Initialize(len(recipes),0)
for j in range(1,len(recipes)+1):
 ne=I();ni=I();tt=D();l.FEOS_Init_Mat(j,9200+j,0,0,None,0,c.byref(ne),c.byref(ni),c.byref(tt))
AR=7.5657e-16

def st(j,rho,T):
 b=I();p=D();e=D();s=D();f=D();q=(D*2)();qt=D();l.FEOS_Get_EOS(j,0,0,rho/1000,T/11604.51812,c.byref(b),c.byref(p),c.byref(e),c.byref(s),c.byref(f),q,c.byref(qt))
 return dict(p=p.value*.1+AR*T**4/3,e=e.value*1e-4+AR*T**4/rho,s=s.value*1e-4/11604.51812+4*AR*T**3/(3*rho),q=list(q))
def shock(z,q):
 j=z['entity'];rho0=z['rho'];s0=st(j,rho0,z['T'])
 def f(x):
  rho,T=np.exp(x);s=st(j,rho,T);return [(s['e']-s0['e']-q)/q,(.5*(s['p']+s0['p'])*(1/rho0-1/rho)-q)/q]
 rr=root(f,np.log([rho0*4,max(1e3,q/1e4)]));rho,T=np.exp(rr.x);err=float(max(abs(np.array(f(rr.x)))))
 if err>1e-5 or rho<=rho0:raise ValueError(f'shockroot failed {z["name"]} q={q:g} err={err:g}')
 return dict(rho_kg_m3=float(rho),T_K=float(T),**st(j,rho,T),q_J_kg=float(q),error=err)
def quant(x,w):
 x=np.asarray(x);ii=np.argsort(x);return np.interp([.1,.5,.9],np.cumsum(w[ii])-.5*w[ii],x[ii]).tolist()
def required_q(z,entropy):
 return math.exp(brentq(lambda x:shock(z,math.exp(x))['s']-entropy,math.log(1e7),math.log(1e12)))
def effective_n(q,k):
 ratio=q/k;return (1+math.sqrt(1-4*ratio/3))/ratio if 0<ratio<=.75 else None
checks=[];results={}
for z in recipes:
 j=z['entity'];stem=Path('/home/tkwa/src/protium-fusion/output/mixture-runs')/z['folder']/f'B1-H{z["tag"]}-fixed_energy-z64-k0p001-t1';src=stem.with_suffix('.initial.json');inp=stem.with_suffix('.inputs.json');meta=stem.with_suffix('.json');p=json.loads(src.read_text());inputs=json.loads(inp.read_text());base=json.loads(meta.read_text());m=np.array(p['mass_fraction_of_total']);rr=np.r_[0,p['outer_radius_m']];rho=base['mass_kg']*m/(4*np.pi/3*np.diff(rr**3));T=np.array(p['T_K']);v=np.array(p['velocity_m_s']);K=.25*(v*v+np.r_[0.,v[:-1]**2]);S=np.array([st(j,r,t)['s'] for r,t in zip(rho,T)]);active=(np.cumsum(m)-.5*m)>=.001;w=m[active]/sum(m[active]);ii=np.where(active)[0]
 rows=[]
 for idx in ii:
  rq={str(f):required_q(z,st(j,rho[idx],T[idx]*f)['s']) for f in [.8,.9,1.,1.1,1.2]}
  rows.append(dict(index=int(idx),specific_hot_radial_K_J_kg=float(K[idx]),hot_density_kg_m3=float(rho[idx]),reference_T_K=float(T[idx]),reference_entropy_J_kg_K=float(S[idx]),required_shock_q_J_kg_by_target_Tfactor=rq,reference_entropy_effective_N=effective_n(rq['1.0'],K[idx]),reference_pair_half_angle_degrees=math.degrees(math.asin(math.sqrt(min(rq['1.0']/K[idx],1.))))))
  checks.append(all(x>0 for x in rq.values()))
 mappings=[]
 for N in [50,100,300,1000,2000]:
  cells=[]
  for idx in ii:
   sh=shock(z,K[idx]*(2/N-4/(3*N*N)));temp=math.exp(brentq(lambda x:st(j,rho[idx],math.exp(x))['s']-sh['s'],math.log(1),math.log(3e9)));cells.append(dict(index=int(idx),mapped_T_K=temp,ratio=temp/T[idx],shock=sh));checks.append(sh['error']<1e-5)
  ratios=[c['ratio'] for c in cells];mappings.append(dict(effective_N=N,temperature_ratio_massmean=float(w@ratios),temperature_ratio_mass_q10_q50_q90=quant(ratios,w),cells=cells))
 Nreq=[c['reference_entropy_effective_N'] for c in rows];angles=[c['reference_pair_half_angle_degrees'] for c in rows];qref=np.array([c['required_shock_q_J_kg_by_target_Tfactor']['1.0'] for c in rows]);results[z['name']]=dict(recipe=z,stage='B1',H_per_catalyst=z['h'],mass_kg=base['mass_kg'],reference_profile=str(src),profile_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),reference_inputs=str(inp),inputs_sha256=hashlib.sha256(inp.read_bytes()).hexdigest(),reference_metadata=str(meta),metadata_sha256=hashlib.sha256(meta.read_bytes()).hexdigest(),initial_Y=inputs['initial_Y'],physical_duration_s=inputs['duration_s'],active_mass_fraction=float(sum(m[active])),reference_entropy_massmean_J_kg_K=float(m@S),reference_q_active_massmean_J_kg=float(w@qref),reference_q_over_active_hot_K=float((w@qref)/(w@K[active])),reference_N_mass_q10_q50_q90=quant(Nreq,w) if all(n is not None for n in Nreq) else None,reference_N_missing_count=sum(n is None for n in Nreq),reference_pair_angle_mass_q10_q50_q90=quant(angles,w),required_cells=rows,fixed_N_mappings=mappings)
 print(json.dumps({k:results[z['name']][k] for k in ['H_per_catalyst','reference_q_active_massmean_J_kg','reference_q_over_active_hot_K','reference_N_mass_q10_q50_q90','reference_pair_angle_mass_q10_q50_q90']}),flush=True)
 print(json.dumps({str(c['effective_N']):c['temperature_ratio_mass_q10_q50_q90'] for c in mappings}),flush=True)
out=dict(scope='Nonideal entropy and collision-energy screen at actual B1 density. Prescribed hot-profile specific K is a proxy for collision energy. Actual cold geometry, radial work and grain homogenization remain unproved; no new external heater. Pair angle assumes equal opposing transverse velocities with q=K sin²(alpha).',cold_reference_scope='H2 rho1062.5 from natural hydrazine1000 times34/32 isotope-volume assumption; T250K and2GPa bulk modulus chosen. Density sensitivity850–1200; H1.5 is a homogenized N/H mixture reference, no qualified phase recipe. FEOS Thomas-Fermi/Cowan lacks molecular phase/fragmentation kinetics.',code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),FEOS_library_sha256=hashlib.sha256(lib.read_bytes()).hexdigest(),material_database_sha256=hashlib.sha256(db.encode()).hexdigest(),cases=results,checks=dict(count=len(checks),passed=int(sum(checks)),all_passed=bool(all(checks))))
a.output.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
