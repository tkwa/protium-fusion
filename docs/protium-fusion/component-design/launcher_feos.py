"""Bounded nonideal merger-entropy screen, Frankfurt FEOS16.7 + primary FPEOS.
Run: python launcher_feos.py --cache /writable/path/launcher-feos
Requires NumPy/SciPy,g++; GPLv3 upstream archive retained unchanged beside data.
No hydrodynamic, grain-contact or yield validation. No extrapolation of FPEOS.
"""
from pathlib import Path
import argparse,subprocess,tarfile,ctypes as c,os,json,hashlib,sys,math
import numpy as np
from scipy.optimize import root,brentq
P=Path(__file__).resolve().parent;ROOT=P.parent
ap=argparse.ArgumentParser();ap.add_argument('--cache',type=Path,required=True);a=ap.parse_args();a.cache=a.cache.resolve();a.cache.mkdir(parents=True,exist_ok=True)
arc=P/'launcher-data/FEOS.tar.gz';build=a.cache/'FEOS';lib=build/'Code/libfeos.so'
if not lib.exists():
 with tarfile.open(arc) as tar:
  for member in tar.getmembers():
   dest=(a.cache/member.name).resolve()
   if not dest.is_relative_to(a.cache) or member.issym() or member.islnk():raise ValueError('unsafe archive member')
  tar.extractall(a.cache)
 files=['COMMON-01_UTILITIES.C','COMMON-02_READFILE.C','LIB-01_TF_SERVICES.C','LIB-02_TF_TABLE_1.C','LIB-02_TF_TABLE_2.C','LIB-03_TF_MIXTURE.C','LIB-04_IONMOD.C','LIB-05_EOS_SERVICES.C','LIB-06_MAXWELL.C','LIB-07_C_INTERFACE.C','LIB-08_FORTRAN_INTERFACE.C']
 done=subprocess.run(['g++','-O2','-fPIC','-shared',*files,'-o','libfeos.so'],cwd=build/'Code',capture_output=True,text=True);(a.cache/'build.log').write_text(done.stdout+done.stderr);done.check_returncode()
db=build/'EOS-Data/FEOS_Material-DB.dat';text=db.read_text().split('\n% launcher additions')[0]+'\n% launcher additions\n'
for num,Z,A,H,rho in [(9001,6,12,2,.9442),(9002,6,13,1,1.12915)]:
 text+=f'\nMaterial-{num}:\n[{num}]_SESAME-Number = 1000\n[{num}]_Treference = 0.025852\n[{num}]_Rhoreference = {rho}\n[{num}]_Bulk-Modulus = 2.0e10\n[{num}]_Number-of-Elements = 2\n[{num}]_A[1] = {A}\n[{num}]_Z[1] = {Z}\n[{num}]_X[1] = 1\n[{num}]_A[2] = 1\n[{num}]_Z[2] = 1\n[{num}]_X[2] = {H}\n'
db.write_text(text);os.chdir(build/'EOS-Data')
l=c.CDLL(str(lib));D=c.c_double;I=c.c_int;DP=c.POINTER(D);IP=c.POINTER(I)
l.FEOS_Initialize.argtypes=[I,I];l.FEOS_Init_Mat.argtypes=[I,I,I,I,DP,I,IP,IP,DP];l.FEOS_Get_EOS.argtypes=[I,I,I,D,D,IP,DP,DP,DP,DP,DP,DP];l.FEOS_Initialize(2,0)
for j,n in enumerate([9001,9002],1):
 ne=I();ni=I();tt=D();l.FEOS_Init_Mat(j,n,0,0,None,0,c.byref(ne),c.byref(ni),c.byref(tt))
AR=7.5657e-16

def state(j,rho,T):
 b=I();p=D();e=D();s=D();f=D();q=(D*2)();qt=D()
 l.FEOS_Get_EOS(j,0,0,rho/1000,T/11604.51812,c.byref(b),c.byref(p),c.byref(e),c.byref(s),c.byref(f),q,c.byref(qt))
 return dict(pressure_Pa=p.value*.1+AR*T**4/3,internal_J_kg=e.value*1e-4+AR*T**4/rho,entropy_J_kg_K=s.value*1e-4/11604.51812+4*AR*T**3/(3*rho),ion_charges=list(q),charge_per_formula=qt.value)

def shock(j,rho0,q):
 base=state(j,rho0,300.)
 def fun(x):
  rho,T=np.exp(x);s=state(j,rho,T)
  return [(s['internal_J_kg']-base['internal_J_kg']-q)/q,(.5*(s['pressure_Pa']+base['pressure_Pa'])*(1/rho0-1/rho)-q)/q]
 rr=root(fun,np.log([rho0*4,max(1e4,q/8000)]));rho,T=np.exp(rr.x);err=float(max(abs(np.array(fun(rr.x)))))
 if err>1e-5 or rho<=rho0:raise ValueError('no shockroot')
 return dict(rho_kg_m3=rho,T_K=T,q_J_kg=q,compression=rho/rho0,**state(j,rho,T),root_residual=err)

sys.path.insert(0,str(P));from launcher_screen import Table
checks=[];out={};sources={}
for stage,j,rho0,path,table in [('B2',1,944.2,ROOT/'b2-completion/B2-64-compression-cap.initial.json','C14H28_EOS_09-18-20.txt'),('C13',2,1129.15,ROOT/'inputs/pulse-profiles/scaled-hot-C13-64-R30-initial.json','C18H18_EOS_09-18-20.txt')]:
 d=json.loads(path.read_text());sources[str(path.relative_to(ROOT))]=hashlib.sha256(path.read_bytes()).hexdigest();m=np.array(d['mass_fraction_of_total']);r=np.r_[0,d['outer_radius_m']];rho=1.1309733552923255e13*m/(4*np.pi/3*np.diff(r**3));T=np.array(d['T_K']);v=np.array(d['velocity_m_s']);nm=np.r_[.5*(m[:-1]+m[1:]),.5*m[-1]];K=float(.5*nm@(v*v));S=np.array([state(j,rr,tt)['entropy_J_kg_K'] for rr,tt in zip(rho,T)]);meanS=float(m@S)
 qmatch=math.exp(brentq(lambda x:shock(j,rho0,math.exp(x))['entropy_J_kg_K']-meanS,math.log(1e7),math.log(1e11)))
 nmatch=brentq(lambda n:K*(2/n-4/(3*n*n))-qmatch,100,1e6)
 rows=[];fpt=Table(table,14)
 for N in [100,1000,2000,3000,10000,100000]:
  q=K*(2/N-4/(3*N*N));sh=shock(j,rho0,q);fp=fpt.shock(q);mapped=[]
  for rr in rho:mapped.append(math.exp(brentq(lambda x:state(j,rr,math.exp(x))['entropy_J_kg_K']-sh['entropy_J_kg_K'],math.log(10),math.log(3e9))))
  mapped=np.array(mapped);row=dict(effective_parallel_directions=N,shock=sh,primary_FPEOS_shock=fp,mapped_initial_temperature_K=mapped.tolist(),mapped_T_massmean_K=float(m@mapped),mapped_T_ratio_massmean=float(m@(mapped/T)))
  if fp['success']:row['FEOS_over_FPEOS_T']=sh['T_K']/fp['temperature_K'];row['FEOS_over_FPEOS_compression']=sh['compression']/fp['compression']
  rows.append(row);checks += [sh['root_residual']<1e-5,np.all(mapped>0),np.isfinite(mapped).all()]
 # Fixed angularcount is not a uniformentropy source: retain actual radialvelocity distribution.
 localK=.25*(v*v+np.r_[0.,v[:-1]**2]);active=(np.cumsum(m)-.5*m)>=.001;w=m[active]/sum(m[active]);local_rows=[];required=[]
 for ii in np.where(active)[0]:
  qq=math.exp(brentq(lambda x:shock(j,rho0,math.exp(x))['entropy_J_kg_K']-S[ii],math.log(1e7),math.log(1e11)))
  # Large-N branch of q/K=2/N-4/(3N²); analytic inversion, if a real branch exists.
  ratio=qq/localK[ii];nn=(1+math.sqrt(1-4*ratio/3))/ratio if ratio<=.75 else None
  required.append(nn)
 def quant(values,weight,probs=(.1,.5,.9)):
  order=np.argsort(values);vv=np.array(values)[order];ww=np.array(weight)[order];return np.interp(probs,np.cumsum(ww)-.5*ww,vv).tolist()
 for N in [1000,2000,3000,10000]:
  ratios=[];temps=[]
  for ii in np.where(active)[0]:
   sh=shock(j,rho0,localK[ii]*(2/N-4/(3*N*N)))
   mt=math.exp(brentq(lambda x:state(j,rho[ii],math.exp(x))['entropy_J_kg_K']-sh['entropy_J_kg_K'],math.log(10),math.log(3e9)));temps.append(mt);ratios.append(mt/T[ii])
  local_rows.append(dict(effective_parallel_directions=N,nonseed_mass_fraction=float(sum(m[active])),local_temperature_ratio_massmean=float(w@ratios),local_temperature_ratio_mass_q10_q50_q90=quant(ratios,w),local_temperature_ratio_minmax=[float(min(ratios)),float(max(ratios))],active_cell_indices=np.where(active)[0].tolist(),mapped_T_K=temps))
 local_requirements=dict(scope='Use local prescribedhotstateK asproxy fordeliveryK; actualcoldstream mustalso supplycompressivework. Omitinner0.1%mass seed. EffectiveN canvarywithradialshell/time; prescribeduniformN producesnonuniformentropy.',specific_K_mass_closure_error=float(m@localK/K-1),required_effective_N_for_reference_entropy_mass_q10_q50_q90=quant(required,w),required_effective_N_range=[min(required),max(required)],fixed_N_profiles=local_rows)
 checks.append(abs(m@localK/K-1)<1e-12)
 # Independent differential identities at contact and hot reference, relativeerrors archived.
 identities=[]
 for rr,tt in [(rho0*4,1e5),(rho0*4,1e6),(float(m@rho),float(m@T))]:
  h=1e-4;s0=state(j,rr,tt);sp=state(j,rr,tt*(1+h));sm=state(j,rr,tt*(1-h));rp=state(j,rr*(1+h),tt);rm=state(j,rr*(1-h),tt)
  dedt=(sp['internal_J_kg']-sm['internal_J_kg'])/(2*h*tt);dsdt=(sp['entropy_J_kg_K']-sm['entropy_J_kg_K'])/(2*h*tt);dpdt=(sp['pressure_Pa']-sm['pressure_Pa'])/(2*h*tt);dedr=(rp['internal_J_kg']-rm['internal_J_kg'])/(2*h*rr);rhs=(s0['pressure_Pa']-tt*dpdt)/rr**2
  serr=abs(tt*dsdt-dedt)/max(abs(dedt),1);eerr=abs(dedr-rhs)/max(abs(dedr),abs(rhs),1)
  identities.append(dict(rho_kg_m3=rr,T_K=tt,relative_TdS_dE_error=serr,relative_first_law_density_error=eerr));checks += [serr<.02,eerr<.02]
 out[stage]=dict(prepared_specific_kinetic_J_kg=K,initial_FEOS_entropy_min_mean_max_J_kg_K=[float(min(S)),meanS,float(max(S))],entropy_matching_shock=shock(j,rho0,qmatch),entropy_matching_effective_direction_count=nmatch,shock_mappings=rows,local_speed_entropy=local_requirements,thermodynamic_identity_checks=identities)
result={'scope':'Optimistic scalar-entropy comparison: equal-area parallel-sector merger shock then isentropic compression to prescribed hot densities. It does not deliver actual positions/velocities, certify grain shock homogenization, or predict burn yield.','model':'Frankfurt FEOS16.7 Thomas-Fermi/Cowan/bonding, no Maxwell/softsphere, plus equilibrium radiation. 300K reference and2GPa bulk modulus arechosenpolymerparameters.','limitations':['Thomas-Fermi misses discrete shell structure; compare FPEOS contact states.','FEOS published source emits compilerwarnings including legacy mismatched new/delete; no library qualification claimed.','Ideal reversible compression is a lower entropy bound; later shocks/3D asymmetry can add entropy.','Mapped uniform entropy is only compared at original spatialdensityprofile; totalenergy restoration and new hydro would be subsequent tests.'],'source_urls':['https://doi.org/10.1016/j.cpc.2018.01.008','https://doi.org/10.17632/6vjsv6v48p.1','https://militzer.berkeley.edu/FPEOS/'],'archive_sha256':hashlib.sha256(arc.read_bytes()).hexdigest(),'material_database_sha256':hashlib.sha256(db.read_bytes()).hexdigest(),'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'profile_sha256':sources,'stages':out,'checks':{'count':len(checks),'passed':int(sum(checks)),'all_passed':bool(all(checks))}}
(P/'launcher-feos.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result['checks']));print(json.dumps({s:{'Nentropy':v['entropy_matching_effective_direction_count'],'q':v['entropy_matching_shock']['q_J_kg'],'T':v['entropy_matching_shock']['T_K']} for s,v in out.items()}))
