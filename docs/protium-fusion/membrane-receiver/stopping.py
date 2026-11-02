"""C13 receiver membrane screen. No hydro rerun, empirical stopping + explicit tails.
Run with PYTHONPATH=/tmp/protium-python python stopping.py --root WORKTREE.
Outputs are local to this file. stopping-data.json preserves selected IAEA rows.
"""
from pathlib import Path
import argparse,json,hashlib,math,re
import numpy as np
from scipy.integrate import cumulative_trapezoid
HERE=Path(__file__).resolve().parent
MU=1.66053906660e-27; EV=1.602176634e-19; C=299792458.; NA=6.02214076e23
# Densities and I are screening choices; compounds use elemental electron-weighted log I.
EL={'H':(1,1.008,19.2),'Be':(4,9.0122,63.7),'B':(5,10.81,76.),'C':(6,12.011,78.),'N':(7,14.007,82.),'O':(8,15.999,95.),'Si':(14,28.0855,173.)}
MAT={'silica':{'rho':2200.,'parts':{'Si':1,'O':2},'db':'SiO2','I':139.2},'carbon':{'rho':2000.,'parts':{'C':1},'db':'Carbon','I':78.},'Be':{'rho':1848.,'parts':{'Be':1},'db':'Be','I':63.7},'BN':{'rho':2100.,'parts':{'B':1,'N':1},'db':None,'I':79.5},'SiC':{'rho':3210.,'parts':{'Si':1,'C':1},'db':'SiC','I':117.}}
DATA=json.loads((HERE/'stopping-data.json').read_text())['selected_data']

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def atom(n):
 if n=='p':return 'H',1,1
 if n=='n':return 'n',0,1
 k=re.match(r'([a-z]+)(\d+)',n);s=k[1].capitalize();return s,{'He':2,'C':6,'N':7,'O':8,'F':9,'Ne':10,'Na':11,'Mg':12}[s],int(k[2])
def target(m):
 p=MAT[m]['parts'];ms=sum(v*EL[k][1] for k,v in p.items());zs=sum(v*EL[k][0] for k,v in p.items())
 return ms,zs,sum(p.values())
def bethe(eu,z,m):
 # Output MeV/(mg/cm2). Only used >=0.2 MeV/u, or matched high tails.
 eu=np.asarray(eu);b2=1-(1+eu/931.49410242)**-2
 ze=z*(1-np.exp(-125*np.sqrt(b2)/z**(2/3)))
 ms,zs,_=target(m);I=MAT[m]['I']*1e-6
 return .307075/1000*(zs/ms)*ze**2/b2*(np.log(1.0219979*b2/(1-b2)/I)-b2)
def analytic(eu,z,m):
 # Low-energy sqrt(E) completion from0.2MeV/u is a model, not Bethe there.
 eu=np.asarray(eu);return np.where(eu>=.2,bethe(np.maximum(eu,.2),z,m),float(bethe(.2,z,m))*np.sqrt(eu/.2))
def measured(eu,s,m):
 rr=[r for r in DATA if r['target_name']==MAT[m]['db'] and r['projectile_name']==s]
 if not rr:return None
 # Average duplicate abscissae; log interpolation of actual measured points.
 raw={}
 for r in rr:
  assert r['energy_unit']=='MeV/u';assert r['stopping_unit_converted']=='MeV/(mg/cm2)'
  raw.setdefault(float(r['energy']),[]).append(float(r['stopping_power_converted']))
 x=np.array(sorted(raw));y=np.array([np.mean(raw[q]) for q in x]);e=np.asarray(eu)
 mid=np.exp(np.interp(np.log(e),np.log(x),np.log(y)));z={'C':6,'N':7,'O':8}[s]
 lo=y[0]*np.sqrt(e/x[0]);hi=y[-1]*analytic(e,z,m)/float(analytic(x[-1],z,m))
 return np.where(e<x[0],lo,np.where(e>x[-1],hi,mid))
def electronic(eu,z,m):
 if z in (6,7,8):
  s={6:'C',7:'N',8:'O'}[z]
  if m=='silica' and z==7:
   # N measurements stop at52keV/u. Blend to neighboring measured C/O by0.1MeV/u.
   low=measured(eu,'N',m);high=np.sqrt(measured(eu,'C',m)*measured(eu,'O',m))
   a=np.clip(np.log(np.asarray(eu)/.05218571)/np.log(.1/.05218571),0,1)
   return low**(1-a)*high**a
  got=measured(eu,s,m)
  if got is not None:return got
 return analytic(eu,z,m)
def nuclear(E,z,a,m):
 # ZBL universal screened elastic stopping, mixture additivity; E inMeV perion.
 result=np.zeros_like(np.asarray(E));ms,_,_=target(m)
 for sym,nu in MAT[m]['parts'].items():
  z2,a2,_=EL[sym];zz=z**.23+z2**.23
  eps=32.53*a2*(np.asarray(E)*1000)/(z*z2*(a+a2)*zz)
  sn=np.where(eps<=30,np.log1p(1.1383*eps)/(2*(eps+.01321*eps**.21226+.19593*np.sqrt(eps))),np.log(eps)/(2*eps))
  cross=8.462e-15*z*z2*a/((a+a2)*zz)*sn #eV cm2 per targetatom
  result+=cross*(NA*nu/ms)*1e-9 #MeVcm2/mg
 return result

def source(root,stage='C13'):
    path=root/'docs/protium-fusion/membrane-receiver/pulse-loading.json';shared=json.loads(path.read_text());d=shared['stages'][stage];p=d['spectral_packets']
    aa=np.asarray(shared['A']);zz=np.asarray(shared['Z']);area=d['area_m2'];vf=np.asarray(p['velocity_m_s']);wf=np.asarray(p['mass_fluence_kg_m2'])*area
    ys=np.asarray(p['Y_at_wall']);ns=wf[:,None]/MU*ys;ee=np.asarray(p['kinetic_MeV_per_u'])[:,None]*aa
    cells=d['cells'];bounds=np.asarray(cells['velocity_bounds_m_s']);ci=np.asarray(p['cell_index'],int)
    # Eight-point velocity quadrature is documented in the shared transport packet.
    g,gw=np.polynomial.legendre.leggauss(8);dvw=np.diff(bounds)[ci]*np.tile(gw/2,len(bounds)-1)
    km=np.sum(ns*ee)*EV*1e6;mass=np.sum(wf)
    compact={'stage':stage,'species':shared['species'],'radius_m':shared['radius_m'],'mass_kg':mass,'final_kinetic_J':float(km),'final_internal_J':d['endpoint']['final_internal_J'] if 'final_internal_J' in d['endpoint'] else 0.,'speed_quadrature_weights_m_s':dvw.tolist(),'checks':{'baryon_mass_relative_error':float(np.sum(ns*aa)*MU/mass-1),'smooth_kinetic_J':float(km)},'hashes':{'pulse-loading.json':sha(path)},'source_scope':'Finite-flight beta-decayed 8-point-per-cell Gauss packets; smooth kinetic energy preserved, no artificial normalization to nodal kinetic energy.'}
    assert abs(compact['checks']['baryon_mass_relative_error'])<1e-10
    return compact,ns,ee,vf,wf,aa,zz

def table(m,z,a,factor=1.):
 E=np.geomspace(1e-10,100.,16000);se=electronic(E/a,z,m)*factor;sn=nuclear(E,z,a,m);st=se+sn
 R=cumulative_trapezoid(1/st,E,initial=0) #mg/cm2
 D=cumulative_trapezoid(sn/st,E,initial=0) #MeV nucleartransfer
 return E,R,D,se,sn

def run_case(m,thick_nm,ns,ee,aa,zz,area,factor=1.,brake_keV_u=None,velocities=None,dv_weights=None,radius=1.1e9):
 sigma=MAT[m]['rho']*thick_nm*1e-9;arealmg=sigma*100 #1kg/m2=100mg/cm2
 totalN=np.sum(ns)/area;inputE=0.;loss=0.;nl=0.;stoppedM=0.;stoppedN=0.;sp={};peaklet=0;loss_samples=np.zeros(len(ns))
 for j,z in enumerate(zz):
  if z==0 or np.sum(ns[:,j])<np.sum(ns)*1e-12:continue
  a=int(aa[j]);E0=ee[:,j] if brake_keV_u is None else np.minimum(ee[:,j],brake_keV_u*a/1000)
  Eg,R,D,se,sn=table(m,int(z),a,factor);r0=np.interp(E0,Eg,R);rem=np.maximum(r0-arealmg,0);E1=np.interp(rem,R,Eg);E1=np.where(rem==0,0,E1);dE=E0-E1;dn=np.interp(E0,Eg,D)-np.interp(E1,Eg,D);n=ns[:,j]/area
  ei=np.sum(n*E0)*EV*1e6;de=np.sum(n*dE)*EV*1e6;dl=np.sum(n*dn)*EV*1e6;sm=np.sum(n*(rem==0))*a*MU
  loss_samples+=n*dE*EV*1e6
  inputE+=ei;loss+=de;nl+=dl;stoppedM+=sm;stoppedN+=np.sum(n*(rem==0));peaklet=max(peaklet,float(np.max(se))*MAT[m]['rho']/1e4) #MeV/(mgcm-2)*rho[gcm-3]*1e-1 keV/nm
  if np.sum(n)>totalN*1e-5:sp[str(j)]={'incident_number_m2':float(sum(n)),'input_J_m2':float(ei),'loss_J_m2':float(de),'stopped_mass_kg_m2':float(sm),'energy_loss_fraction':float(de/ei),'stopped_number_fraction':float(np.sum(n*(rem==0))/sum(n))}
 ms,_,nt=target(m);atoms=sigma/(ms*MU)*nt
 out={'material':m,'thickness_nm':thick_nm,'areal_mass_kg_m2':sigma,'electronic_multiplier':factor,'residual_cap_keV_u':brake_keV_u,'incident_J_m2':inputE,'projectile_energy_loss_J_m2':loss,'energy_loss_fraction':loss/inputE,'energy_transmitted_fraction':1-loss/inputE,'elastic_transfer_J_m2':nl,'stopped_mass_kg_m2':stoppedM,'stopped_mass_fraction':stoppedM/(np.sum(ns*aa)*MU/area),'stopped_number_fraction':stoppedN/totalN,'target_atoms_m2':atoms,'incident_ions_per_target_atom':totalN/atoms,'NRT_proxy_dpa_Ed40_eV':.8*nl/EV/(2*40*atoms),'uniform_sputter_yield_for_1pct_atom_loss':.01*atoms/totalN,'maximum_electronic_stopping_keV_nm':peaklet,'species':sp}
 if velocities is not None:
  flux=loss_samples/dv_weights*velocities**2/radius;i=int(np.argmax(flux));out.update({'peak_projectile_loss_flux_W_m2':float(flux[i]),'peak_projectile_loss_time_s':float(radius/velocities[i])})
 return out

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ar=ap.parse_args();d,ns,ee,v,w,a,z=source(ar.root)
 area=4*np.pi*d['radius_m']**2;kw={'velocities':v,'dv_weights':np.asarray(d['speed_quadrature_weights_m_s']),'radius':d['radius_m']};mf=d['mass_kg']/area;ef=d['final_kinetic_J']/area;N=ns.sum()/area
 rows=[]
 for m in MAT:
  for t in [10.,30.,100.,300.,1000.,3000.,10000.]:
   rows.append(run_case(m,t,ns,ee,a,z,area,**kw))
 sens=[run_case('silica',t,ns,ee,a,z,area,f,**kw) for t in [10.,100.,1000.,10000.] for f in [.5,2.]]
 brakes=[run_case(m,t,ns,ee,a,z,area,1.,e) for m in ['silica','carbon'] for t in [10.,100.,1000.] for e in [.01,.1,1.,10.]]
 q=[]
 for wt,label in [(w,'mass'),(np.sum(ns*ee,axis=1),'energy')]:
  ix=np.argsort(v);cs=np.cumsum(wt[ix]);q.append({'weight':label,'velocity_quantiles_1_50_99_m_s':np.interp(np.array([.01,.5,.99])*cs[-1],cs,v[ix]).tolist(),'energy_per_u_quantiles_MeV':(.5*MU*np.interp(np.array([.01,.5,.99])*cs[-1],cs,v[ix])**2/EV/1e6).tolist()})
 coverage=[]
 for m in MAT:
  for k in ['C','N','O']:
   rr=[r for r in DATA if r['projectile_name']==k and r['target_name']==MAT[m]['db']]
   coverage.append({'material':m,'projectile':k,'measured_min_MeV_u':min((float(r['energy']) for r in rr),default=None),'measured_max_MeV_u':max((float(r['energy']) for r in rr),default=None),'datasets':sorted({r['ref_id'] for r in rr})})
 out={'scope':'Normal-incidence CSDA projectile energy removal by intact membranes; no electron/recoil escape or membrane survival claim. Finite-flight beta-decayed wall isotope fractions from shared pulse-loading packet. Species-dependent measured electronic stopping with labeled interpolation/extrapolation; ZBL elastic transfer.','pulse':{'mass_fluence_kg_m2':mf,'kinetic_fluence_J_m2':ef,'raw_positive_charge_if_fully_stripped_C_m2':float(np.sum(ns*z)/area*EV),'cold_comoving_electron_kinetic_J_m2':float(np.sum(ns*z*ee/a)/area*EV*1e6*(9.1093837139e-31/MU)),'ion_fluence_m2':N,'heavy_ion_fluence_m2':float(ns[:,z>=6].sum()/area),'ions_per_nm2':N/1e18,'quantiles':q,'composition_mass':{n:float(ns[:,i].sum()*a[i]*MU/d['mass_kg']) for i,n in enumerate(d['species']) if ns[:,i].sum()*a[i]*MU/d['mass_kg']>1e-7}},'materials':MAT,'measurement_coverage':coverage,'direct_rows':rows,'silica_electronic_factor_sensitivity':sens,'upstream_braking_residual_cap_rows':brakes,'checks':d['checks'],'input_sha256':{**d['hashes'],'stopping-data.json':sha(HERE/'stopping-data.json')},'code_sha256':sha(__file__)}
 allrows=[]
 for stage in ['A','B1','B2']:
  sd,sn,se,sv,sw,sa,sz=source(ar.root,stage);skw={'velocities':sv,'dv_weights':np.asarray(sd['speed_quadrature_weights_m_s']),'radius':sd['radius_m']}
  for m in MAT:
   for t in [10.,100.,1000.,10000.]:
    row=run_case(m,t,sn,se,sa,sz,area,**skw);row['stage']=stage;allrows.append(row)
 out['other_stages_direct_rows']=allrows
 (HERE/'stopping-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
 print(json.dumps(out['pulse'],indent=2));print('silica',[(r['thickness_nm'],r['energy_loss_fraction'],r['stopped_mass_fraction'],r['NRT_proxy_dpa_Ed40_eV']) for r in rows if r['material']=='silica'])
if __name__=='__main__':main()
