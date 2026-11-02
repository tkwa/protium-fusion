"""Positive-prestress three-family support at reduced receiver masses.
Static construction plus global affine energy/control screen; no local network proof.
The comparison masses keep archived current128 rate and99% energy capture solely
for reference. No capture qualification or selected thinner plant is inferred.
"""
from pathlib import Path
import json,math,sys,hashlib
import numpy as np
P=Path(__file__).parent;BASE=P.parent;sys.path.insert(0,str(BASE/'mixture-design'))
import plant_compare as pc
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
pkt=BASE/'mixture-design/current-rate128-sigma2000.json';red=BASE/'mixture-design/current-rate128-reduced-inputs.json'
j=json.loads(pkt.read_text())['top_detailed'][0];archive=json.loads(red.read_text())['cases'];selected=[next(x for x in archive.values() if x['name']==name) for name in j['plant']['selected_cases'].values()]
R=j['plant']['radius_m'];area=4*math.pi*R*R;G=6.67430e-11;n2=G*1.98847e30/(5*149597870700.)**3;n=math.sqrt(n2);omega=1e4/R;E2D=1e8;nu=.2
# Equal-area sphere quadrature for nonlinear prestress area.
z,w=np.polynomial.legendre.leggauss(40);ph=2*math.pi*np.arange(80)/80
normals=np.array([[math.sqrt(1-zz*zz)*math.cos(p),math.sqrt(1-zz*zz)*math.sin(p),zz]for zz in z for p in ph]);weights=np.repeat(w/2,80)/80
qz,qw=np.polynomial.legendre.leggauss(96);qz=(qz+1)/2;qw=qw/2;theta=2*math.pi*np.arange(1024)/1024

def energy2(E,a,sigma,T):
 ap=2*T/(sigma*R*R);b=np.array([(a+ap+4*n2)/2,(a+ap-2*n2)/2,(a+ap-4*n2)/2]);t=np.array([3*n2,0.,-n2]);tr=np.trace(E);norm=np.sum(E*E);trE2=np.trace(E@E)
 qg=-1.5*a*(trE2/5-2*norm/15+tr*tr/5);qk=0.
 for i in range(3):
  u,v=[x for x in range(3)if x!=i];Q=-.5*(np.sum(E[:,u]**2)+np.sum(E[:,v]**2))+9/8*(E[u,u]**2+E[v,v]**2)+7/4*E[u,u]*E[v,v]+(E[u,v]+E[v,u])**2/8;qk+=b[i]*Q
 qt=-.5*np.sum(t[:,None]*E*E)
 qa=(2*tr*tr+2*norm-3*trE2)/15
 qp=3*T/(sigma*R*R)*qa
 sym=(E+E.T)/2;ts=np.trace(sym);sq=np.trace(sym@sym);lam=E2D*nu/(1-nu*nu);mu=E2D/(2*(1+nu))
 qe=3/(sigma*R*R)*((lam/5+mu/15)*ts*ts+(lam/15+7*mu/15)*sq)
 return qg+qk+qt+qp+qe

def nonlinear(F,a,sigma,T):
 ap=2*T/(sigma*R*R);b=np.array([(a+ap+4*n2)/2,(a+ap-2*n2)/2,(a+ap-4*n2)/2]);t=np.array([3*n2,0.,-n2]);C=F@F.T
 grav=-3*a*sum(ww/math.sqrt(np.linalg.det(np.eye(3)+zz*zz*(C-np.eye(3))))for zz,ww in zip(qz,qw));kin=0.
 for i in range(3):
  u,v=[x for x in range(3)if x!=i];tt=np.zeros((3,len(theta)));tt[u]=np.cos(theta);tt[v]=np.sin(theta);kin+=b[i]/np.mean(np.linalg.norm(F@tt,axis=0))**2
 tide=-.5*np.sum(t[:,None]*F*F)
 ar=np.linalg.det(F)*np.linalg.norm(np.linalg.solve(F.T,normals.T),axis=0)@weights
 return grav+kin+tide+3*T/(sigma*R*R)*(ar-1)

def modes(a,sigma,T,days):
 basis=np.eye(9).reshape(9,3,3);H=np.zeros((9,9))
 for i in range(9):
  H[i,i]=2*energy2(basis[i],a,sigma,T)
  for k in range(i):H[i,k]=H[k,i]=energy2(basis[i]+basis[k],a,sigma,T)-.5*(H[i,i]+H[k,k])
 C=np.zeros((9,9));Pa=np.zeros((9,9))
 for c in range(3):C[c,3+c]=-2*n;C[3+c,c]=2*n
 for k,E in enumerate(basis):Pa[:,k]=((E-E.T)/2).ravel()
 tau=days*86400;H+=Pa/tau**2;C+=1.4*Pa/tau
 ev=np.linalg.eigvals(np.block([[np.zeros((9,9)),np.eye(9)],[-H/a,-C/math.sqrt(a)]]))*math.sqrt(a)
 return {'attitude_days':days,'maximum_growth_per_s':float(max(ev.real)),'growing':bool(max(ev.real)>1e-7*math.sqrt(a)),'H_eigenvalues_over_a':np.linalg.eigvalsh(H/a).tolist()}

rows=[];checks=[];rng=np.random.default_rng(85);old=json.loads((P/'support-topology.json').read_text());oldbase=json.loads((P/'sector-account.json').read_text());oldeps=max(old['geometry']['axis_surface_stream_kg_m2'].values())/oldbase['support']['total_areal_mass_kg_m2']
for sigma_receiver in [500.,1000.,1500.,2000.]:
 plant=pc.station(selected,j['cycle'],j['plant']['helium_per_s'],absorber_sigma=sigma_receiver)
 M=plant['mass_kg'];sigma=M/area;a=G*M/(2*R**3);T=max(1e4,(4.4*n2-a)*sigma*R*R/2);ap=2*T/(sigma*R*R);b=np.array([(a+ap+4*n2)/2,(a+ap-2*n2)/2,(a+ap-4*n2)/2]);sigstream=sigma*b/omega**2
 assert min(sigstream)>0
 # Independent nonlinear energy checks; elastic contribution is subtracted
 # from polynomial because nonlinear() includes only selfgravity/tides/streams/prestress.
 saved=E2D;E2D=0.;h=2e-4;f0=nonlinear(np.eye(3),a,sigma,T);errs=[];first=[]
 for k in range(8):
  E=rng.normal(size=(3,3));E/=np.linalg.norm(E);fp=nonlinear(np.eye(3)+h*E,a,sigma,T);fm=nonlinear(np.eye(3)-h*E,a,sigma,T)
  errs.append(abs((fp+fm-2*f0)/(2*h*h)-energy2(E,a,sigma,T))/a);first.append(abs((fp-fm)/(2*h))/a)
 E2D=saved;assert max(errs)<2e-6 and max(first)<2e-6
 mode=[modes(a,sigma,T,d)for d in [30.,100.,300.]]
 force=sum(b)-b;target=np.array([a-3*n2,a,a+n2])+ap;assert max(abs(force-target))<1e-28
 streammass=float(sigstream.sum()*area);energy=streammass*1e8/3;loss=.001*energy/(2*math.pi/omega)
 row={'receiver_kg_m2_comparison':sigma_receiver,'mass_kg_same_rate99pct_capture':M,'body_kg_m2':sigma,'a_over_n2':a/n2,'effective_a_over_n2':(a+ap)/n2,'pretension_N_m':T,'added_inward_pressure_Pa':2*T/R,'frame_requirement':{'existing_package_kg_m2':20.,'extensional_modulus_N_m':E2D,'Poisson_ratio':nu,'small_biaxial_prestrain':(1-nu)*T/E2D,'ideal_100GPa_2000density_sheet_kg_m2':2.,'stress_in_ideal_1mm_load_layer_Pa':T/.001,'stress_if_8kg_m2_is_loadbearing_at2000density_Pa':T/.004,'scope':'Cold load-path and joint requirements. The20kg/m² frame is retained inside services. An8kg/m² equivalent load layer leaves12kg/m² for joints/attachments, but geometry and material allowable must actually achieve E2D and tension together.'},'stream_axis_kg_m2':sigstream.tolist(),'stream_mass_kg':streammass,'stream_mass_relative_to_existing_scalar_budget':streammass/plant['mass_breakdown_kg']['support_streams'],'stream_loss_W_m2_at_0p001_per_revolution':loss/area,'heaviest_local_epsilon_ratio_to_old_tested':max(sigstream)/sigma/oldeps,'affine_cases':mode,'nonlinear_check_max_second_order_error_over_a':max(errs),'nonlinear_check_max_first_order_residual_over_a':max(first)}
 rows.append(row)
out={'scope':__doc__,'code_sha256':sha(__file__),'input_sha256':{str(p):sha(p)for p in [pkt,red,Path(pc.__file__),P/'support-topology.json',P/'sector-account.json']},'radius_m':R,'n2_s2':n2,'target_minimum_effective_a_over_n2':4.4,'positive_family_limit':4.,'prestress_area_quadratic':'(2 tr(E)^2 + 2 ||E||² - 3 tr(E²))/15','prestress_normalized_energy':'3 T/(sigma R²) times surface-area expansion','rows':rows,'qualification':'Static topology and9DOF affine screen only. Local flexible-guide instability/control remains separate; smaller positive epsilon is compared with its tested range, not proof of a full coupled network. No new plant performance/capture credit. Source, guide hardware, local frame and processing allocations do not shrink with stream mass.'}
(P/'receiver-prestress.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps(rows,indent=2))
