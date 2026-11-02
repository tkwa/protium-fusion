"""Conditional radial graphite recapture honeycomb: ray transport and interfaces.
Sticking is irreversible incorporation over the required operating interval;
transient adsorption followed by desorption is not counted as capture.
"""
from pathlib import Path
import json,hashlib,math
import numpy as np
P=Path(__file__).resolve().parent;seed=104721;root=Path('/home/tkwa/.codex/worktrees/protium-fusion/website');base=root/'docs/protium-fusion/component-design/catcher-followup';rho=2000.;pitch=.01;H=.1
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()

def cosdir(rng,n,axis,sgn):
 u=rng.random(n);ph=rng.uniform(0,2*np.pi,n);a=np.sqrt(u);b=np.sqrt(1-u);d=np.empty((n,3));d[:,axis]=a*sgn;j=(axis+1)%3;k=(axis+2)%3;d[:,j]=b*np.cos(ph);d[:,k]=b*np.sin(ph);return d

def mc(N,stick,diffuse,seed,source="floor"):
 rng=np.random.default_rng(seed);pos=np.column_stack([rng.random(N),rng.random(N),np.full(N,0. if source=='floor' else 10.)]);v=cosdir(rng,N,2,1. if source=='floor' else -1.);idx=np.arange(N);escaped=0;captured=0;direct=0;iters=0;collision_sum=0;exit_coll=[];angular_sum=0.;angular2_sum=0.;escaped_cos=[]
 while len(idx):
  # distance to each enclosing plane; topz=10 is open, floorz=0 closed.
  bound=np.where(v>0,np.array([1.,1.,10.]),np.zeros(3));dist=np.divide(bound-pos,v,out=np.full_like(v,np.inf),where=np.abs(v)>1e-15);dist=np.where(dist>1e-12,dist,np.inf);ax=np.argmin(dist,axis=1);steps=dist[np.arange(len(pos)),ax];pos+=v*steps[:,None];out=(ax==2)&(v[:,2]>0)
  mu=v[out,2];sn=np.sqrt(np.maximum(1-mu*mu,1e-30));aw=np.arcsin(mu)/(2*sn);angular_sum+=float(sum(aw));angular2_sum+=float(sum(aw*aw));escaped_cos.extend(mu.tolist())
  escaped+=int(sum(out));direct+=int(sum(out)) if iters==0 else 0;exit_coll.extend([iters]*int(sum(out)))
  ids=np.flatnonzero(~out);idx=idx[ids];pos=pos[ids];v=v[ids];ax=ax[ids]
  if not len(idx):break
  collision_sum+=len(idx);absorbed=rng.random(len(idx))<stick;captured+=int(sum(absorbed));ids=np.flatnonzero(~absorbed);idx=idx[ids];pos=pos[ids];v=v[ids];ax=ax[ids]
  if not len(idx):break
  rnd=rng.random(len(idx))<diffuse
  for axis in range(3):
   aa=np.flatnonzero((ax==axis)&rnd)
   if len(aa):
    # inward signs opposite incidentnormalcomponent (floor always +z).
    signs=-np.sign(v[aa,axis]);v[aa]=cosdir(rng,len(aa),axis,signs)
  aa=np.flatnonzero(~rnd);v[aa,ax[aa]]*=-1
  pos+=1e-11*v # numerical push inward, negligible versus unitcell.
  iters+=1
  if iters>10000:raise RuntimeError('uncaptured residual')
 assert escaped+captured==N
 p=escaped/N;z=1.959963984540054;den=1+z*z/N;center=(p+z*z/(2*N))/den;half=z*math.sqrt(p*(1-p)/N+z*z/(4*N*N))/den
 return {'source':source,'mean_firstflight_radial_column_weight_per_emitted_atom':angular_sum/N,'escaped_mean_radial_column_weight_per_particle':angular_sum/max(escaped,1),'escaped_angle_deg_quantiles_10_50_90':np.degrees(np.arccos(np.quantile(escaped_cos,[.9,.5,.1]))).tolist() if escaped_cos else [],'rays':N,'irreversible_sticking_probability':stick,'diffuse_fraction_among_unstuck':diffuse,'escape_probability':p,'escape_binomial95percent_interval':[center-half,center+half],'direct_escape_probability_MC':direct/N,'captured_probability':captured/N,'mean_wall_encounters':collision_sum/N,'maximum_encounters':iters,'scope':'Ideal squarecell, uniformLambertianfloor emission; diffuse cosine or specular reflection. Wall/floor have same irreversible sticking. Isothermal speed is irrelevant to this geometric eventualescape probability.'}
# Independently integrated exact direct square-to-square view factor.
x,w=np.polynomial.legendre.leggauss(80);x=(x+1)/2;w=w/2;xx,yy=np.meshgrid(x,x,indexing='ij');ww=np.outer(w,w);view=float(4/math.pi*np.sum(ww*(1-xx)*(1-yy)*100/(100+xx*xx+yy*yy)**2))
# Direct paths have a computable carbon-column lower bound, even perfectcapture.
angles_direct=[]
xxq,wq=np.polynomial.legendre.leggauss(160);ph=(xxq+1)*np.pi/8;wp=wq*np.pi/8
for aspect in [10.,30.,100.]:
 integ=0.;fv=0.
 for az,waz in zip(ph,wp):
  thmax=np.arctan(1/(aspect*np.cos(az)));th=(xxq+1)*thmax/2;wt=wq*thmax/2;accept=(1-aspect*np.tan(th)*np.cos(az))*(1-aspect*np.tan(th)*np.sin(az))
  integ+=8/(2*np.pi)*waz*np.sum(wt*np.cos(th)*(np.pi/2-th)*accept)
  fv+=8/(2*np.pi)*waz*np.sum(wt*2*np.cos(th)*np.sin(th)*accept)
 angles_direct.append({'aspect_H_over_p':aspect,'direct_escape_probability':float(fv),'direct_radial_column_weight_per_floor_emission':float(integ),'mean_direct_column_weight_per_escaped_atom':float(integ/fv),'scope':'Independentangularquadrature ofLambertianfloor andfiniteopeningacceptance. Normalizedtoonebarecosine-emittedsphericalflight. It is anunavoidablefirstflight contributionfor thisstraightgeometry; returningflights add.'})
rows=[]
for s in [.1,.25,.5,1.]:
 for a in ([0.,.5,1.] if s<1 else [1.]):rows.append(mc(1_000_000,s,a,seed+int(s*1000)+int(a*10)))
incoming=[]
for st in [.1,.25,.5,1.]:
 for a in ([0.,.5,1.] if st<1 else[1.]):incoming.append(mc(500000,st,a,seed+20000+int(st*1000)+int(a*10),'opening'))
geom=[];angles=[];escaped=[]
for t in [5e-6,1e-5,1e-4]:
 f=1-(pitch/(pitch+t))**2;mass=rho*H*f
 geom.append({'wall_thickness_m':t,'clear_cell_width_m':pitch,'depth_m':H,'projected_tip_fraction':f,'added_graphite_mass_kg_m2':mass,'sidewall_area_per_receiver_area':4*pitch*H/(pitch+t)**2,'axial_mass_column_kg_m2_on_tips':rho*H,'linear_tip_to_base_resistance_m2_K_W_per_tip_area_at_k50':H/50,'total_tip_conductance_per_receiver_area_W_m2_K_at_k50':f*50/H})
 for th in [0.,1e-5,1e-4,1e-3]:
  for phi in [0.,math.pi/4]:
   dx=H*math.tan(th)*math.cos(phi);dy=H*math.tan(th)*math.sin(phi);intercept=1-max(0,pitch-abs(dx))*max(0,pitch-abs(dy))/(pitch+t)**2
   angles.append({'wall_thickness_m':t,'beam_to_cell_axis_angle_rad':th,'azimuth_rad':phi,'intercepted_primary_mass_fraction':intercept,'additional_interception_above_tips':intercept-f,'unrecovered_H_He_fraction_upper_if_all_intercepted_reemits':intercept,'scope':'Uniform parallel radial beam and straight aligned walls. Reemissionfraction andwallgasrecovery multiply this upper routingfraction; lateralions returningtobacking mayreduce it.'})
 for r in rows:
  for ratio in [1.,3.,10.]:
   esc=(1-f)*r['escape_probability']+f*ratio
   escaped.append({'wall_thickness_m':t,'sticking_probability':r['irreversible_sticking_probability'],'diffuse_fraction':r['diffuse_fraction_among_unstuck'],'tip_yield_relative_to_floor':ratio,'whole_face_escaping_yield_relative_to_bare_floor':esc,'suppression_factor':1/esc,'scope':'Outward Lambertian tip emission escapes directly. Equal incident fluence per projected floor/tip area; tilt/otherirradiatedsidewalls areadditional.'})
checks=[]
for r in rows:
 assert abs(r['captured_probability']+r['escape_probability']-1)<1e-15;checks.append('probability_closure')
 assert r['escape_probability']>=view-5*math.sqrt(view*(1-view)/r['rays']);checks.append('direct_lower_bound')
r=rows[-1];assert abs(r['escape_probability']-view)<5*math.sqrt(view*(1-view)/r['rays']);checks.append('direct_view_independent_integral')
assert all(x['added_graphite_mass_kg_m2']<5 for x in geom);checks.append('mass')
out={'scope':__doc__,'code_sha256':sha(__file__),'seed':seed,'direct_opening_view_factor_exact_quadrature':view,'far_opening_approximation':pitch*pitch/(math.pi*H*H),'direct_flight_radial_column_bounds':angles_direct,'returning_cosine_carbon':incoming,'ray_transport':rows,'geometry':geom,'tilt_and_gas_interception':angles,'whole_face_yield_sensitivity':escaped,'checks_passed':len(checks),'checks':checks,'limits':['Sticking is a conditional permanentincorporation probability, not a measured irradiatedgraphitevalue.','Finite carbonresidence,chemistry,clusters,sidewallerosion,dust andtemperaturedependentsticking areunmodeled.','Primarytip/tiltedwallirradiation anditsH/Hegasrouting remain explicit separate sources.','The hoodcollimatesescapingcarbon. Escapeprobabilityalonecannotmultiplythebarecosinecavitycolumn; exactfirstflightangularweightsandreturncapturemustbeaccounted.','The geometry cannot be placed over a molecularpump mouth without recomputingconductance.','Statistical95percentintervals quantifyMonteCarlosamplingonly; they are not subjectiveengineeringcredibleintervals.']}
(P/'results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'direct':view,'rays':rows,'geometry':geom,'checks':len(checks)},indent=2))
