"""A/C13 preparatory nuclear and radiation screen on exact hot source fields.
Ballistic backward geometry plus FEOS constant-entropy temperatures estimates
transport scales; it omits pressure work and is not a delivered trajectory.
"""
from pathlib import Path
import json,sys,time
import numpy as np
from scipy.integrate import simpson
P=Path(__file__).resolve().parent;ROOT=P.parent.parent
sys.path[:0]=[str(ROOT),str(ROOT/'b2-completion'),str(ROOT/'mixture-design')]
from initial_screen import load,sha,quant
from feos_contact import FEOS
from backward_isentrope import Adiabats
import radial_completion as radial
import pulse_model as pm
import updated_network
feos=FEOS();rate=updated_network.install(radial,quantile='median');out={};begin=time.monotonic()
for stage in ['A','C13']:
 profile,inp,base,r,v,m,T0,rho0,ye,yi=load(stage);n=len(m);w=m/sum(m);ad=Adiabats(stage,rho0,T0,m,feos,301);t=np.linspace(0,.0005,501);radii=r[None,:]-t[:,None]*v[None,:];rho=m/(4*np.pi/3*np.diff(radii**3,axis=1));T=np.array([ad.temperature(row) for row in rho]);parent=np.array(inp['source_refinement']['parent_cell_for_each_new_cell']);src=json.loads(Path(inp['source_refinement']['source_initial']).read_text());sw=np.array(src['mass_fraction_of_total']);active=np.isin(parent,np.where(np.cumsum(sw)-.5*sw>=.001)[0])
 cv=np.array([(feos.state(stage,rho0[j],T0[j]*1.0001)['u_J_kg']-feos.state(stage,rho0[j],T0[j]*.9999)['u_J_kg'])/(T0[j]*.0002) for j in range(n)]);rad=[];er=7.5657e-16*T**4;rc=.5*(radii[:,1:]+radii[:,:-1]);rf=.5*(rho[:,1:]+rho[:,:-1]);ef=.5*(er[:,1:]+er[:,:-1]);grad=np.diff(er,axis=1)/np.diff(rc,axis=1)
 for kappa in [.001,.01,.03]:
  rr=np.abs(grad)/(kappa*rf*np.maximum(ef,1e-300));lim=(2+rr)/(6+3*rr+rr**2);flux=-299792458*lim*grad/(kappa*rf);fouter=299792458*er[:,-1]/(4+3*kappa*rho[:,-1]*(radii[:,-1]-rc[:,-1]));lum=np.c_[np.zeros(len(t)),4*np.pi*radii[:,1:-1]**2*flux,4*np.pi*radii[:,-1]**2*fouter];q=np.diff(lum,axis=1)/m;signed=simpson(-q/T,x=t,axis=0)/cv;absolute=simpson(np.abs(q)/T,x=t,axis=0)/cv;loss=simpson(lum[:,-1],x=t);ftrans=simpson(np.sum(np.abs(lum[:,1:-1]),axis=1),x=t)
  rad.append({'opacity_m2_kg':kappa,'outer_escape_J':float(loss),'outer_escape_fraction_prepared':float(loss/base['initial_energy_J']),'sum_absolute_internal_face_transfer_fraction_prepared':float(ftrans/base['initial_energy_J']),'max_abs_equivalent_hot_dT_over_T_all':float(max(absolute)),'max_abs_equivalent_hot_dT_over_T_nonseed':float(max(absolute[active])),'nonseed_mass_quantiles_abs_equivalent_hot_dT_over_T':quant(absolute[active],w[active]),'nonseed_mass_quantiles_signed_equivalent_hot_dT_over_T':quant(signed[active],w[active]),'total_mass_above_1pct_abs_equivalent_hot_dT':float(sum(w[absolute>.01])),'total_mass_above_10pct_abs_equivalent_hot_dT':float(sum(w[absolute>.1])),'cells':[{'index':j,'mass_fraction':float(w[j]),'seed':bool(not active[j]),'signed_equivalent_hot_dT_over_T':float(signed[j]),'abs_equivalent_hot_dT_over_T':float(absolute[j])} for j in range(n)]})
 y=np.array([inp['initial_Y'].get(n,0.) for n in pm.SPECIES]);cat='n14' if stage=='A' else 'c13';ci=pm.SPECIES.index(cat);tempmax=float(max(T0));rhomax=float(max(rho0));dy,qq,qnu,qc=radial.network(np.array([tempmax]),np.array([rhomax/1000]),y[None,:]);duration=.001;frac=abs(dy[0,ci])*duration/y[ci];qe=(qq[0]-qnu[0])*1e-4*duration*sum(m)
 result={'scope':__doc__,'stage':stage,'source128_initial_sha256':sha(Path('/home/tkwa/src/protium-fusion/output/source128')/f'{stage}-H{3 if stage=="A" else 4}-source128-carbon-median.initial.json'),'profile_hydrogen_per_catalyst':3 if stage=='A' else 4,'duration_s':float(t[-1]),'time_samples':len(t),'adiabat_nodes':301,'pressure_work_omitted':True,'density_minmax_kg_m3':[float(rho.min()),float(rho.max())],'temperature_minmax_K':[float(T.min()),float(T.max())],'source_seed_mass_fraction':float(sum(w[~active])),'radiation':rad,'nuclear_hot_state_envelope':{'T_K':tempmax,'rho_kg_m3':rhomax,'held_time_s':duration,'catalyst_changed_fraction':float(frac),'nuclear_deposition_fraction_prepared':float(qe/base['initial_energy_J']),'scope':'All material held for1ms at the largest actual hot-profile temperature and density, using fixed initial composition and median-carbon rates. A deliberately hot/dense initial-rate envelope conditional on preparatory states remaining below those maxima; not a full network evolution or an absolute bound on unmodeled seed formation.'},'adiabat_consistency_max_relative_energy_difference':max(g['work_vs_FEOS_delta_u_relative'] for g in ad.groups)}
 raw=Path('/home/tkwa/src/protium-fusion/output/assembly-audit')/f'{stage}-ballistic-isentropic-screen.npz';np.savez_compressed(raw,time_s=t,radius_m=radii,density_kg_m3=rho,temperature_K=T,mass_kg=m);result['trace_path']=str(raw);result['trace_sha256']=sha(raw);out[stage]=result
 print(stage,{k:v for k,v in result.items() if k not in ['radiation']},flush=True);print([{k:v for k,v in d.items() if k!='cells'} for d in rad],flush=True)
record={'scope':__doc__,'code_sha256':sha(__file__),'rate_update':rate,'FEOS_library_sha256':sha(feos.library),'cases':out,'wall_seconds':time.monotonic()-begin};(P/'other-stages-screen.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
