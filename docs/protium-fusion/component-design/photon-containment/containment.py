"""Carbon photon-containment bounds plus a clearly limited Compton-only transport sensitivity."""
from pathlib import Path
import json, hashlib, math
import numpy as np
# NIST X-Ray Mass Attenuation Coefficients, carbon; E MeV, mu and mu_en cm2/g.
NIST=np.array([
[.001,2211,2209],[.0015,700.2,699],[.002,302.6,301.6],[.003,90.33,89.63],[.004,37.78,37.23],[.005,19.12,18.66],[.006,10.95,10.54],[.008,4.576,4.242],[.010,2.373,2.078],[.015,.8071,.5627],[.020,.442,.2238],[.030,.2562,.06614],[.040,.2076,.03343],[.050,.1871,.02397],[.060,.1753,.02098],[.080,.1610,.02037],[.100,.1514,.02147],[.150,.1347,.02449],[.200,.1229,.02655],[.300,.1066,.02870],[.400,.09546,.02950],[.500,.08715,.02969],[.600,.08058,.02956],[.800,.07076,.02885],[1,.06361,.02792],[1.25,.05690,.02669],[1.5,.05179,.02551],[2,.04442,.02345],[3,.03562,.02048],[4,.03047,.01849],[5,.02708,.01710],[6,.02469,.01607],[8,.02154,.01468],[10,.01959,.01380],[15,.01698,.01258],[20,.01575,.01198]])
SIGMAS=[1000.,2000.,5000.]
ME=.51099895; RE=2.8179403262e-15; NE=6*6.02214076e23/.012011

def coeff(E,col):return np.exp(np.interp(np.log(E),np.log(NIST[:,0]),np.log(.1*NIST[:,col])))
def floor(E):
 # Rounded down from the minimum of the explicit piecewise log-linear NIST curve.
 v=min(float(coeff(E,2)),float(np.min(.1*NIST[NIST[:,0]<=E,2])))
 return math.floor(v*100000)/100000

def kn(E):
 a=E/ME
 return 2*np.pi*RE**2*NE*((1+a)/a**2*(2*(1+a)/(1+2*a)-np.log1p(2*a)/a)+np.log1p(2*a)/(2*a)-(1+3*a)/(1+2*a)**2)

def sample_compton(E,rng):
 n=len(E);eps=np.empty(n);ids=np.arange(n)
 while len(ids):
  e0=1/(1+2*E[ids]/ME);a1=np.log(1/e0);a2=(1-e0**2)/2;r=rng.random((3,len(ids)))
  er=np.where(r[0]<a1/(a1+a2),np.exp(-a1*r[1]),np.sqrt(e0**2+(1-e0**2)*r[1]))
  ct=1-ME/E[ids]*(1/er-1);acc=r[2] <= 1-er/(1+er**2)*(1-ct**2)
  eps[ids[acc]]=er[acc];ids=ids[~acc]
 return eps

def transport(E0,S,N,seed):
 """Normal incidence on thin spherical closed shell; inner exits return elsewhere.
 Free cold electrons only. No photoelectric, Rayleigh, pair or electron transport.
 All residual photons below1keV counted as leakage, not as deposited heat.
 """
 rng=np.random.default_rng(seed);E=np.full(N,E0);x=np.zeros(N);u=np.ones(N);alive=np.ones(N,bool);leak=np.zeros(N);tail=np.zeros(N);deposited=np.zeros(N);direct=np.zeros(N,bool);returns=0;step=0
 while np.any(alive):
  ids=np.flatnonzero(alive);step+=1
  if step>10000:raise RuntimeError('unexpected long Compton tail')
  z=x[ids]-u[ids]*np.log(rng.random(len(ids)))/kn(E[ids]);back=z<0;returns+=int(back.sum());z[back]*=-1;u[ids[back]]*=-1
  out=z>=S;escaped=ids[out];leak[escaped]=E[escaped];alive[escaped]=False
  if step==1:direct[escaped]=True
  ids=ids[~out];x[ids]=z[~out]
  if not len(ids):continue
  eps=sample_compton(E[ids],rng);ct=np.clip(1-ME/E[ids]*(1/eps-1),-1,1);phi=rng.uniform(0,2*np.pi,len(ids));u[ids]=np.clip(u[ids]*ct+np.sqrt(np.maximum(0,1-u[ids]**2))*np.sqrt(np.maximum(0,1-ct**2))*np.cos(phi),-1,1)
  deposited[ids]+=E[ids]*(1-eps);E[ids]*=eps
  done=ids[E[ids]<.001];tail[done]=E[done];alive[done]=False
 assert np.max(abs(leak+tail+deposited-E0))<1e-12*E0
 w=(leak+tail)/E0;se=float(np.std(w,ddof=1)/np.sqrt(N));mean=float(np.mean(w))
 return {'incident_energy_MeV':E0,'sigma_kg_m2':S,'histories':N,'seed':seed,'photon_energy_outer_leakage_fraction':float(np.mean(leak)/E0),'untransported_subkeV_energy_fraction_counted_as_leakage':float(np.mean(tail)/E0),'combined_leakage_fraction':mean,'deposited_fraction':float(np.mean(deposited)/E0),'sampling_standard_error':se,'sampling_80pct_interval_normal_approx':[max(0,mean-1.281551566*se),min(1,mean+1.281551566*se)],'uncollided_fraction_sample':float(np.mean(direct)),'uncollided_fraction_exact':float(np.exp(-kn(E0)*S)),'max_steps':step,'inner_cavity_return_events':returns,'max_energy_balance_error_MeV':float(np.max(abs(leak+tail+deposited-E0)))}

def main():
 sp=json.loads(Path(__file__).with_name('spectrum.json').read_text());hot=json.loads(Path(__file__).with_name('hot-bypass.json').read_text());hotreserve=max(x['all_positive_as_photons_D1p10_lab_bypass_MeV_per_He'] for x in hot['aggregate_cycle']);thermal=sp['cycle']['H_nonnu_MeV_per_He']+sp['cycle']['U_prepared_MeV_per_He'];lines=[]
 for line in sp['primary_photon_lines']:
  E=line['energy_MeV'];mu=float(coeff(E,1));muen=float(coeff(E,2));k=floor(E)
  lines.append({**line,'mu_m2_kg':mu,'mu_en_m2_kg':muen,'minimum_degraded_energy_mu_en_m2_kg':k,'free_electron_KN_m2_kg':float(kn(E)),'thickness_rows':[{'sigma_kg_m2':S,'uncollided_escape_fraction':math.exp(-mu*S),'first_collision_local_deposited_fraction_lower_reference':muen/mu*(-math.expm1(-mu*S)),'localized_cascade_energy_leakage_upper':math.exp(-k*S)} for S in SIGMAS]})
 aggregate=[]
 maxD=max(v['max_final_bulk_Doppler_factor'] for v in sp['stages'].values())
 for S in SIGMAS:
  g=sum(l['MeV_per_He']*math.exp(-l['minimum_degraded_energy_mu_en_m2_kg']*S) for l in lines)
  # Reserve all primary positron kinetic energy as <=5.14304MeV radiation. This is deliberately above physical beta endpoint.
  b=sp['main_beta_kinetic_MeV_per_He']*math.exp(-floor(5.14304)*S)
  nuclear_gamma=sum(l['MeV_per_He']*math.exp(-l['minimum_degraded_energy_mu_en_m2_kg']*S) for l in lines if l['energy_MeV']!=.511)
  anni=sum(l['MeV_per_He'] for l in lines if l['energy_MeV']==.511)
  flight=nuclear_gamma+(anni+sp['main_beta_kinetic_MeV_per_He'])*math.exp(-floor(5.14304)*S)
  soft=sp['cycle']['escaped_prompt_MeV_per_He']*math.exp(-floor(1.)*S)
  hard=sp['cycle']['escaped_prompt_MeV_per_He']*math.exp(-floor(20.)*S)
  boosts=[]
  for D in [1.,maxD,1.10]:
   delayed=sum(D*l['MeV_per_He']*math.exp(-floor(D*l['energy_MeV'])*S) for l in lines if l['energy_MeV']!=.511)+D*(anni+sp['main_beta_kinetic_MeV_per_He'])*math.exp(-floor(D*5.14304)*S)
   boosts.append({'source_Doppler_factor':D,'delayed_hard_reserve_loss_fraction_thermal':delayed/thermal,'plus_saved_prompt_at20MeV_lab_loss_fraction_thermal':(delayed+hard)/thermal,'plus_entire_D1p10_hot_direct_bypass_loss_fraction_thermal':(delayed+hard+hotreserve)/thermal,'remaining_nonphoton_after_full_hot_reserve_fraction_for_99pct_capture':.01-(delayed+hard+hotreserve)/thermal,'remaining_nonphoton_loss_fraction_for_99pct_capture':.01-(delayed+hard)/thermal,'remaining_nonphoton_loss_MeV_per_He_for_99pct_capture':.01*thermal-delayed-hard})
  guards=[{'removed_boundary_column_kg_m2':G,'delayed_primary_plus_beta_leakage_MeV_per_He':sum(l['MeV_per_He']*math.exp(-l['minimum_degraded_energy_mu_en_m2_kg']*(S-G)) for l in lines)+sp['main_beta_kinetic_MeV_per_He']*math.exp(-floor(5.14304)*(S-G))} for G in [0,20,50]]
  aggregate.append({'sigma_kg_m2':S,'Doppler_boost_scenarios':boosts,'primary_delayed_photon_leakage_upper_MeV_per_He':g,'primary_delayed_photon_leakage_upper_fraction_thermal':g/thermal,'all_beta_as_hard_photon_leakage_upper_MeV_per_He':b,'delayed_primary_plus_all_beta_leakage_upper_fraction_thermal':(g+b)/thermal,'inflight_annihilation_all_annihilation_plus_beta_below5p143MeV_fraction_thermal':flight/thermal,'inflight_plus_all_grey_prompt_below1MeV_fraction_thermal':(flight+soft)/thermal,'inflight_plus_all_grey_prompt_at20MeV_fraction_thermal':(flight+hard)/thermal,'plus_all_grey_prompt_below1MeV_fraction_thermal':(g+b+soft)/thermal,'plus_all_grey_prompt_at20MeV_fraction_thermal':(g+b+hard)/thermal,'boundary_column_sensitivity_NOT_electron_transport_bound':guards})
 mc=[transport(E,S,50000,190601+int(S)+int(E*1000)) for S in SIGMAS for E in [.511,2.312593]]
 return {'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'hot_bypass_sha256':hashlib.sha256(Path(__file__).with_name('hot-bypass.json').read_bytes()).hexdigest(),'hot_direct_bypass_full_loss_reserve_MeV_per_He':hotreserve,'spectrum_sha256':hashlib.sha256(Path(__file__).with_name('spectrum.json').read_bytes()).hexdigest(),'NIST_source':'https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z06.html','NIST_definition_source':'https://physics.nist.gov/PhysRefData/XrayMassCoef/chap3.html','KN_sampling_source':'https://geant4.web.cern.ch/documentation/dev/prm_html/PhysicsReferenceManual/electromagnetic/gamma_incident/compton/compton.html','nist_table_units':['MeV','cm2/g','cm2/g'],'nist_table':NIST.tolist(),'primary_lines':lines,'thermal_boundary_MeV_per_He':thermal,'MC_thermal_weighted_primary_leakage_fraction':[{'sigma_kg_m2':S,'primary_two_line_fraction_thermal':sum(next(l['MeV_per_He'] for l in lines if l['energy_MeV']==m['incident_energy_MeV'])*m['combined_leakage_fraction'] for m in mc if m['sigma_kg_m2']==S)/thermal,'sampling_standard_error_thermal':math.sqrt(sum((next(l['MeV_per_He'] for l in lines if l['energy_MeV']==m['incident_energy_MeV'])*m['sampling_standard_error'])**2 for m in mc if m['sigma_kg_m2']==S))/thermal} for S in SIGMAS],'aggregate_localized_cascade_bounds':aggregate,'KN_only_transport_sensitivity':mc,'bound_assumptions':['Doppler scenarios boost both delayed-photon energy and its total source energy. The actual maximum resolved final bulk factor is recorded;1.10 is an additional finite envelope, not a mathematically absolute thermal-particle speed limit. Saved grey-prompt energy is already a model escape-energy account and its20MeV scenario is specified in the laboratory frame.','Piecewise log-linear NIST coefficient interpolation, with no smaller energy-absorption coefficient below1keV than the stated minimum; low-energy photons ultimately absorbed rather than upgraded back into hard photons.','Local charged-secondary energy deposition and photon regeneration at the collision location. Nonlocal electron and positron transport is not solved.','No photon energy upscatter by a hot material. Receiver at~1000K, apart from transient heating, versus keV/MeV photons.','Closed shell returns inner-directed photons to another receiver patch; no ports/gaps and no credit for heat received in cold components.','Energy-weighted exponential transport bound uses the minimum mu_en over all down-scattered photon energies, not incident mu_en.','Removing20/50kg/m2 from the bound illustrates sensitivity to an assumed boundary layer; it is not a rigorous replacement for coupled electron-photon transport.'],'MC_limitations':['Compton-only free cold electron model deposits recoil energy locally. Missing photoelectric absorption, coherent scattering, pair production, charged-particle ranges and bremsstrahlung prevent treating these results as an actual carbon capture certification.','Untracked<1keV photon energy is charged entirely to leakage. Exact energy balance checked for every simulated history.','Sampling intervals describe Monte Carlo noise only; they are not engineering credible intervals.','Incident MC line energies are rest-frame sensitivity values; the analytic tables separately include Doppler scenarios. Normal incidence plus cavity returns is a specified geometry; delayed photons created in the wall itself have a different source-depth problem.']}
if __name__=='__main__':
 d=main();Path(__file__).with_name('containment.json').write_text(json.dumps(d,indent=2)+'\n');print(json.dumps({'bounds':d['aggregate_localized_cascade_bounds'],'MC':d['KN_only_transport_sensitivity']},indent=2))
