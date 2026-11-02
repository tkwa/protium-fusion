"""Independent cold-tail ReacLib audit; no hydro or historical-file changes."""
import argparse,json,hashlib,sys,math
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--trace',type=Path,required=True);p.add_argument('--record',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.root=a.root.resolve();sys.path.insert(0,str(a.root))
import pulse_model as pm
r=json.loads(a.record.read_text());tr=np.load(a.trace);nz=int(tr['zones']);nv=int(tr['variables_per_zone']);nn=nv-5
z=tr['state'].reshape(nz,nv,-1);T=z[:,nn,:]*r['T0_K'];Y=z[:,:nn,:];times=tr['time_s'];w=np.diff(np.linspace(0,1,nz+1)**3)
rad=np.concatenate([np.zeros((1,T.shape[1])),z[:,nn+1,:]],axis=0)
rho=r['rho0_kg_m3']/1000*w[:,None]/np.diff(rad**3,axis=0)
NAv=6.02214076e23;HBAR=1.054571817e-34;KB=1.380649e-23;MU=1.66053906660e-27
nidx=pm.SPECIES.index('n');isn=np.any(pm.REACT==nidx,axis=1);isweak=np.array([bool(x.weak_type) for x in pm.RATES]);two=pm.DEXP==1
regularized=isn&~isweak&two&(pm.Q>0)

def vals(temp):
 t=float(temp)/1e9;tf=np.array([1.,1/t,t**(-1/3),t**(1/3),t,t**(5/3),np.log(t)])
 with np.errstate(over='ignore',under='ignore',invalid='ignore'):
  return np.bincount(pm.RATE_INDEX,weights=np.exp(pm.COEFF@tf),minlength=pm.NR)

def finite(x):return float(x) if np.isfinite(x) else str(float(x))
anchor=vals(1e7);scans=[]
for temp in [1e2,1e3,1e4,1e5,3.0136540238366614e5,1e6,1e7,1e8,1e9,2.8e9]:
 v=vals(temp)
 bad=[]
 for j,rate in enumerate(pm.RATES):
  if v[j]>max(anchor[j]*1e8,1e15):
   bad.append({'index':j,'reaction':str(rate),'coefficient':finite(v[j]),'anchor10MK':finite(anchor[j]),'is_regularized_exothermic_neutron_binary':bool(regularized[j])})
 scans.append({'T_K':temp,'large_extrapolations':bad,'maximum_rate':finite(v.max()),'unfixed_large_extrapolation_count':sum(not b['is_regularized_exothermic_neutron_binary'] for b in bad)})

rates=[]
for j in np.where(regularized)[0]:
 rate=pm.RATES[j];targets=[int(k) for k in pm.REACT[j] if k>=0 and k!=nidx];target=targets[0];A=pm.AA[target];mu=MU*A/(A+1)
 # Conservative full elastic s-wave 4pi/k^2 unitarity envelope (capture's
 # absorptive partial-wave bound is smaller); not an evaluated reaction rate.
 def envelope(temp):return NAv*1e6*4*np.sqrt(2*np.pi)*HBAR**2/(mu**1.5*np.sqrt(KB*temp))
 rev=[k for k in range(pm.NR) if np.array_equal(pm.STOICH[:,k],-pm.STOICH[:,j])]
 rates.append({'index':int(j),'reaction':str(rate),'label':str(getattr(rate,'label','')),'fname':str(rate.fname),'Q_MeV':float(pm.Q[j]),'anchor10MK':float(anchor[j]),'rate_at_301365K':finite(vals(301365.40238366614)[j]),'conservative_swave_unitarity_at301365K':float(envelope(301365.40238366614)),'coefficients':[list(map(float,s.a)) for s in rate.sets],'reverse_indices':rev})

cold=T<1e7;yc=np.maximum(Y,0);negative=np.maximum(-Y,0)
# Conservative instantaneous isotope-energy artifact bound from discarded
# negative baryon fraction, c^2 times fraction / actual input J/kg.
negative_baryons=(negative*pm.AA[None,:,None]).sum(axis=1)
maxnegative=float((w[:,None]*negative_baryons).sum(axis=0).max())
peak_n=float(np.max(yc[:,nidx,:]));cold_peak_n=float(np.max(np.where(cold,yc[:,nidx,:],0)))
weighted_cold_n=np.sum(w[:,None]*np.where(cold,yc[:,nidx,:],0),axis=0)
# Neutron source channels without a neutron reactant, conservatively evaluated
# at 10MK in every cold sample. This overstates barrier-suppressed cold source
# rates; their actual rates are also recorded.
source=(pm.STOICH[nidx]>0)&~isn
actual_source=np.zeros_like(T);anchor_source=np.zeros_like(T)
flux_rows=[]
for i,k in zip(*np.where(cold)):
 real=vals(T[i,k]);counter=[]
 for j in np.where(source)[0]:
  factor=pm.PREF[j]*rho[i,k]**pm.DEXP[j]
  for kk in pm.REACT[j]:
   if kk>=0:factor*=yc[i,kk,k]
  f=float(factor*real[j])*pm.STOICH[nidx,j];fa=float(factor*anchor[j])*pm.STOICH[nidx,j]
  actual_source[i,k]+=f;anchor_source[i,k]+=fa
  if fa>0:counter.append((fa,int(j),f))
 if counter and k==T.shape[1]-1:
  flux_rows.append({'zone':int(i),'T_K':float(T[i,k]),'rho_g_cm3':float(rho[i,k]),'source_channels':[{'index':j,'reaction':str(pm.RATES[j]),'actual_dY_n_s':f,'at10MK_dY_n_s':fa} for fa,j,f in sorted(counter,reverse=True)[:5]]})
# Use cold duration times the largest anchor source, more conservative than
# integrating a sparse trace over just cold intervals.
source_budget=times[-1]*float(np.max((w[:,None]*anchor_source).sum(axis=0)))
initial_cold_budget=float(np.sum(w*np.max(np.where(cold,yc[:,nidx,:],0),axis=1)))
# Ten MeV per neutron is a diagnostic capture-energy bound, not a universal
# cascade theorem; inventory values allow other energy-per-neutron choices.
Einperkg=r['initial_energy_J']/r['mass_kg'];MeV=1.602176634e-13
inventory=(initial_cold_budget+source_budget)
atol=float(r.get('atol',1e-11));cold_weight=float(w[np.any(cold,axis=1)].sum());atol_neutron_budget=atol*cold_weight

summary={'trace_samples':len(times),'cold_cells_ever':int(np.any(cold,axis=1).sum()),'cold_zone_ids':np.where(np.any(cold,axis=1))[0].tolist(),'min_temperature_K':float(T.min()),'cold_neutron_max_Y':cold_peak_n,'all_temperature_neutron_max_Y':peak_n,'maximum_mass_weighted_cold_neutron_Y':float(weighted_cold_n.max()),'sum_percell_cold_neutron_max_Y_weighted':initial_cold_budget,'upper_source_neutrons_Y_over_entire_saved_duration_at10MK':source_budget,'ten_MeV_per_available_cold_neutron_energy_over_initial':inventory*10*MeV/MU/Einperkg,'cold_material_fraction_ever':cold_weight,'cold_neutron_inventory_at_one_absolute_tolerance_Y':atol_neutron_budget,'ten_MeV_per_atol_scale_cold_neutron_energy_over_initial':atol_neutron_budget*10*MeV/MU/Einperkg,'abundance_error_note':'The very small positive Y from the saved trace is below atol and cannot establish physical precision. The atol-scale calculation is a sensitivity scale, not a rigorous global ODE error bound.','maximum_negative_baryon_fraction_massweighted':maxnegative,'c2_negative_inventory_energy_over_initial_bound':maxnegative*(299792458.**2)/Einperkg,'maximum_positive_cold_neutron_source_dY_s_at_actual_T':float(actual_source.max()),'maximum_positive_cold_neutron_source_dY_s_at10MK':float(anchor_source.max())}
# Whole-cell RHS invariant identities hold for any common set of reaction flows.
stoich_b=float(np.max(np.abs(pm.AA@pm.STOICH)));energy=float(np.max(np.abs(pm.Q+pm.MASS@pm.STOICH)))
validation={'nuclear_stoichiometry_baryon_identity':stoich_b==0.,'reaction_Q_mass_identity':energy<1e-12,'source_budget_finite':bool(np.isfinite(source_budget)),'all_selected_regularizations_are_binary_exothermic_neutron_channels':bool(np.all(isn[regularized]&two[regularized]&(pm.Q[regularized]>0)))}
assert all(validation.values())
out={'scope':'Diagnostic contingency only: first correct the isotropic artificial-viscosity activation and retain original rates if the corrected trajectory stays hot. Cold nuclear-rate audit only. The anchor continuation is a physically motivated approximation, not an evaluated low-energy rate. No completed hydrodynamics is claimed.','summary':summary,'regularized_reactions':rates,'temperature_scans':scans,'neutron_source_at_cold_endpoint':flux_rows,'validation':validation,'input_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [a.trace,a.record,a.root/'pulse_model.py',a.root/'radial_pulse.py']},'proposed_rule':{'threshold_K':1e7,'classification':'binary exothermic neutron-induced reactions only','below_threshold':'constant N_A<sigma v> equal to10MK value, matching s-wave1/v cross-section asymptote','weak_charged_and_reverse_reactions':'retain actual-temperature evaluation except diagnosed C12+O16 charged-fit overflow: if needed below10MK, use a matched Gamow-suppressed continuation. Reverse-pair correction must be consistent if its cold flow is material','sensitivities':'zero cold neutron-induced rates versus capped upper-rate stress, with all hot coefficients unchanged','physical_status':'Approximate asymptotic continuation and bounded trace-inventory sensitivity, not a measured rate replacement.'}}
a.output.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps(summary,indent=2));print('regularized',len(rates),'checks',validation);print('remaining divergences',[(s['T_K'],s['unfixed_large_extrapolation_count']) for s in scans])
