"""Convolve primary low-energy BCA coefficients and explicit high-energy sensitivities."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
P=Path(__file__).resolve().parent;S=P.parent
sys.path.insert(0,str(S));import calculate as stopping;import spectra as sp
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
mu=1.66053906660e-27;e=1.602176634e-19

def loginterp(x,xp,yp):return np.exp(np.interp(np.log(x),np.log(xp),np.log(yp)))
def normal(rows):return sorted([r for r in rows if r['angle_deg']==0],key=lambda r:r['E_eV'])
def cdf_bracket(row,cut_eV):
 q=np.array(row['energy_fraction_quantiles_10_30_50_70_90'])*row['E_eV'];p=np.array([.1,.3,.5,.7,.9]);j=int(np.searchsorted(q,cut_eV));lo=0. if j==0 else float(p[j-1]);hi=1. if j==5 else float(p[j])
 if 'energy_fraction_bounds'in row:
  mn,mx=np.array(row['energy_fraction_bounds'])*row['E_eV']
  if cut_eV<mn:lo=hi=0.
  if cut_eV>=mx:lo=hi=1.
 return lo,hi
def main():
 data=json.loads((P/'sources/primary-tables.json').read_text());single=json.loads((P/'single-collision.json').read_text());surface=json.loads((S/'surface-spectra.json').read_text());hext=json.loads((P/'sources/eckstein-He-normal.json').read_text())
 _,_,old,_,_,_=sp.setup();tables=stopping.Tables(old);rows=[];checks=[]
 for stage,x in surface['stages'].items():
  for species,A,Z in [('p',1,1),('he4',4,2)]:
   E=np.array(x['packet_energy_MeV_u'])*A*1e6;n=np.array(x['species'][species]['packet_nuclei_m2']);N=n.sum();ker=single['kernel'][species];ke=np.array([r['E_MeV']*1e6 for r in ker]);kr=np.array([r['R']for r in ker]);kre=np.array([r['RE']for r in ker]);t=tables.get(Z,A);rnuc=np.interp(E/1e6,t['E'],t['R'])/100
   # A column-times-cross-section sensitivity; every event is counted as escaping.
   # No claim that 1 or 10 barn bounds the physical differential cross section.
   col_per_barn=rnuc/(12.011*mu)*1e-28
   row={'stage':stage,'species':species,'nuclei_m2':float(N),'mean_incident_flux_m2_s':float(N/x['period_s']),'unknown_incident_below1eV_fraction':float(n[E<1].sum()/N),'incident_KE_J_m2':float(n@E)*e,'column_event_probability_per_barn':float(n@col_per_barn/N),'scenarios':{}}
   for dataset in ['old','new']:
    fn=('H_on_C' if species=='p' else 'He_on_C')+('' if dataset=='old' else ('_1.00' if species=='p' else '_0.00'))+'.txt';p=normal(data['datasets'][fn]['rows']);pe=np.array([r['E_eV']for r in p]);pr=np.array([r['RN']for r in p]);pre=np.array([r['RE']for r in p]);low=(E>=1)&(E<=5000);extension=(E>5000)&(E<=20000) if species=='he4' else np.zeros_like(low);high=(E>20000) if species=='he4' else E>5000
    vals=np.zeros_like(E);evals=np.zeros_like(E);vals[low]=loginterp(E[low],pe,pr);evals[low]=loginterp(E[low],pe,pre)*E[low]
    if species=='he4':
     ex=np.array(hext['rows']);vals[extension]=loginterp(E[extension],ex[:,0],ex[:,1]);evals[extension]=loginterp(E[extension],ex[:,0],ex[:,2])*E[extension]
    rv=loginterp(E[high],ke,kr);ev=loginterp(E[high],ke,kre)*E[high]
    for factor in [1.,2.,3.,10.]:
     R=vals.copy();RE=evals.copy();R[high]=np.minimum(1,factor*rv);RE[high]=np.minimum(E[high],factor*ev);key=f'{dataset}_primary_high_x{factor:g}';r=float(n@R/N);en=float(n@RE)*e
     row['scenarios'][key]={'R':r,'low_primary_R':float(n[low]@vals[low]/N),'He_5to20keV_primary_R':float(n[extension]@vals[extension]/N),'high_model_R':float(n[high]@R[high]/N),'mean_first_reflected_E_eV':en/(N*r*e),'reflected_KE_J_m2':en,'all_reflected_eventually_warm_fraction_ceiling':r+row['unknown_incident_below1eV_fraction']}
     for cut,modelkey in [(1.,'Rout_below_1e-06MeV'),(100.,'Rout_below_0.0001MeV')]:
      bracket=np.array([cdf_bracket(v,cut)for v in p]);bound=np.array([float(n[low]@np.interp(np.log(E[low]),np.log(pe),pr*bracket[:,j])/N)for j in [0,1]])
      # No energy quantiles in He5–20keV extension: its entire RN is only an upper allowance.
      bound[1]+=float(n[extension]@vals[extension]/N)+row['unknown_incident_below1eV_fraction']
      highcut=float(n[high]@np.interp(np.log(E[high]),np.log(ke),[v[modelkey]for v in ker])/N)*factor
      row['scenarios'][key][f'exit_below_{cut:g}eV_fraction_bracket']=(bound+highcut).tolist()
     checks.extend([0<=r<=1,en>=0,en<=row['incident_KE_J_m2']])
   rows.append(row)
 means={}
 for species in ['p','he4']:
  selected=[r for r in rows if r['species']==species];flux=sum(r['mean_incident_flux_m2_s']for r in selected)
  mean={k:sum(r[k]*r['mean_incident_flux_m2_s']for r in selected)/flux for k in ['unknown_incident_below1eV_fraction','column_event_probability_per_barn']};mean['mean_incident_flux_m2_s']=flux;mean['scenarios']={}
  for key in rows[0]['scenarios']:
   mean['scenarios'][key]={k:sum(r['scenarios'][key][k]*r['mean_incident_flux_m2_s']for r in selected)/flux for k in ['R','low_primary_R','He_5to20keV_primary_R','high_model_R','all_reflected_eventually_warm_fraction_ceiling']}
   q=sum(r['scenarios'][key]['reflected_KE_J_m2']/surface['stages'][r['stage']]['period_s']for r in selected);mean['scenarios'][key]['mean_reflected_power_W_m2']=q;mean['scenarios'][key]['mean_first_reflected_E_eV']=q/(flux*mean['scenarios'][key]['R']*e)
   for cut in [1,100]:
    k=f'exit_below_{cut}eV_fraction_bracket';mean['scenarios'][key][k]=[sum(r['scenarios'][key][k][j]*r['mean_incident_flux_m2_s']for r in selected)/flux for j in [0,1]]
  baseline=.00103 if species=='p' else .000875;tau=45253465. if species=='p' else 63769871.8
  mean['single_collision_first_flight_inventory_relative_to_baseline_warm_inventory']=single['cadence_weighted'][species]['weighted_first_flight_s']/(baseline*tau)
  mean['effective_backward_cross_section_barn_to_add_0p001_if_all_column_events_escape']=.001/mean['column_event_probability_per_barn']
  means[species]=mean
 comparisons={}
 for species,fn in [('p','H_on_C_1.00.txt'),('he4','He_on_C_0.00.txt')]:
  p=normal(data['datasets'][fn]['rows']);ke=np.array([r['E_MeV']*1e6 for r in single['kernel'][species]]);kr=np.array([r['R']for r in single['kernel'][species]]);comparisons[species]=[{'E_eV':r['E_eV'],'primary_R':r['RN'],'model_R':float(loginterp(r['E_eV'],ke,kr)),'primary_over_model':float(r['RN']/loginterp(r['E_eV'],ke,kr))}for r in p if r['E_eV']>=1000]
 angles={name:[{'E_eV':r['E_eV'],'angle_deg':r['angle_deg'],'RN':r['RN'],'mean_reflected_energy_eV':r['E_eV']*r['RE']/r['RN']}for r in data['datasets'][name]['rows']if r['E_eV']in[1000,5000]]for name in ['H_on_C_1.00.txt','He_on_C_0.00.txt']}
 out={'scope':__doc__,'inputs_sha256':{str(p):sha(p)for p in [Path(__file__),P/'sources/primary-tables.json',P/'sources/eckstein-He-normal.json',P/'single-collision.json',S/'surface-spectra.json',S/'calculate.py']},'stages':rows,'cadence_weighted':means,'primary_vs_single_collision':comparisons,'angle_diagnostic':angles,'method':['Log-log interpolation of primary RN/RE inside actual table domain. Native old/new TRIM split retained. He5–20keV uses the separate Eckstein2002 table; its5keV normalization differs modestly from newer EIRENE.','Above5keV H /20keV He, use the explicit single-collision Coulomb model times1,2,3,10. These are model sensitivities, not fitted posterior or confidence limits. No source table is extrapolated silently.','Incident below1eV is carried as zero-to-one reflection uncertainty. Multiplicity-free ultimate warm fraction from first-reflected incident atoms cannot exceed first reflection, even with arbitrarily many wall encounters.','Column-event sensitivity counts all collisions across the full stopping column as escaping. It is intentionally generous for an assumed integrated backward cross section. It neither validates a nuclear cross-section cap nor includes sputtered pre-existing target gas.','Normal-incidence family assumes a smooth macroscopic spherical surface and no trajectory-deflecting field. Angle tables diagnose local facet sensitivity; no roughness distribution is inferred.'],'checks':{'number':len(checks),'passed':sum(checks)}}
 (P/'comparison.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps({'means':means,'checks':out['checks'],'comparison':comparisons},indent=2));assert all(checks)
if __name__=='__main__':main()
