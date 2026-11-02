"""Collection checkpoint: source-resolved CSDA and nuclear-loss requirements.
No historical input changed. Constant-cross-section rows are acceptance
sensitivity tests, not predictions or rigorous upper bounds. Only the C+C
fusion model is tied to a measured nuclear fit; other species need closure.
"""
import argparse,hashlib,importlib.util,json,math
from pathlib import Path
import numpy as np
from scipy.integrate import cumulative_trapezoid
HERE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def quant(x,w,ps=(.5,.9,.99,.999,.9999)):
 ix=np.argsort(x);x=x[ix];w=w[ix];s=np.cumsum(w)
 return {str(p):float(x[min(np.searchsorted(s,p*s[-1]),len(x)-1)]) for p in ps}
def main(root,out):
 base=root/'docs/protium-fusion';path=base/'membrane-receiver/stopping.py';spec=importlib.util.spec_from_file_location('component_stopping',path);st=importlib.util.module_from_spec(spec);spec.loader.exec_module(st)
 extra=HERE/'collection-stopping-data.json';st.DATA+=json.loads(extra.read_text())['selected_data']
 st.EL.update({'Al':(13,26.9815385,166.),'Sn':(50,118.710,488.)})
 original_analytic=st.analytic
 def positive_analytic(eu,z,m):
  if m!='tin':return original_analytic(eu,z,m)
  eu=np.asarray(eu);anchor=.5
  return np.where(eu>=anchor,st.bethe(np.maximum(eu,anchor),z,m),float(st.bethe(anchor,z,m))*np.sqrt(eu/anchor))
 st.analytic=positive_analytic
 st.MAT.update({'aluminium':{'rho':2300.,'parts':{'Al':1},'db':'Al','I':166.},'tin':{'rho':6700.,'parts':{'Sn':1},'db':'Sn','I':488.},'alumina':{'rho':3900.,'parts':{'Al':2,'O':3},'db':'Al2O3','I':145.2},'hydrogen':{'rho':.001,'parts':{'H':1},'db':'H2','I':19.2}})
 wp=base/'membrane-receiver/pulse-loading.json';cp=base/'open-loop-redesign/station-inputs.json';w=json.loads(wp.read_text());cy=json.loads(cp.read_text())['cycle'];A=np.array(w['A']);Z=np.array(w['Z']);area=4*np.pi*w['radius_m']**2
 rows=[];coverage=[]
 for material in ['carbon','aluminium','alumina','tin','hydrogen']:
  at=st.target(material)[0]*st.MU/st.target(material)[2]
  for ion in ['C','N','O']:
   rr=[r for r in st.DATA if r['target_name']==st.MAT[material]['db'] and r['projectile_name']==ion]
   coverage.append({'material':material,'projectile':ion,'measured_MeV_u_range':[min(float(r['energy']) for r in rr),max(float(r['energy']) for r in rr)] if rr else None,'primary_DOIs':sorted({r['doi'] for r in rr if r['doi']})})
  for stage,s in w['stages'].items():
   p=s['spectral_packets'];n=np.array(p['mass_fluence_kg_m2'])[:,None]/st.MU*np.array(p['Y_at_wall']);eu=np.array(p['kinetic_MeV_per_u']);den=n.sum();cno=n[:,(Z>=6)&(Z<=8)].sum();cat=next(x for x in cy['stage_pulse_data'].values() if x['stage']==stage);perHe=area*cat['pulses_per_He'];species=[];all_range=[];all_n=[];cno_range=[];cno_n=[];coating_loss={sig:0. for sig in [.1,1.,3.]};cno_loss=coating_loss.copy();cfusion=0.
   for j in range(len(A)):
    if Z[j]==0 or n[:,j].sum()<den*1e-10:continue
    eg,rg,dg,se,sn=st.table(material,int(Z[j]),int(A[j]));assert np.all(se+sn>0) and np.all(np.diff(rg)>=0),(material,j)
    ranges=np.interp(A[j]*eu,eg,rg)/100
    all_range.extend(ranges);all_n.extend(n[:,j]);iscno=6<=Z[j]<=8
    if iscno:cno_range.extend(ranges);cno_n.extend(n[:,j])
    pp={sig:float(np.sum(n[:,j]*(-np.expm1(-sig*1e-28*ranges/at)))) for sig in coating_loss}
    for sig in pp:
     coating_loss[sig]+=pp[sig]
     if iscno:cno_loss[sig]+=pp[sig]
    fc=None
    if material=='carbon' and Z[j]==6:
     ec=eg*12/(12+A[j]);ew=np.maximum(ec,4.9);sig=2.1*7.4**2/(2*ew)*np.logaddexp(0,2*np.pi*(ew-5.7)/2.1)/100
     g=.9895*6*6*np.sqrt(12*A[j]/(12+A[j]));sig*=np.where(ec<4.9,4.9/ec*np.exp(np.maximum(-700.,-g*(1/np.sqrt(ec)-1/np.sqrt(4.9)))),1.)
     tau=cumulative_trapezoid(sig/(se+sn)/100/at*1e-28,eg,initial=0);fc=float(np.sum(n[:,j]*(-np.expm1(-np.interp(A[j]*eu,eg,tau)))));cfusion+=fc
    species.append({'species':w['species'][j],'nuclei_per_He':float(n[:,j].sum()*perHe),'mean_stopping_column_kg_m2':float(n[:,j]@ranges/n[:,j].sum()),'max_stopping_column_kg_m2':float(max(ranges)),'one_barn_reaction_probability':pp[1.]/n[:,j].sum(),'C_fusion_fit_reactions_per_He':None if fc is None else fc*perHe})
   cr=np.array(cno_range);cn=np.array(cno_n);rho=st.MAT[material]['rho']
   rows.append({'material':material,'stage':stage,'nuclei_per_He':float(den*perHe),'CNO_nuclei_per_He':float(cno*perHe),'kinetic_fluence_J_m2':float(np.sum(n*A*eu[:,None])*st.EV*1e6),'max_source_MeV_u':float(eu.max()),'CNO_range_quantiles_kg_m2':quant(cr,cn),'CNO_range_quantiles_um_at_assumed_density':{q:x/rho*1e6 for q,x in quant(cr,cn).items()},'all_nuclei_range_quantiles_kg_m2':quant(np.array(all_range),np.array(all_n)),'assumed_density_kg_m3':rho,'constant_cross_section_acceptance_tests':[{'sigma_barn':sig,'target_reactions_per_He':coating_loss[sig]*perHe,'CNO_projectile_reaction_probability':cno_loss[sig]/cno,'scope':'If each selected nonelastic event destroys one target atom and one incidentCNO, these are loss sensitivities. Actual channels can retain catalyst; this is not a cross-section prediction.'} for sig in coating_loss],'C_fusion_fit_reactions_per_He':cfusion*perHe,'species':species})
 summaries=[]
 wall_budget=cy['net_structural_carbon_allowance_kg_per_He']/(12.011*st.MU)
 for mat in ['carbon','aluminium','alumina','tin','hydrogen']:
  r=[x for x in rows if x['material']==mat];cno=sum(x['CNO_nuclei_per_He'] for x in r);one=sum(x['constant_cross_section_acceptance_tests'][1]['target_reactions_per_He'] for x in r)
  summaries.append({'material':mat,'all_projectiles_per_He':sum(x['nuclei_per_He'] for x in r),'CNO_per_He':cno,'one_barn_target_events_per_He':one,'CNO_one_barn_reaction_fraction':sum(x['CNO_nuclei_per_He']*x['constant_cross_section_acceptance_tests'][1]['CNO_projectile_reaction_probability'] for x in r)/cno,'carbon_full_budget_barn_equivalent':wall_budget/one if mat=='carbon' else None,'carbon_fusion_fit_events_per_He':sum(x['C_fusion_fit_reactions_per_He'] for x in r)})
 station=json.loads(cp.read_text());rate=.8/(station['receiver']['occupied_s_per_He_at_R1e9']['mass999']*w['radius_m']/1e9);year=365.25*86400;gross=rate*cy['gross_fuel_kg_per_He'];batch=[]
 for material,thick in [('carbon',10.),('carbon',25.),('alumina',50.),('aluminium',100.)]:
  sigma=st.MAT[material]['rho']*thick*1e-6
  batch.append({'material':material,'thickness_um':thick,'whole_skin_kg_m2':sigma,'annual_processing_kg_s':sigma*area/year,'annual_processing_W_m2_at2MJkg':sigma*2e6/year,'fuel_fraction_after_one_year_if_uniformly_mixed':gross*year/area/sigma,'scope':'Annualwhole-skinprocessing sensitivity; no everypulse release requirement. Chemicalregenerationenergy and yield need qualification.'})
 result={'scope':'First component checkpoint; material-stopping depth and loss acceptance screen, not achievedreturn or nucleartransportclosure. Solid-state stopping used for molten-metal depth; density/phase/charge and straggling uncertainties require qualification.','source_hashes':{('docs/protium-fusion/component-design/'+p.name if p==extra else str(p.relative_to(root))):sha(p) for p in [path,wp,cp,extra,base/'membrane-receiver/stopping-data.json']},'code_sha256':sha(Path(__file__)),'source':'https://www-nds.iaea.org/stopping/','C_fusion_model':{'primary':'https://doi.org/10.1103/PhysRevC.85.014607','Rb_fm':7.4,'Vb_MeV':5.7,'hbaromega_MeV':2.1,'measured_fit_Ecm_MeV':[4.9,8.0],'limitations':'C12+C13 measured fit appliedtoC12andC13 onnaturalcarbon asillustrativeenvelope; matchedGamowtail below4.9 and Wongextension above8 are labeledextrapolations. NoN/Omodelborrowed. Allcompoundfusion countedlost; nosecondarytransport or transferchannels.'},'budget':{'CNO_return_requirement':1-cy['carbon_loss_per_cold_stage_output'],'total_C_H':.007,'fuel_C_H':.0056,'remaining_material_C_H':.0014,'carbon_atoms_all_material_budget_per_He':wall_budget},'batch_extraction':{'helium_per_s':rate,'gross_fuel_kg_s':gross,'gross_fuel_kg_m2_year':gross*year/area,'synchronized_one_year_fuel_inventory_upper_kg':gross*year,'staggered_uniform_age_mean_fuel_inventory_kg':.5*gross*year,'scenarios':batch},'coverage':coverage,'rows':rows,'summary':summaries,'checks':{'all_ranges_positive':all(x['CNO_range_quantiles_kg_m2']['0.5']>0 for x in rows),'all_probabilities_0_1':all(0<=z['CNO_projectile_reaction_probability']<=1 for x in rows for z in x['constant_cross_section_acceptance_tests'])}}
 out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps({'summary':summaries,'depths':[{'material':x['material'],'stage':x['stage'],'q':x['CNO_range_quantiles_um_at_assumed_density']} for x in rows],'checks':result['checks']},indent=2))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--output',type=Path,default=HERE/'collection.json');a=ap.parse_args();main(a.root,a.output)
