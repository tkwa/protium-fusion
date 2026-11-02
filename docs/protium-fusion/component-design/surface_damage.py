"""Bounded selected64 graphite skin damage, heat resistance and recycled restoration.
No hydro rerun. NRT uses total ZBL elastic transfer as a damage-energy PROXY;
it is neither surviving defects nor a conductivity/lifetime law.
"""
from pathlib import Path
import argparse, hashlib, importlib.util, json, math
import numpy as np
YEAR=365.25*86400.; EV=1.602176634e-19
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load_module(p,name):
 s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--plant',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--receivers',type=Path);a=ap.parse_args()
 root=a.root.resolve();base=root/'docs/protium-fusion';receiver_dir=a.receivers or base/'mixture-design/receiver-results';plant=json.loads(a.plant.read_text())['top_detailed'][0]['plant']
 radius=plant['radius_m'];area=4*math.pi*radius**2;rho=2000.;cp=1800.;kbulk=50.;assert abs(area/plant['area_m2']-1)<1e-12
 code=base/'mixture-design/receiver_stopping.py';rs=load_module(code,'damage_stopping');pm,st,decay,channels,hashes=rs.setup(root)
 hashes.update({str(code):sha(code),str(a.plant):sha(a.plant)})
 from fast_eos import FastEOS
 eos=FastEOS();hashes[str(base/'fast_eos.py')]=sha(base/'fast_eos.py');hashes[str(base/'electron-eos.npz')]=sha(base/'electron-eos.npz')
 edges=np.array([0.,1.,5.,10.,25.,50.,100.,250.,1000.]); factors=[.5,1.,2.];eds=[12.,25.,40.,60.];stages=[];cache={};stage_cat={'A':'n14','B1':'n15','B2':'c12','C13':'c13'}
 for stage,case in plant['selected_cases'].items():
  cpath=receiver_dir/f'{case}.receiver.json';rec=json.loads(cpath.read_text());hashes[str(cpath)]=sha(cpath)
  record=Path(rec['source_paths']['record']);d,n,eu,ncat,checks,hh,vmax=rs.spectrum(record,pm,decay,radius,64)
  hashes[str(record)]=hh['record'];hashes[d['trace_path']]=hh['trace'];period=plant['pulse_interval_s'][stage_cat[stage]]
  qpeak=rec['ballistic']['peak_kinetic_flux_W_m2']*(rec['radius_m']/radius)**3
  tr=np.load(d['trace_path']);nz=int(tr['zones']);nv=int(tr['variables_per_zone']);names=d['state_variables'];last=tr['state'].reshape(nz,nv,-1)[:,:,-1];mass=tr['cell_mass_kg']
  ro=np.r_[0.,last[:,names.index('r/R0')]*d['R0_m']];vo=np.r_[0.,last[:,names.index('v/cs0')]*float(tr['velocity_unit_cm_s'])/100];y=np.maximum(last[:,:pm.NN],0.);bar=y@pm.AA;rho_cell=mass/(4*math.pi/3*np.diff(ro**3))
  u=eos.state(rho_cell/1000,last[:,names.index('T/T0')]*d['T0_K'],(y@pm.ZZ)/bar,y.sum(axis=1)/bar).e*1e-4
  assert abs(float(mass@u)/d['final_internal_J']-1)<1e-6
  b=radius-ro[:-1]+np.diff(ro)/np.diff(vo)*vo[:-1];qfull=float(np.max(mass/np.diff(vo)*vo[1:]**2/b/area*(.5*vo[1:]**2+u)))
  # This is a maximum of per-cell arrival fluxes, matching the monotone ballistic map.
  assert qfull>=qpeak*(1-1e-5)
  active=[j for j in range(pm.NN) if pm.ZZ[j]>0 and n[:,j].sum()>n.sum()*1e-12]
  rows=[]
  for factor in factors:
   cum=np.zeros(len(edges));total=0.;species={};skipped=0.;kin=0.
   for j in range(pm.NN):
    en=pm.AA[j]*eu;count=n[:,j]/area;kin+=float(count@en*pm.MEV)
    if j not in active: skipped+=float(count@en*pm.MEV);continue
    key=('carbon',int(pm.ZZ[j]),int(pm.AA[j]),factor)
    if key not in cache: cache[key]=st.table(*key)
    eg,rg,D,se,sn=cache[key];assert en.max()<eg[-1]
    rr=np.interp(en,eg,rg);dn=np.interp(en,eg,D);dep=float(count@dn*pm.MEV);total+=dep
    species[pm.SPECIES[j]]={'incident_number_m2':float(count.sum()),'elastic_transfer_J_m2':dep}
    for i,z in enumerate(edges):
     rem=np.maximum(rr-rho*z*1e-6*100,0);exitE=np.where(rem==0,0,np.interp(rem,rg,eg));delta=dn-np.interp(exitE,eg,D)
     cum[i]+=float(count@delta*pm.MEV)
   assert np.min(np.diff(cum))>-1e-8 and abs(cum[0])<1e-8 and total<=kin
   bins=[]
   for lo,hi,e in zip(edges[:-1],edges[1:],np.diff(cum)):
    atoms=rho*(hi-lo)*1e-6/(12.011*pm.MU)
    bins.append({'depth_um':[lo,hi],'elastic_transfer_J_m2':float(e),'NRT_proxy_dpa_per_pulse':{str(ed):float(.4*e/EV/ed/atoms) for ed in eds}})
   skins=[]
   for thick in [50.,100.]:
    e=float(cum[list(edges).index(thick)]);atoms=rho*thick*1e-6/(12.011*pm.MU)
    skins.append({'thickness_um':thick,'elastic_transfer_J_m2':e,'fraction_elastic_transfer_in_skin':e/total,'NRT_proxy_dpa_per_pulse':{str(ed):.4*e/EV/ed/atoms for ed in eds}})
   rows.append({'electronic_stopping_multiplier':factor,'elastic_transfer_J_m2':total,'elastic_fraction_of_kinetic':total/kin,'unmodeled_kinetic_fraction':skipped/kin,'species':species,'depth_bins':bins,'skins':skins})
  stages.append({'stage':stage,'case':case,'period_s':period,'pulses_per_year':YEAR/period,'CNO_fluence_m2':float(n[:,(pm.ZZ>=6)&(pm.ZZ<=8)].sum()/area),'all_nuclei_fluence_m2':float(n.sum()/area),'peak_kinetic_W_m2':qpeak,'peak_kinetic_plus_advected_internal_W_m2':qfull,'prompt_fluence_J_m2':d['radiation_escaped_J']/area,'energy_1_to_99_arrival_duration_s':(rec['ballistic']['arrival_1_50_99_percent_energy_s'][2]-rec['ballistic']['arrival_1_50_99_percent_energy_s'][0])*radius/rec['radius_m'],'rows':rows,'source_checks':checks})
 annual=[]
 for factor in factors:
  rr=[next(r for r in s['rows'] if r['electronic_stopping_multiplier']==factor) for s in stages]
  annual.append({'electronic_stopping_multiplier':factor,'elastic_transfer_J_m2_year':sum(r['elastic_transfer_J_m2']*s['pulses_per_year'] for r,s in zip(rr,stages)),
   'skins':[{'thickness_um':thick,'NRT_proxy_dpa_year':{str(ed):sum(next(x for x in r['skins'] if x['thickness_um']==thick)['NRT_proxy_dpa_per_pulse'][str(ed)]*s['pulses_per_year'] for r,s in zip(rr,stages)) for ed in eds}} for thick in [50.,100.]],
   'depth_bins':[{'depth_um':list(edges[i:i+2]),'NRT_proxy_dpa_year_Ed25':sum(r['depth_bins'][i]['NRT_proxy_dpa_per_pulse']['25.0']*s['pulses_per_year'] for r,s in zip(rr,stages))} for i in range(len(edges)-1)]})
 qmax=max(s['peak_kinetic_plus_advected_internal_W_m2'] for s in stages)
 thermal=[]
 for thick in [50.,100.]:
  for k in [1.,4.,10.,25.,50.]:
   length=thick*1e-6
   thermal.append({'thickness_um':thick,'skin_k_W_m_K':k,'material_K_plus_U_peak_W_m2':qmax,'full_skin_resistance_rise_K':qmax*length/k,'increment_above_bulk_k_skin_K':qmax*length*(1/k-1/kbulk),'skin_diffusion_time_s':rho*cp*length**2/k,'maximum_permitted_contact_resistance_for_50K_total_m2_K_W':max(50/qmax-length/k,0.),'minimum_k_for_50K_skin_only_W_m_K':qmax*length/50})
 radpath=base/'component-design/mixture-radiation.json';rad=json.loads(radpath.read_text());hashes[str(radpath)]=sha(radpath);prompt=[]
 prompt_reference_case='C13-H4-fixed_energy-z64-k0p001-t1';source=rad['traces'][prompt_reference_case];assert abs(rad['radius_m']/radius-1)<1e-12
 for proxy in ['outer_T','surface_flux_equivalent_Teff','outer_T_quarter']:
  rr=source['cases'][proxy];m=rr['materials']['carbon']
  for k in [1.,4.,10.,50.]:
   resolved=m['resolved_surface_no_conduction_sensible_rise_K']*m['cp_J_kg_K']/cp
   soft=m['below1keV_all_surface_conduction_peak_K']*math.sqrt(m['k_W_m_K']*m['cp_J_kg_K']/(k*cp))
   prompt.append({'spectrum_proxy':proxy,'fluence_weighted_kT_keV':rr['weighted_kT_keV'],'skin_k_W_m_K':k,'resolved_no_conduction_rise_K':resolved,'soft_tail_half_space_peak_K':soft,'sum_of_separate_maxima_K':resolved+soft,'first50um_energy_transfer_J_m2':m['first50um_energy_transfer_J_m2']})
 repair=[];factory_per_flow=plant['factory_mass_per_repaired_kg_s'];baseflow=plant['routine_repaired_kg_s']/area
 for thick in [50.,100.]:
  sigma=rho*thick*1e-6
  for name,period in [('annual',YEAR),('quarterly',YEAR/4),('monthly',YEAR/12),('weekly',7*86400.),('daily',86400.)]:
   flow=sigma/period
   repair.append({'thickness_um':thick,'schedule':name,'period_s':period,'recycled_kg_m2_year':flow*YEAR,'recycled_kg_s_whole_station':flow*area,'flow_fraction_of_baseline_routine_capacity':flow/baseflow,'process_W_m2_at_MJ_kg':{str(e):flow*e*1e6 for e in [2.,3.1,10.]},'factory_kg_m2_at_existing_capacity_coefficient':flow*factory_per_flow,'fresh_carbon_required_if_closed_reuse_kg_s':0.})
 out={'scope':'Selected64 source and cadence damage sensitivity, not measured component lifetime or qualification. Planar CSDA path identified with depth; straggling/roughness/recoil cascades omitted. Ed sensitivities and stopping multipliers are scenarios, not credible intervals. Total ZBL elastic transfer is used as a conservative-leaning NRT damage-energy proxy, without recoil electronic partition or dynamic annealing; electronic excitation damage is separate.',
 'code_sha256':sha(__file__),'input_sha256':hashes,'radius_m':radius,'area_m2':area,'selected_cases':plant['selected_cases'],'graphite':{'density_kg_m3':rho,'cp_J_kg_K':cp,'bulk_k_W_m_K':kbulk,'carbon_atomic_mass_u':12.011},'NRT_formula':'0.8 * elastic_transfer / (2 * Ed * number_target_atoms)',
 'stages':stages,'annual_CNO_fluence_m2':sum(s['CNO_fluence_m2']*s['pulses_per_year'] for s in stages),'annual_all_nuclei_fluence_m2':sum(s['all_nuclei_fluence_m2']*s['pulses_per_year'] for s in stages),'annual':annual,'material_thermal':thermal,'prompt_sensitivities':prompt,'prompt_reference_case':prompt_reference_case,'prompt_reference_matches_selected':prompt_reference_case==plant['selected_cases']['C13'],'selected_C13_prompt_fluence_relative_to_reference':next(s for s in stages if s['stage']=='C13')['prompt_fluence_J_m2']/source['prompt_fluence_J_m2'],
 'prompt_scope':'Imported homogeneous first-interaction depth proxy. Resolved-depth no-conduction maximum scales inversely with cp and is independent of conductivity. Only the sub-keV surface tail is rescaled with effusivity. Treating the low-k skin as a half-space is conservative relative to an ideally contacted higher-k backing for surface heat. A true layered photon/electron transient and interfaces remain uncomputed. Material and prompt peaks occur at different times.',
 'restoration':repair,'restoration_scope':'Gross internally recycled skin throughput. Full reworking is a sensitivity, not a prescribed action at each period. Existing factory capacity is held fixed; capacity above it must be reassigned/added. At 2 MJ/kg this is a chosen complete process cost, with no extra atomization energy. Retained matter is not fresh feed. No achieved annual skin life is inferred.',
 'baseline_routine_capacity_kg_m2_year':baseflow*YEAR,'baseline_factory_kg_m2':plant['mass_breakdown_kg']['replacement_factory_capacity']/area,'baseline_net_W_m2':plant['net_W']/area,
 'validation':{'source_baryon_max_relative':max(abs(s['source_checks']['baryon_mass_relative']) for s in stages),'source_smooth_energy_max_relative':max(abs(s['source_checks']['smooth_kinetic_relative']) for s in stages),'all_damage_nonnegative':all(r['elastic_transfer_J_m2']>=0 for s in stages for r in s['rows']),'all_damage_less_than_kinetic':all(r['elastic_fraction_of_kinetic']<=1 for s in stages for r in s['rows'])}}
 a.output.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
 print(json.dumps({'annual_CNO_m2':out['annual_CNO_fluence_m2'],'annual_nominal':annual[1],'peak_W_m2':qmax,'validation':out['validation']},indent=2))
if __name__=='__main__':main()
