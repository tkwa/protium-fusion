"""Bounded thin-cover and incoming-gas stopping comparators; no integrity law."""
from pathlib import Path
import sys,json,math
import numpy as np
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P));import calculate as ca;import spectra as sp

def main():
 rs,pm,st,decay,weak,hashes=sp.setup();tab=ca.Tables(st);mapping=json.loads((sp.R/'component-design/final-family/delivered-B1-C13H2-stage-map.json').read_text());cache={};rows=[]
 for cat,decl in mapping['stage_map'].items():
  x=sp.packets(sp.R/decl['record'],Path(decl['trace']),pm,decay,order=64);n=x['n'];eu=x['eu'];stage=x['record']['stage']
  for material in ['carbon','aluminium']:
   rho=st.MAT[material]['rho']
   for thick_nm in [100.,1000.]:
    sigma=rho*thick_nm*1e-9;out={};loss=elastic=incident=0.;stopped_mass=0.;total_charged_mass=0.
    for j in np.where(pm.ZZ>0)[0]:
     N=n[:,j];z=int(pm.ZZ[j]);a=int(pm.AA[j]);E=a*eu
     if N.sum()==0:continue
     if material=='carbon':t=tab.get(z,a,1.)
     else:
      key=(material,z,a)
      if key not in cache:
       eg,rg,D,se,sn=st.table(material,z,a);cache[key]={'E':eg,'R':rg,'D':D}
      t=cache[key]
     rr=np.interp(E,t['E'],t['R']);rem=np.maximum(rr-sigma*100,0);exitE=np.where(rem>0,np.interp(rem,t['R'],t['E']),0);lossE=E-exitE;dn=np.interp(E,t['E'],t['D'])-np.where(exitE>0,np.interp(exitE,t['E'],t['D']),0)
     inc=float(N@E*pm.MEV);ls=float(N@lossE*pm.MEV);el=float(N@dn*pm.MEV);stop=float(N[rem==0].sum());mass=float(N.sum()*a*pm.MU);incident+=inc;loss+=ls;elastic+=el;stopped_mass+=stop*a*pm.MU;total_charged_mass+=mass
     if pm.SPECIES[j] in ['p','he4'] or 6<=z<=8:
      out[pm.SPECIES[j]]={'incident_nuclei_m2':float(N.sum()),'stopped_primary_number_fraction_CSDA':stop/float(N.sum()),'primary_transmitted_number_fraction_CSDA':float(N[rem>0].sum()/N.sum()),'incident_kinetic_J_m2':inc,'primary_energy_loss_J_m2':ls,'elastic_transfer_J_m2':el}
    atom_mass,_,atomsper=st.target(material);target_atoms=sigma/(atom_mass*pm.MU)*atomsper
    rows.append({'stage':stage,'material':material,'assumed_density_kg_m3':rho,'thickness_nm':thick_nm,'sigma_kg_m2':sigma,'incident_kinetic_J_m2':incident,'primary_energy_loss_J_m2':loss,'primary_energy_transmitted_fraction':1-loss/incident,'stopped_charged_mass_fraction_CSDA':stopped_mass/total_charged_mass,'elastic_transfer_J_m2':elastic,'NRT_proxy_events_m2':{str(ed):.4*elastic/ca.EV/ed for ed in [25.,40.]},'NRT_proxy_dpa':{str(ed):.4*elastic/ca.EV/ed/target_atoms for ed in [25.,40.]},'incident_charged_particles_per_target_atom':float(n[:,pm.ZZ>0].sum()/target_atoms),'species':out})
 sys.path.insert(0,str(sp.R/'component-design'));import material_cycle as mc
 gas=[]
 for v in [.3e6,2.4e6,4.6e6]:
  for name,z,a in [('p',1,1),('he4',2,4)]:
   m=mc.atomic_kg(name);E=.5*m*v*v/pm.MEV;t=tab.get(z,a,1.);rg=float(np.interp(E,t['E'],t['R'])/100);rp=float(np.interp(E,t['E'],t['projected'])/100);gas.append({'species':name,'relative_speed_m_s':v,'atomic_mass_kg':m,'kinetic_MeV':E,'graphite_CSDA_depth_um_at2000kg_m3':rg/2000*1e6,'graphite_mean_projected_depth_um_at2000kg_m3':rp/2000*1e6,'graphite_CSDA_column_kg_m2':rg})
 for row in rows:assert 0<=row['primary_energy_transmitted_fraction']<=1 and row['elastic_transfer_J_m2']<=row['primary_energy_loss_J_m2']*(1+1e-8)
 result={'scope':'Thin intact normal-incidence cover comparator. Primary projectile transmission and energy removal only; no gas impermeability, structural integrity, track-to-hole conversion, secondary electron retention, thermal response or stopping-to-retention inference. Carbon uses current updated model; aluminium uses archived analytic H/He and measured CNO stopping plus ZBL. Aluminium density2300kg/m3 is the earlier collection-layer assumption; support/phase compatibility remains separate.','source_spectra_sha256':sp.sha(P/'surface-spectra.json'),'code_sha256':sp.sha(__file__),'rows':rows,'incoming_gas_comparator':gas,'incoming_gas_scope':'Atomic-mass kinetic energy in pellet frame, NIST graphite stopping. These sub-MeV gas impacts differ from outgoing source ions; full energy retention remains an upper deposition assumption.','NRT_scope':'Total elastic recoil transfer proxy atEd25/40, not surviving defects, holes or failed gas sealing. Ionization-track overlap alone does not imply perforation.'}
 sp.write(P/'cover-and-gas.json',result);print(json.dumps({'incoming_gas':gas,'carbon_light_cover':[{'stage':r['stage'],'nm':r['thickness_nm'],'H_stopped':r['species']['p']['stopped_primary_number_fraction_CSDA'],'He_stopped':r['species']['he4']['stopped_primary_number_fraction_CSDA'],'elastic_J_m2':r['elastic_transfer_J_m2'],'dpa25':r['NRT_proxy_dpa']['25.0']} for r in rows if r['material']=='carbon']},indent=2))
if __name__=='__main__':main()
