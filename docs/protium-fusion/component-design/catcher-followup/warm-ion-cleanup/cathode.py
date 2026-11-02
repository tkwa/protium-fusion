"""Controlled-energy residual-gas implantation into the existing routed graphite concept.
Monoenergetic mean-range and conditional mobile first-passage calculation only.
"""
from pathlib import Path
import sys,json,math,hashlib,importlib.util
import numpy as np
P=Path(__file__).resolve().parent;S=P.parent/'stopping'
sys.path.insert(0,str(S));import calculate as st;import spectra as sp
spec=importlib.util.spec_from_file_location('channel_routing',P.parent/'channel-routing/calculate.py');cr=importlib.util.module_from_spec(spec);spec.loader.exec_module(cr)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
rs,pm,old,decay,weak,hashes=sp.setup();tables=st.Tables(old)
rows=[];checks=[]
for eta in [1.,100.,1000.,10000.,100000.]:
 z,u,t=cr.solve(20.,.1,150.,eta,320,128)
 for V in [3e4,1e5,3e5,1e6]:
  for species,Z,A,nfragment in [('H+',1,1,1),('H2+',1,1,2),('He+',2,4,1)]:
   E=V/1e6/nfragment;tab=tables.get(Z,A);path=float(np.interp(E,tab['E'],tab['R'])/100/2000*1e6);depth=float(np.interp(E,tab['E'],tab['projected'])/100/2000*1e6)
   coll=float(np.interp(depth,np.r_[0,z,150.],np.r_[0,u,u[-1]]));tt=float(np.interp(depth,np.r_[0,z,150.],np.r_[0,t,t[-1]]))
   rows.append({'incident_ion':species,'acceleration_V':V,'fragment_energy_MeV':E,'fragments_per_incident_ion':nfragment,'graphite_CSDA_um':path,'mean_projected_depth_um':depth,'Dparallel_over_Dnormal':eta,'mobile_collection_fraction_at_mean_projected_depth':coll,'mobile_front_fraction_at_mean_projected_depth':1-coll,'exit_time_s_if_Dnormal_1um2_s':tt,'implant_below100nm_cap_in_mean_range_model':depth>.1})
   checks.append(-1e-9<=coll<=1+1e-9 and depth<=path and path<150)
 # Directly invert the numerical routing curve, rather than an exponential guess.
 good=z[u>=.99]
 if len(good):
  xmin=float(np.interp(.99,u,z));
  for name,Z,A,nfragment in [('H+',1,1,1),('H2+',1,1,2),('He+',2,4,1)]:
   tab=tables.get(Z,A);depth=tab['projected']/100/2000*1e6;emin=float(np.interp(xmin,depth,tab['E']))
   rows.append({'required_for_mobile_collection_fraction':.99,'incident_ion':name,'Dparallel_over_Dnormal':eta,'required_mean_projected_depth_um':xmin,'minimum_acceleration_V_in_mean_range_model':emin*1e6*nfragment})
e=1.602176634e-19;mu=1.66053906660e-27;eps=8.8541878128e-12;area=4*math.pi*(1.1e9)**2
# Exact routed-cap chemical-form handoff; preserve its stated collection scenario.
gaspath=P.parent/'gas/chemistry.json';gas=json.loads(gaspath.read_text())['baseline']
loads=[]
for name,Q,m in [('He+',gas['He_cavity_atomic_kg_s'],4.00260325413*mu),('H2+',gas['H_cavity_atomic_kg_s'],2*1.00782503223*mu),('H+',gas['H_cavity_atomic_kg_s'],1.00782503223*mu)]:
 J=Q/m/area*e
 loads.append({'species':name,'residual_cavity_mass_kg_s':Q,'one_pass_current_A_m2':J,'300kV_one_pass_electric_W_m2':J*3e5,'two_attempts_300kV_electric_W_m2':J*6e5,'scope':'H+ and H2+ are alternative chemical forms of the same H mass. Work excludes escaped secondary electrons and auxiliary electron-cloud refresh.'})
gap=[]
for L in [.1,.3,1.]:
 for name,m in [('H+',1.00782503223*mu),('H2+',2*1.00782503223*mu),('He+',4.00260325413*mu)]:
  V=3e5;v=math.sqrt(2*e*V/m);gap.append({'species':name,'gap_m':L,'mean_field_V_m':V/L,'cold_start_transit_s':2*L/v,'planar_unneutralized_Child_Langmuir_A_m2':4/9*eps*math.sqrt(2*e/m)*V**1.5/L**2,'exit_Larmor_radius_at0p203823T_m':m*v/e/.2038231872994329})
out={'scope':__doc__,'code_sha256':sha(__file__),'input_sha256':{str(p):sha(p)for p in[S/'calculate.py',S/'sources/nist-graphite-tables.json',P.parent/'channel-routing/calculate.py']},'geometry':{'pitch_um':20.,'cap_um':.1,'active_depth_um':150.,'boundary':'Ideal absorbing side channel, absorbing cavity front, reflecting rear, uniform lateral source; graphite diffusion tensor spatially constant.'},'rows':rows,'residual_load_scenarios':loads,'300kV_gap_scenarios':gap,'checks':{'number':len(checks),'passed':sum(checks)},'limits':['NIST mean projected range is used at one energy; it is not a range-straggling or reflection distribution. Consequently a calculated99% mobile routing is not an experimentally supported99% incident-ion collection efficiency.','H2+ is assumed to dissociate into equal-speed H fragments in the stopping material. Molecular stopping, dissociative ionization before acceleration and mixed charge states require their actual spectra.','The controlled cathode is at negative potential behind grounded shielding. Positive ions enter through extraction optics; ionizer electrons remain in their separate return circuit. A negatively biased secondary-suppression electrode around the cathode must return emitted electrons. Escaped secondary electrons add300keV each to the electrical debit.','The candidate shares the hot irradiated barrier, channels, mobile/trapped populations and regeneration qualifications of the main catcher. Ion optical transfer is still an independent factor in eta_net.','The electron-cloud B field does not by itself guide the accelerated heavy ions: the reported ion gyroradii require electrostatic optics and actual trajectories. High voltage vacuum insulation, desorbed gas, joints and pulse shielding remain hardware requirements.','The Child-Langmuir values are plane-diode scales, not a designed extraction aperture. Actual currents are many orders smaller even for sparse cathodes.','A shielded cathode can occupy a small fraction of wall area because gas throughput is low; this does not reduce the broad molecular admission area or the electron-cloud magnetic-field requirement.']}
out['input_sha256'][str(gaspath)]=sha(gaspath)
assert all(checks);(P/'cathode.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
print(json.dumps({'selected':[r for r in rows if r.get('acceleration_V')==3e5 and r['Dparallel_over_Dnormal'] in[100.,1000.,10000.]],'loads':loads,'checks':out['checks']},indent=2))
