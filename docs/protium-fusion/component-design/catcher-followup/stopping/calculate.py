"""Current-family charged-particle stopping, intact-slab transmission and elastic proxy.
H/He use updated NIST evaluated stopping; CNO uses archived measured electronic
stopping plus ZBL. Ranges are distinct from retention/release probabilities.
"""
from pathlib import Path
import sys,json,math,hashlib
import numpy as np
from scipy.integrate import cumulative_trapezoid
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P));import spectra as sp
R=sp.R;EV=1.602176634e-19;YEAR=365.25*86400.;RHO=2000.;SIGMA=350.
sha=sp.sha;write=sp.write
EDGES=np.array([0,.01,.1,.25,.5,1,2,5,10,25,50,100,250,1000,175000.])
PROBS=[.001,.01,.1,.5,.9,.99,.999,.9999,1.]
def q(x,w):
 ix=np.argsort(x);c=np.cumsum(w[ix]);return {str(p):float(x[ix[min(np.searchsorted(c,p*c[-1]),len(ix)-1)]]) for p in PROBS}
def loginterp(x,xp,yp):return np.exp(np.interp(np.log(np.maximum(x,1e-300)),np.log(xp),np.log(yp)))
class Tables:
 def __init__(self,st):self.st=st;self.data=json.loads((P/'sources/nist-graphite-tables.json').read_text())['tables'];self.cache={};self.diagnostics={}
 def get(self,z,a,factor=1.):
  key=(z,a,factor)
  if key in self.cache:return self.cache[key]
  if z not in [1,2]:
   E,Rg,D,se,sn=self.st.table('carbon',z,a,factor);out=dict(E=E,R=Rg,D=D,se=se,sn=sn,projected=None,published_R=None,source='IAEA measured CNO electronic / archived analytic other ions + ZBL; straight-CSDA depth proxy')
  else:
   assert (z,a) in [(1,1),(2,4)]
   name='p' if z==1 else 'he4';d=np.array(self.data[name]['rows']);E=np.unique(np.r_[np.geomspace(1e-12,100.,24000),d[d[:,0]<=100,0]])
   se=loginterp(E,d[:,0],d[:,1])/1000;sn=loginterp(E,d[:,0],d[:,2])/1000;low=E<d[0,0];se[low]=d[0,1]/1000*np.sqrt(E[low]/d[0,0]);sn[low]=self.st.nuclear(E[low],z,a,'carbon')
   raw0=cumulative_trapezoid(1/(se+sn),E,initial=0);raw=cumulative_trapezoid(1/(factor*se+sn),E,initial=0);anchor=d[0,4]*1000;Rnom=anchor+raw0-np.interp(d[0,0],E,raw0);slope=np.log(d[1,4]/d[0,4])/np.log(d[1,0]/d[0,0]);Rnom[low]=anchor*(E[low]/d[0,0])**slope
   ratio=np.divide(raw,raw0,out=np.ones_like(raw),where=raw0>0);ratio[0]=ratio[1];Rg=Rnom*ratio;D=cumulative_trapezoid(sn/(factor*se+sn),E,initial=0);det=loginterp(E,d[:,0],d[:,5]);published=loginterp(E,d[:,0],d[:,4])*1000;published[low]=Rnom[low]
   assert np.all(np.diff(Rg)>0) and np.all(D>=0) and np.all(D<=E*(1+1e-12)),(z,a,factor,float(np.min(np.diff(Rg))),float(np.max(D/E)))
   out=dict(E=E,R=Rg,D=D,se=factor*se,sn=sn,projected=Rg*det,published_R=published,source='Updated NIST ICRU90 electronic+elastic stopping integrated above1keV; published1keV CSDA anchor. Evaluated detour factor gives mean-projected-depth proxy; sub1keV completion explicit.')
   if factor==1:
    self.diagnostics[name]=[{'energy_MeV':e,'integrated_CSDA_kg_m2':float(np.interp(e,E,Rg)/100),'published_CSDA_kg_m2':float(np.interp(e,E,published)/100),'relative_difference':float(np.interp(e,E,Rg)/np.interp(e,E,published)-1)} for e in [.001,.003,.01,.03,.1,.3,1.,3.,10.]]
  self.cache[key]=out;return out

def report_group(ranges,counts,mass,E,projected=None):
 total=float(counts.sum());out={'nuclei_m2':total,'zero_CSDA_subgrid_number_fraction':float(counts[ranges==0].sum()/total),'mass_kg_m2_baryon':float(mass.sum()),'kinetic_J_m2':float((counts*E).sum()*EV*1e6),'CSDA_column_kg_m2_quantiles_number':q(ranges,counts),'CSDA_depth_um_quantiles_number':q(ranges/RHO*1e6,counts),'CSDA_depth_um_quantiles_mass':q(ranges/RHO*1e6,mass),'CSDA_depth_um_mean_number':float(counts@ranges/total/RHO*1e6),'stopped_primary_fraction_CSDA_at_sigma350':float(counts[ranges<=SIGMA].sum()/total),'unscattered_transmission_number_fraction_CSDA':float(counts[ranges>SIGMA].sum()/total),'straight_depth_absorbing_boundaries_rear_probability':{str(L):float(counts@np.minimum(ranges/RHO/(L*1e-6),1)/total) for L in [25.,50.,100.]} ,'CSDA_stopping_depth_bins':[]}
 for lo,hi in zip(EDGES[:-1],EDGES[1:]):
  mask=(ranges/RHO*1e6>=lo)&(ranges/RHO*1e6<hi);out['CSDA_stopping_depth_bins'].append({'depth_um':[float(lo),float(hi)],'nuclei_m2':float(counts[mask].sum()),'mass_kg_m2_baryon':float(mass[mask].sum())})
 if projected is not None:
  out['mean_projected_depth_um_quantiles_number']=q(projected/RHO*1e6,counts);out['mean_projected_depth_um_mean_number']=float(counts@projected/total/RHO*1e6);out['mean_projected_depth_absorbing_boundaries_rear_probability']={str(L):float(counts@np.minimum(projected/RHO/(L*1e-6),1)/total) for L in [25.,50.,100.]}
 return out

def main():
 rs,pm,st,decay,weak,hashes=sp.setup();tables=Tables(st);mapping=json.loads((R/'component-design/final-family/delivered-B1-C13H2-stage-map.json').read_text());interface=json.loads((P/'surface-spectra.json').read_text());checks=[];rows=[];surface={'source_surface_spectra_sha256':sha(P/'surface-spectra.json'),'density_kg_m3':RHO,'scope':'Packet arrays match the32-point surface spectrum exactly. CSDA is path length; projected is evaluated mean-depth proxy, not a straggling/retention distribution.','stages':{}}
 for cat,decl in mapping['stage_map'].items():
  record=R/decl['record'];x=sp.packets(record,Path(decl['trace']),pm,decay,order=64);d=x['record'];stage=d['stage'];n=x['n'];eu=x['eu'];A=pm.AA;Z=pm.ZZ;total_ke=float(np.sum(n*A*eu[:,None])*pm.MEV);period=interface['stages'][stage]['period_s'];charged=Z>0;stage_rows=[]
  for factor in [1.,.5,2.]:
   ranges=np.zeros_like(n);project=np.zeros_like(n);D=np.zeros_like(n);cum=np.zeros(len(EDGES));species={};max_bound=0.;transE=0.;lowDupper=0.
   for j in np.where(charged)[0]:
    count=n[:,j];N=float(count.sum())
    if N==0:continue
    tab=tables.get(int(Z[j]),int(A[j]),factor);E=A[j]*eu;assert E.max()<tab['E'][-1];rr=np.interp(E,tab['E'],tab['R']);dd=np.interp(E,tab['E'],tab['D']);ranges[:,j]=rr/100;D[:,j]=dd
    bound=float(np.interp(A[j]*.5*pm.MU*x['vmax']**2/pm.MEV,tab['E'],tab['R'])/100);max_bound=max(max_bound,bound)
    rem=np.maximum(rr-SIGMA*100,0);exitE=np.where(rem>0,np.interp(rem,tab['R'],tab['E']),0);transE+=float(count@exitE*pm.MEV)
    proj=None
    if tab['projected'] is not None:proj=np.interp(E,tab['E'],tab['projected'])/100;project[:,j]=proj;lowDupper+=N*.001*pm.MEV
    group=report_group(ranges[:,j],count,count*A[j]*pm.MU,E,proj);group.update(source=tab['source'],maximum_endpoint_velocity_CSDA_column_kg_m2=bound,elastic_transfer_J_m2=float(count@dd*pm.MEV),below1keV_incident_number_fraction=float(count[E<.001].sum()/N),NRT_proxy_events_m2={str(ed):float(.4*(count@dd)*1e6/ed) for ed in [25.,40.]})
    species[pm.SPECIES[j]]=group
    for k,depth in enumerate(EDGES):
     residual=np.maximum(rr-RHO*depth*1e-6*100,0);outE=np.where(residual>0,np.interp(residual,tab['R'],tab['E']),0);dout=np.where(outE>0,np.interp(outE,tab['E'],tab['D']),0);cum[k]+=float(count@(dd-dout)*pm.MEV)
   elastic=float(np.sum(n*D)*pm.MEV);groups={}
   for name,mask in [('H',Z==1),('He',Z==2),('CNO',(Z>=6)&(Z<=8)),('all_charged',charged)]:
    c=n[:,mask].ravel();rg=ranges[:,mask].ravel();mass=(n[:,mask]*A[mask]*pm.MU).ravel();E=np.broadcast_to(eu[:,None]*A[mask],n[:,mask].shape).ravel();proj=project[:,mask].ravel() if name in ['H','He'] else None;groups[name]=report_group(rg,c,mass,E,proj)
   bins=[]
   for lo,hi,en in zip(EDGES[:-1],EDGES[1:],np.diff(cum)):
    atoms=RHO*(hi-lo)*1e-6/(12.011*pm.MU);bins.append({'depth_um':[float(lo),float(hi)],'elastic_transfer_J_m2':float(en),'NRT_proxy_events_m2':{str(ed):float(.4*en/EV/ed) for ed in [25.,40.]},'NRT_proxy_dpa':{str(ed):float(.4*en/EV/ed/atoms) for ed in [25.,40.]}})
   skins=[]
   for depth in [25.,50.,100.]:
    en=float(cum[list(EDGES).index(depth)]);skins.append({'depth_um':depth,'elastic_transfer_J_m2':en,'NRT_proxy_events_m2':{str(ed):.4*en/EV/ed for ed in [25.,40.]},'NRT_proxy_dpa':{str(ed):.4*en/EV/ed/(RHO*depth*1e-6/(12.011*pm.MU)) for ed in [25.,40.]}})
   stage_rows.append({'electronic_stopping_multiplier':factor,'groups':groups,'species':species,'all_charged_elastic_transfer_J_m2':elastic,'all_charged_NRT_proxy_events_m2':{str(ed):.4*elastic/EV/ed for ed in [25.,40.]},'all_charged_NRT_proxy_dpa_whole_slab':{str(ed):.4*elastic/EV/ed/(SIGMA/(12.011*pm.MU)) for ed in [25.,40.]},'elastic_transfer_depth_bins':bins,'skins':skins,'maximum_endpoint_velocity_CSDA_column_kg_m2':max_bound,'slab_to_maximum_CSDA_column_ratio':SIGMA/max_bound,'intact_slab_primary_transmitted_kinetic_J_m2_CSDA':transE,'HHe_unresolved_sub1keV_elastic_energy_interval_width_J_m2_upper':lowDupper})
   checks.append({'stage':stage,'factor':factor,'nonnegative_ranges':bool(np.all(ranges[:,charged][n[:,charged]>0]>=0)),'elastic_less_than_incident':elastic<=total_ke,'depth_energy_closure':bool(abs(cum[-1]/elastic-1)<1e-10),'nonnegative_depth_bins':bool(np.min(np.diff(cum))>-1e-7),'zero_primary_CSDA_transmission':transE==0})
  periodrow={'stage':stage,'record_path':str(record),'record_sha256':sha(record),'trace_sha256':sha(Path(decl['trace'])),'period_s':period,'pulses_per_year':YEAR/period,'incident_kinetic_J_m2':total_ke,'rows':stage_rows};rows.append(periodrow)
  s32=interface['stages'][stage];surf={'period_s':period,'species':{}}
  for s,z,a in [('p',1,1),('he4',2,4)]:
   tab=tables.get(z,a,1.);E=np.array(s32['packet_energy_MeV_u'])*a;rr=np.interp(E,tab['E'],tab['R'])/100;pr=np.interp(E,tab['E'],tab['projected'])/100;pub=np.interp(E,tab['E'],tab['published_R'])/100;surf['species'][s]={'packet_CSDA_depth_um':(rr/RHO*1e6).tolist(),'packet_mean_projected_depth_um':(pr/RHO*1e6).tolist(),'packet_published_NIST_CSDA_depth_um':(pub/RHO*1e6).tolist(),'packet_nuclei_m2':s32['species'][s]['packet_nuclei_m2'],'range_model':tab['source']}
  surface['stages'][stage]=surf
  print(json.dumps({'stage':stage,'groups':{k:{q:v[q] for q in ['mass_kg_m2_baryon','CSDA_depth_um_quantiles_number','CSDA_depth_um_mean_number']} for k,v in stage_rows[0]['groups'].items() if k in ['H','He','CNO']},'elastic_J_m2':stage_rows[0]['all_charged_elastic_transfer_J_m2'],'max_column':stage_rows[0]['maximum_endpoint_velocity_CSDA_column_kg_m2']}),flush=True)
 annual=[]
 for f in [1.,.5,2.]:
  selected=[next(r for r in s['rows'] if r['electronic_stopping_multiplier']==f) for s in rows];annual.append({'electronic_stopping_multiplier':f,'elastic_transfer_J_m2_year':sum(s['pulses_per_year']*r['all_charged_elastic_transfer_J_m2'] for s,r in zip(rows,selected)),'NRT_proxy_events_m2_year':{str(ed):sum(s['pulses_per_year']*r['all_charged_NRT_proxy_events_m2'][str(ed)] for s,r in zip(rows,selected)) for ed in [25.,40.]}})
 for check in checks:assert all(v for k,v in check.items() if k not in ['stage','factor']),check
 out={'scope':__doc__,'radius_m':1.1e9,'graphite_density_kg_m3':RHO,'slab_sigma_kg_m2':SIGMA,'slab_thickness_m':SIGMA/RHO,'source_surface_spectra_sha256':sha(P/'surface-spectra.json'),'primary_HHe_source':'https://www.nist.gov/system/files/documents/2017/04/26/newstar.pdf','nist_tables_sha256':sha(P/'sources/nist-graphite-tables.json'),'source_code_sha256':hashes,'code_sha256':sha(__file__),'rows':rows,'annual':annual,'checks':checks,'NIST_published_vs_integrated_range_diagnostic':tables.diagnostics,'range_scope':'Intact solid normal-incidence CSDA path model. H/He projected depths are mean projected-range proxies using evaluated detour factors, not stochastic implantation profiles. No retention, backscatter, secondary yield, pore/gap transmission, diffusion or release law inferred. Electronic multipliers0.5/2 are sensitivity choices, not probabilistic intervals.','NRT_scope':'0.8*total elastic transfer/(2*Ed); Ed25/40eV. Total nuclear elastic transfer is a damage-energy proxy without electronic partition of recoils, cascade efficiency or thermal annealing. Counts are not surviving defects or mandatory atoms to replace. Electronic-excitation damage and neutral neutron transport are separate.','low_energy_scope':'Below1keV: electronic sqrt(E) matched to NIST; ZBL elastic stopping; CSDA completion matched to published1keV range using lowest published range slope. Below1keV elastic energy has at most1keV/incidentHHe nucleus uncertainty; its explicit aggregate interval width is reported. Heavy-ion packets below the archived1e-10MeV stopping grid have unresolved zero path in this model; their number fraction is reported explicitly. Incoming charge-state equilibration, material temperature/microstructure and range straggling remain separate.'}
 write(P/'stopping-results.json',out);write(P/'surface-stopping.json',surface)
if __name__=='__main__':main()
