"""Independent saved-record, units, integral, and provenance checks; no physics rerun."""
from pathlib import Path
import json,math,hashlib
P=Path(__file__).resolve().parent;R=P.parents[2];ROOT=R.parents[1]
MU=1.66053906660e-27;EV=1.602176634e-19;YEAR=365.25*86400
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
load=lambda p:json.loads(Path(p).read_text())
checks=[]
def ck(label,ok):
 checks.append({'check':label,'passed':bool(ok)})
 assert ok,label
def close(label,a,b,tol=1e-10):ck(label,abs(a-b)<=tol*max(abs(a),abs(b),1e-290))
s=load(P/'surface-spectra.json');d=load(P/'stopping-results.json');i=load(P/'surface-stopping.json');a=load(P/'cover-and-gas.json')
ck('source map hash',sha(s['source_map'])==s['source_map_sha256']);ck('cadence hash',sha(s['cadence_source'])==s['cadence_source_sha256'])
for path,h in s['helper_sha256'].items():ck('helper '+path,sha(ROOT/path)==h)
for obj,key in [(d,'source_surface_spectra_sha256'),(i,'source_surface_spectra_sha256'),(a,'source_spectra_sha256')]:ck('surface interface hash '+key,obj[key]==sha(P/'surface-spectra.json'))
for name,obj in [('spectra.py',s),('calculate.py',d),('auxiliary.py',a)]:ck('code '+name,sha(P/name)==obj['code_sha256'])
ns=load(P/'sources/nist-source.json');ck('NIST PDF hash',sha(P/'sources/nist-newstar-2017.pdf')==ns['sha256'])
ck('NIST parsed table hash',sha(P/'sources/nist-graphite-tables.json')==d['nist_tables_sha256'])
close('sphere area',s['area_m2'],4*math.pi*(1.1e9)**2);close('slab thickness',d['slab_thickness_m'],350/2000)
for stage,row in s['stages'].items():
 label=stage+' '
 ck(label+'source record hash',sha(row['record_path'])==row['record_sha256']);ck(label+'raw trace hash',sha(row['trace_path'])==row['trace_sha256'])
 rec=load(row['record_path']);close(label+'total mass/area',row['total_mass_fluence_kg_m2_baryon'],rec['mass_kg']/s['area_m2']);close(label+'pulse cadence',row['pulses_per_year'],YEAR/row['period_s'])
 ck(label+'closure diagnostics',max(abs(v) for v in row['checks'].values())<1e-6)
 ck(label+'user rounded H+He check',abs(row['H_plus_He_ug_m2_baryon']-row['user_rounded_H_plus_He_ug_m2'])<(.005 if stage=='B1' else .05))
 close(label+'all-species baryon closure',sum(row['whole_pulse_species_wall_nuclei_m2'][name]*mass*MU for name,mass in zip(s['nuclear_species'],s['mass_numbers'])),row['total_mass_fluence_kg_m2_baryon'])
 close(label+'packet baryon mass',sum(row['packet_mass_kg_m2_baryon']),row['total_mass_fluence_kg_m2_baryon'])
 close(label+'smooth bulk kinetic energy',sum(m*E*EV*1e6/MU for m,E in zip(row['packet_mass_kg_m2_baryon'],row['packet_energy_MeV_u'])),row['smooth_kinetic_fluence_J_m2'])
 for name,A in [('p',1),('he4',4)]:
  x=row['species'][name];di=i['stages'][stage]['species'][name];c=x['packet_nuclei_m2'];prefix=label+name+' '
  close(prefix+'count integral',sum(c),x['nuclei_m2']);close(prefix+'mass conversion',x['nuclei_m2']*A*MU*1e9,x['mass_ug_m2_baryon']);close(prefix+'kinetic integral',sum(n*E*A*EV*1e6 for n,E in zip(c,row['packet_energy_MeV_u'])),x['incident_kinetic_J_m2'])
  ck(prefix+'identical packet ordering',di['packet_nuclei_m2']==c and len(c)==len(di['packet_CSDA_depth_um'])==len(row['packet_arrival_after_ignition_s']))
  ck(prefix+'mean projection bounded by path',all(0<=pr<=cs*(1+1e-12) for pr,cs in zip(di['packet_mean_projected_depth_um'],di['packet_CSDA_depth_um'])))
 for rr in next(x for x in d['rows'] if x['stage']==stage)['rows']:
  prefix=label+str(rr['electronic_stopping_multiplier'])+' '
  elastic=rr['all_charged_elastic_transfer_J_m2'];close(prefix+'depth-bin elastic closure',sum(b['elastic_transfer_J_m2'] for b in rr['elastic_transfer_depth_bins']),elastic)
  ck(prefix+'intact slab CSDA stop',rr['intact_slab_primary_transmitted_kinetic_J_m2_CSDA']==0 and rr['slab_to_maximum_CSDA_column_ratio']>900)
  for ed in [25.,40.]:
   close(prefix+'NRT identity '+str(ed),rr['all_charged_NRT_proxy_events_m2'][str(ed)],.4*elastic/(ed*EV))
   for skin in rr['skins']:
    atoms=2000*skin['depth_um']*1e-6/(12.011*MU)
    close(prefix+'skin dpa '+str(ed)+' '+str(skin['depth_um']),skin['NRT_proxy_dpa'][str(ed)]*atoms,.4*skin['elastic_transfer_J_m2']/(ed*EV))
for ar in d['annual']:
 selected=[next(r for r in st['rows'] if r['electronic_stopping_multiplier']==ar['electronic_stopping_multiplier']) for st in d['rows']]
 close('annual elastic '+str(ar['electronic_stopping_multiplier']),ar['elastic_transfer_J_m2_year'],sum(st['pulses_per_year']*r['all_charged_elastic_transfer_J_m2'] for st,r in zip(d['rows'],selected)))
 for ed in [25.,40.]:close('annual NRT '+str(ar['electronic_stopping_multiplier'])+' '+str(ed),ar['NRT_proxy_events_m2_year'][str(ed)],.4*ar['elastic_transfer_J_m2_year']/(ed*EV))
for row in a['rows']:
 label=row['stage']+' '+row['material']+' '+str(row['thickness_nm'])+' '
 close(label+'cover energy ledger',row['primary_energy_loss_J_m2']+row['primary_energy_transmitted_fraction']*row['incident_kinetic_J_m2'],row['incident_kinetic_J_m2'])
 ck(label+'elastic below removed energy',0<=row['elastic_transfer_J_m2']<=row['primary_energy_loss_J_m2'])
 for name in ['p','he4']:
  x=row['species'][name];close(label+name+' number balance',x['stopped_primary_number_fraction_CSDA']+x['primary_transmitted_number_fraction_CSDA'],1)
for row in a['incoming_gas_comparator']:
 label='gas '+row['species']+' '+str(row['relative_speed_m_s'])+' '
 close(label+'energy',row['kinetic_MeV']*EV*1e6,.5*row['atomic_mass_kg']*row['relative_speed_m_s']**2)
 close(label+'depth column',row['graphite_CSDA_depth_um_at2000kg_m3']*1e-6*2000,row['graphite_CSDA_column_kg_m2'])
result={'scope':'Saved-record algebra, units and hashes; no claim of experimental range, retention, damage or hardware validation.','passed':sum(x['passed'] for x in checks),'total':len(checks),'code_sha256':sha(__file__),'checks':checks}
(P/'validation.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'passed':result['passed'],'total':result['total']}))
