"""Actual outgoing blackout duty, magnetic volume and H/H2 extension requirements."""
from pathlib import Path
import json,csv,math,hashlib
from scipy.optimize import brentq
P=Path(__file__).resolve().parent;R=P.parents[2];evpath=R/'component-design/final-family-dispatch/sigma350/events.csv';schedpath=evpath.with_name('schedule.json')
s=json.loads((P/'screen.json').read_text());b=s['baseline'];sched=json.loads(schedpath.read_text());ev=list(csv.DictReader(evpath.open()));horizon=sched['count_horizon_s'];e=1.602176634e-19;me=9.1093837139e-31;eps=8.8541878188e-12;mu0=1.25663706212e-6;kb=1.380649e-23;mu=1.66053906660e-27;sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
duties=[]
for pad in [0.,1.,10.,100.]:
 spans=sorted((max(0,float(x['focus_s'])-1),min(horizon,float(x['main_capture_end_s'])+pad)) for x in ev if float(x['focus_s'])<horizon);merged=[]
 for lo,hi in spans:
  if merged and lo<=merged[-1][1]:merged[-1][1]=max(merged[-1][1],hi)
  else:merged.append([lo,hi])
 off=sum(hi-lo for lo,hi in merged);duties.append({'pre_focus_guard_s':1.,'post_main_capture_recovery_s':pad,'blackout_s':off,'live_fraction':1-off/horizon})
D=next(x['live_fraction'] for x in duties if x['post_main_capture_recovery_s']==10.)
ve=math.sqrt(2*100*e/me);vhe=b['He_mean_speed_m_s'];khe=3.5e-21*ve;a=1.;L=.01;r=.001;C=.5;eta=.5;g=.01

def conductance(nu,v):
 dk=2*r*v/3;c=C*v/4;si=math.sqrt(dk*nu)*math.tanh(L*math.sqrt(nu/dk));return si*c/(si+c)

def row(duty):
 nu=brentq(lambda n:4*a*eta*duty*conductance(n,vhe)/vhe-g,1e-12,1e8);ne=nu/khe;omega=math.sqrt(2*40*e/me)/r;B=max(.2,math.sqrt(4*me*ne/eps),me/e*(omega+ne*e*e/(2*eps*me*omega)));N=ne*a*L
 gaps=[]
 for gap in [.01,.05,.2,.3,.5]:
  UB=gap*B*B/(2*mu0);gaps.append({'effective_field_volume_per_shell_area_m':gap,'field_energy_J_m2':UB,'ideal_magnet_mass_kg_m2_at_BH100kJ_m3':2*UB*7500/1e5,'mass_fraction_of_existing_plant':2*UB*7500/1e5/b['plant_total_mass_per_shell_area_kg_m2']})
 return {'duty':duty,'He_electron_density_m3':ne,'B_T':B,'electron_column_m2':N,'electron_refresh_W_m2_at1keV_0p1s_no_recovery':duty*N*1000*e/.1,'self_potential_V':ne*e*r*r/(4*eps),'field_volume_sensitivities':gaps}
selected=row(D);rows=[row(d) for d in [1.,.8,D,.5,.25]]

def beb(T,B,U,N,Q):
 t=T/B;u=U/B;S=4*math.pi*(.52918e-10)**2*N*(13.6057/B)**2
 return 0. if t<=1 else S/(t+u+1)*(Q*.5*math.log(t)*(1-1/t**2)+(2-Q)*(1-1/t-math.log(t)/(t+1)))
const={'H':{'B_eV':13.6057,'U_eV':13.6057,'N':1.,'Q':.5668,'mass_u':1.00782503223},'H2':{'B_eV':15.43,'U_eV':15.98,'N':2.,'Q':1.,'mass_u':2*1.00782503223}}
chem=[]
for name,c in const.items():
 sig=beb(100,c['B_eV'],c['U_eV'],c['N'],c['Q']);v=math.sqrt(8*kb*1000/(math.pi*c['mass_u']*mu));nu=selected['He_electron_density_m3']*sig*ve;Seff=conductance(nu,v);per_eta=4*a*D*Seff/v
 chem.append({'species':name,'NIST_BEB_constants':c,'cross_section_at100eV_m2':sig,'mean_speed_1000K_m_s':v,'same_cloud_g_if_eta05':per_eta*.5,'required_net_transfer_eta_for_g001':.01/per_eta,'scope':'Electron ionization and neutral access only. H2 fragmentation,H/H2 chemistry,reflection,hot trapping and eventual transfer remain distinct; this is not a qualifiedH pump.'})
checks=[]
for rr in rows:checks.append({'duty':rr['duty'],'positive_B':rr['B_T']>0,'positive_mass':all(z['ideal_magnet_mass_kg_m2_at_BH100kJ_m3']>0 for z in rr['field_volume_sensitivities'])})
assert all(c['positive_B'] and c['positive_mass'] for c in checks)
result={'scope':__doc__,'code_sha256':sha(__file__),'sources_sha256':{str(evpath):sha(evpath),str(schedpath):sha(schedpath),str(P/'screen.json'):sha(P/'screen.json')},'horizon_s':horizon,'dispatch_events':len(ev),'blackout_scope':'Close1s beforefocus through actual99.9%main outgoing capture end, plus recovery time. Incoming launch is not globally blacked out. Omittedslowtails and localsourceapertures remain separate. Offintervals are clipped and unioned across exact370day schedule.','duties':duties,'selected_actual_duty':selected,'duty_sensitivities':rows,'H_extension':chem,'H_primary_sources':['https://physics.nist.gov/cgi-bin/Ionization/table.pl?ionization=H2','https://physics.nist.gov/cgi-bin/Ionization/ion_data.php?differential=Y&id=HI&initial=&ision=I','https://pml.nist.gov/PhysRefData/Ionization/intro.html'],'placement':'Candidate warm1cm shallowcell array on chamber-facing side of existing150um collectionskin, with retractable/shielded front grids during prompt/outgoing main windows. C=.5 is molecular admission, not opaque primary-interceptionfraction. The actual exposedsurface optical/particlegeometry and primaryheat/gasreturn maps are unclosed. Cooler magnets behindbacking/insulation use0.2–0.5m fieldvolume sensitivities; nearfront cooledislands require their own insulation and pulseprotection ledger.','magnet_qualification':'Fieldenergy/permanentmagnetloadline lowerdiagnostics at assumedBHmax100kJ/m3,density7500. No1000K permanentmagnet or superconductingmaterial assumed. Realfringing/yokes/coercivity/thermalpath add requirements.','double_counting':'g acts only on residualatoms alreadyin cavity after localrouting. Any interceptionof newlydesorbedgas by thesamecells changes upstreamresidualfraction andmustnot also be credited as removalofthatgasfromthecavity.','checks':checks}
(P/'placement.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'duties':duties,'selected':selected,'H':chem},indent=2))
