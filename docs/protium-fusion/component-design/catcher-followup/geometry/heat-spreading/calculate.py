"""Fixed-mass heat-spreading sensitivity on actual final350 dispatch.
Effectivebulkconductivity changes the nonuniformmode, not the nuclear/source
energy,meanpickup,storageheatcapacity,skinorpostmass. Highk mayrepresent a
future heatspreader; no qualifiedhardware or materialproperty is asserted.
"""
from pathlib import Path
import json,hashlib,csv,math
import numpy as np
from scipy.signal import fftconvolve
from scipy.optimize import brentq
ROOT=Path('/home/tkwa/.codex/worktrees/protium-fusion/website');COMP=ROOT/'docs/protium-fusion/component-design';G=COMP/'catcher-followup/geometry';P=Path(__file__).resolve().parent;read=lambda p:json.loads(Path(p).read_text());sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
h=read(G/'temperature-histories.json');hp=G/'temperature_histories.py';src=hp.read_text().split('coarse_dt=60.',1)[0];ns={'__file__':str(hp),'__name__':'heatspread_readonly'};exec(compile(src,str(hp),'exec'),ns)
L=.175;rho=2000.;cp=1800.;sigma=350.;CA=sigma*cp;dt=.25;Rskin=.00015/4;Rpost=.00005/(50*.25);Rtot=Rskin+Rpost;qmean=ns['qmean'];edges=ns['edges'];times=ns['time'];curves=ns['curves'];basefloor=h['model']['initial_bulk_offset_K'];coarse=np.load(h['coarse_history']['npz']);bulk0=coarse['bulk_K'];profiles=[]
for pr in h['profiles']:
 old=np.genfromtxt(pr['csv'],delimiter=',',names=True);assert sha(pr['csv'])==pr['sha256'];tt=old['relative_to_prompt_arrival_s'];ab=pr['prompt_arrival_s']+tt;b,n,m,de=ns['evaluate'](ab);profiles.append((pr,tt,ab,b,m,de,old))
# Sourceheat and uniformstoredenergy are frozen beforechanging anykernels.
qstage={name:(curves[name]['q_material'][1:]+curves[name]['q_decay'][1:]).copy() for name in curves};n=np.arange(1,600);results_cache={};checks=[]

def evaluate_k(k,save=False):
 key=float(k)
 if key in results_cache and not save:return results_cache[key]
 if math.isinf(k):ker=np.zeros(len(times));hpulse=np.zeros(len(times));floor=basefloor-qmean*L/(2*50);offset=0.
 else:
  ker=np.diff(ns['ns']['front_step'](edges,L,k,cp,'N'));lam=k/(rho*cp)*(n*np.pi/L)**2;hpulse=np.zeros_like(times)
  for i in range(0,len(times),4096):hpulse[i:i+4096]=2/CA*np.sum(np.exp(-lam[:,None]*times[i:i+4096]),axis=0)
  floor=basefloor-qmean*L/(2*50)+qmean*L/(2*k);offset=qmean*L/(6*k)
 for name in curves:curves[name]['nonuniform']=np.r_[0.,fftconvolve(qstage[name],ker)[:len(times)]+ns['stages'][name]['prompt_fluence']*hpulse]
 rows=[]
 for pr,tt,ab,b,m,de,old in profiles:
  bb,nn,mm,dd=ns['evaluate'](ab);assert max(abs(bb-b))<1e-10 and max(abs(mm-m))<1e-8 and max(abs(dd-de))<1e-8
  bulk=b+floor;body=bulk+nn+offset;T=body+(m+de)*Rtot;im=int(np.argmax(T));mask=tt>=.25
  row={'stage':pr['stage'],'selection':pr['selection'],'max_front_material_phase_K':float(max(T[mask])),'max_front_time_s':float(tt[im]),'body_front_at_peak_K':float(body[im]),'bulk_at_front_peak_K':float(bulk[im]),'skin_and_post_rise_at_peak_K':float((m[im]+de[im])*Rtot),'max_bulk_K':float(max(bulk)),'max_material_plus_decay_flux_W_m2':float(max(m+de)),'duration_above_K_s':{str(th):float(sum(T[mask]>th)*dt) for th in [1200,1250,1300,1400,1500]},'scope':'Material/decay andresolvedprompttail; microsecondprompt interval is separate.'}
  rows.append(row)
  if save:
   f=P/f"k{'inf' if math.isinf(k) else int(k)}-{pr['stage']}-{pr['selection']}.csv"
   with f.open('w') as fp:
    w=csv.writer(fp);w.writerow(['relative_s','bulk_K','body_front_K','skin_plus_posts_front_K','material_plus_decay_W_m2']);w.writerows(zip(tt,bulk,body,T,m+de))
   row['history_file']=str(f);row['history_sha256']=sha(f)
  if k==50:checks.append({'name':'reference_'+pr['stage']+'/'+pr['selection'],'max_body_temperature_replay_abs_K':float(max(abs(body-old['front_bulk_k50_K'])))});assert max(abs(body-old['front_bulk_k50_K']))<1e-8
 maxrow=max(rows,key=lambda x:x['max_front_material_phase_K']);bulkmax=float(max(bulk0)+(floor-basefloor));out={'bulk_effective_k_W_m_K':None if math.isinf(k) else k,'isothermal_limit':math.isinf(k),'through_depth_conductance_W_m2_K':None if math.isinf(k) else k/L,'nonuniform_fundamental_decay_s':0. if math.isinf(k) else L*L*rho*cp/(math.pi**2*k),'new_initial_bulk_offset_K':floor,'maximum_sampled_whole_schedule_bulk_K':bulkmax,'maximum_material_front_over_profiles_K':maxrow['max_front_material_phase_K'],'controlling_profile':maxrow['stage']+'/'+maxrow['selection'],'profiles':rows};results_cache[key]=out;return out

rows=[evaluate_k(k,save=True) for k in [50.,100.,200.,500.,1000.,5000.,float('inf')]]
limits=rows[-1];screens=[]
for target in [1200.,1250.]:
 if limits['maximum_material_front_over_profiles_K']>=target:
  screens.append({'chosen_screen_K':target,'finite_bulk_k_solution':False,'isothermal_limit_front_K':limits['maximum_material_front_over_profiles_K'],'reason':'Fixedskin/postresistanceanduniformstoragealreadyexceedthechosencomparison.'})
 else:
  k=math.exp(brentq(lambda x:evaluate_k(math.exp(x))['maximum_material_front_over_profiles_K']-target,math.log(50),math.log(1e7),xtol=1e-8));r=evaluate_k(k);screens.append({'chosen_screen_K':target,'finite_bulk_k_solution':True,'required_effective_k_W_m_K':k,'through_depth_conductance_W_m2_K':k/L,'controlling_profile':r['controlling_profile'],'front_K':r['maximum_material_front_over_profiles_K'],'fundamental_decay_s':r['nonuniform_fundamental_decay_s'],'status':'Chosencomparison,notfailuretemperature;allotherparametersfixed.'})
# Peak heat deliveredto a hypotheticalfront spreader; its condensing/storage
# side can remain buffered. Do not demand instantaneousheatenginepower.
flux=[]
for name,s in ns['stages'].items():
 q=qstage[name];peak=float(max(q));fluence=s['fluence'];flux.append({'stage':name,'material_plus_delayed_peak_W_m2':peak,'total_prompt_plus_material_plus_delayed_fluence_J_m2':fluence,'constant_cycle_mean_contribution_W_m2':fluence/s['period'],'evaporation_mass_flux_kg_m2_s_for_latent_Jkg':{str(lat):peak/lat for lat in [1e6,4e6,20e6]},'phase_moved_kg_m2_if_all_stageheat_transported_for_latent_Jkg':{str(lat):fluence/lat for lat in [1e6,4e6,20e6]},'scope':'Latentheatsaredeclaredenergy-scale sensitivities,notassignedNa/Liproperties. A distributedspreader needtransport onlythenonuniformpart; thesearefullincomingheatrequirements.'})
# SurfaceR screen at infinitebulk k, preservingthe actualtimecorrelation.
Rlimits=[]
infloor=limits['new_initial_bulk_offset_K']
for target in [1200.,1250.]:
 req=[]
 for pr,tt,ab,b,m,de,old in profiles:
  pos=(m+de)>0;allow=(target-(b+infloor))/(m+de+1e-300);ii=int(np.argmin(np.where(pos,allow,np.inf)));req.append({'stage':pr['stage'],'selection':pr['selection'],'maximum_total_skin_post_resistance_m2_K_W':float(allow[ii]),'controlling_time_s':float(tt[ii]),'bulk_K':float(b[ii]+infloor),'flux_W_m2':float(m[ii]+de[ii])})
 Rlimits.append({'chosen_screen_K':target,'maximum_total_resistance_m2_K_W':min(x['maximum_total_skin_post_resistance_m2_K_W'] for x in req),'actual_total_resistance_m2_K_W':Rtot,'profiles':req})
assert all(rows[i]['maximum_material_front_over_profiles_K']>=rows[i+1]['maximum_material_front_over_profiles_K'] for i in range(len(rows)-1));checks.append({'name':'monotonic_k_curve','passed':True})
out={'scope':__doc__,'code_sha256':sha(__file__),'input_sha256':{str(p):sha(p) for p in [G/'temperature-histories.json',hp,COMP/'source128_thin_thermal.py',COMP/'receiver_thickness.py',COMP/'final-receiver-thermal-input.json',COMP/'final-family-dispatch/sigma350/events.csv',COMP.parent/'mixture-design/selected-plant/point.json']},'parameters':{'sigma_kg_m2':sigma,'thickness_m':L,'rho_kg_m3':rho,'cp_J_kg_K':cp,'dt_s':dt,'mean_pickup_W_m2':qmean,'skin_depth_m':.00015,'skin_k_W_m_K':4.,'skin_resistance_m2_K_W':Rskin,'post_resistance_m2_K_W':Rpost,'total_skin_post_resistance_m2_K_W':Rtot,'min_pickup_temperature_K':1000.,'whole_schedule_events':4228,'background_nonuniform_cutoff_s':20000.},'conductivity_cases':rows,'chosen_screens':screens,'infinite_bulk_skin_resistance_screens':Rlimits,'spreader_energy_flux_requirements':flux,'checks':checks,'checks_passed':len(checks),'limits':['Effectivek is a bulkresponse requirement,notatensorpropertyclaimforordinaryirradiatedgraphite.','Uniformstoredenergy,meanheatpickupandmassremainfixed; a heatpipeorphase-change networkmustfititsownmassandworking-fluidinventoryintotheplant.','Fixed150umskinandpostsretain theirseriesresistance; increasingonlybulkconductivitycannotremoveit.','Microsecondprompt finite-depth absorption is not recomputed andnot includedin1200/1250material-phasecomparisons.','SevenactualrepresentativeprofilesincludeallC13eventsandthewarmestA/B1/B2;thisisnotanexhaustivefine-timeintegrationofall4228events.','The entireincomingheat isretainedforthis conservative thermalcomparison; no extraescapingradiationcredit.','Selectedconstantk/cp,contactresistance,irradiationdamage,capillary/sonic/boiling limitsandmechanicalcycling requirecomponentqualification.']};(P/'results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps({'curve':[{k:r[k] for k in ['bulk_effective_k_W_m_K','isothermal_limit','maximum_material_front_over_profiles_K','maximum_sampled_whole_schedule_bulk_K','controlling_profile']} for r in rows],'screens':screens,'Rlimits':Rlimits,'flux':flux,'checks':len(checks)},indent=2))
