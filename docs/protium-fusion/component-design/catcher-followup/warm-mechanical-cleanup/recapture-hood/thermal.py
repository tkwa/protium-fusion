"""Conditional honeycomb tip: finite100mm graphite fin and actual pulse history.
One-dimensional axial heat conduction, full nominaltip energy retention,
insulated sidefaces (no radiativecooling credit), bottom on actualreceiverfront.
"""
from pathlib import Path
import json,hashlib,csv
import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.signal import lfilter
ROOT=Path('/home/tkwa/.codex/worktrees/protium-fusion/website');COMP=ROOT/'docs/protium-fusion/component-design';P=Path(__file__).resolve().parent;G=COMP/'catcher-followup/geometry';read=lambda p:json.loads(Path(p).read_text());sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
h=read(G/'temperature-histories.json');hp=G/'temperature_histories.py';prefix=hp.read_text().split('coarse_dt=60.',1)[0];ns={'__file__':str(hp),'__name__':'hood_readonly_history'};exec(compile(prefix,str(hp),'exec'),ns)
H=.1;rho=2000.;cp=1800.;dt=.25;floor=h['model']['initial_bulk_offset_K'];offset=h['model']['steady_back_extraction_front_offset_K'];checks=[]

def solve(q,B,k,N):
 dx=H/N;diff=k/(rho*cp*dx*dx);diag=np.full(N,-2*diff);diag[0]=-diff;diag[-1]=-3*diff;off=np.full(N-1,diff);lam,U=eigh_tridiagonal(diag,off);r=np.exp(lam*dt);fac=np.expm1(lam*dt)/lam
 vq=U[0,:]/(rho*cp*dx);vB=U[-1,:]*2*diff;x=(np.arange(N)+.5)*dx;init=B[0]+q[0]*(H-x)/k;z0=U.T@init;tip=np.full(len(q),0.);baseflux_integral=0;finalenergy=0.;wallinitenergy=sum(init)*rho*cp*dx
 # Only the tip row and integrated/baseflux facts need reconstruction.
 endpoint_last=0.;modalintegral=np.zeros(N)
 for j in range(N):
  force=vq[j]*q[:-1]+vB[j]*B[:-1];zz,_=lfilter([fac[j]],[1.,-r[j]],force,zi=[r[j]*z0[j]])
  zs=np.r_[z0[j],zz];tip+=U[0,j]*zs;endpoint_last+=U[-1,j]*zs[-1]
  modalintegral[j]=fac[j]*sum(zs[:-1])+(fac[j]-dt)/lam[j]*sum(force)
  finalenergy+=sum(U[:,j])*zs[-1]*rho*cp*dx
 integrated_bottom=float(U[-1,:]@modalintegral);integrated_B=sum(B[:-1])*dt;baseheat=2*k/dx*(integrated_bottom-integrated_B);Ein=sum(q[:-1])*dt;res=(finalenergy-wallinitenergy)+baseheat-Ein
 checks.append({'name':f'fin_energy_k{k}_N{N}','relative_residual':float(res/max(1,abs(Ein)))});assert abs(res)<1e-6*max(1,abs(Ein))
 return tip+q*dx/(2*k)

rows=[];profileout=[];selected={'A':'n14','B1':'n15','B2':'c12','C13':'c13'}
for pr in h['profiles']:
 t=np.arange(-20000.,4000.+dt,dt);absolute=pr['prompt_arrival_s']+t;b,f,m,de=ns['evaluate'](absolute);B=b+floor+f+offset;q=m+de
 # Tip itself gets a separate promptenergy impulse. The microsecond maximum
 # is evaluated withthe savedfinite-depthTeffproxy, not the timebin peak.
 for eventtime,name in zip(ns['et'],ns['enames']):
  ix=int(np.searchsorted(absolute,eventtime))
  if 0<=ix<len(q):q[ix]+=ns['stages'][name]['prompt_fluence']/dt
 for k,N in [(50.,100),(50.,200)]:
  tip=solve(q,B,k,N);mask=t>=-600;qt=q.copy();nonprompt=m+de;Tdam=tip+nonprompt*(1e-4/4-1e-4/k);im=np.argmax(np.where(t>=.25,tip,-np.inf));ip=int(np.searchsorted(t,0.))-1
  row={'stage':pr['stage'],'selection':pr['selection'],'wall_axial_k_W_m_K':k,'cells':N,'dx_m':H/N,'maximum_tip_material_phase_K':float(tip[im]),'peak_relative_time_s':float(t[im]),'receiver_front_at_tip_peak_K':float(B[im]),'incident_flux_at_tip_peak_W_m2':float(nonprompt[im]),'maximum_tip_plus100um_k4_damage_resistance_K':float(max(Tdam[t>=.25])),'tip_pre_prompt_K':float(tip[ip]),'prompt_Teff_microsecond_tip_peak_K_proxy':float(tip[ip]+pr['prompt_interval_upper']['upper_increment_bulk_K']),'prompt_duration_s':pr['prompt_interval_upper']['duration_s'],'max_base_front_saved_window_K':float(max(B[mask]))}
  rows.append(row)
  if N==200:
   # Can be used by the yield agent without assuming onepeaklasts allpulse.
   ff=P/f"tip-{pr['stage']}-{pr['selection']}.csv"
   with ff.open('w') as fp:
    w=csv.writer(fp);w.writerow(['relative_s','tip_K_k50','tip_K_100um_damage_k4','base_front_K','material_plus_decay_W_m2']);w.writerows(zip(t[mask],tip[mask],Tdam[mask],B[mask],nonprompt[mask]))
   profileout.append({'stage':pr['stage'],'selection':pr['selection'],'file':str(ff),'sha256':sha(ff)})
for pr in h['profiles']:
 rr=[r for r in rows if r['stage']==pr['stage'] and r['selection']==pr['selection']];err=rr[1]['maximum_tip_material_phase_K']-rr[0]['maximum_tip_material_phase_K'];checks.append({'name':'mesh:'+pr['stage']+'/'+pr['selection'],'peak_K_difference_200minus100':err});assert abs(err)<2
out={'scope':__doc__,'code_sha256':sha(__file__),'input_sha256':{str(p):sha(p) for p in [G/'temperature-histories.json',hp,COMP/'source128_thin_thermal.py',COMP/'final-receiver-thermal-input.json',COMP/'final-receiver-prompt-proxy.json']},'parameters':{'depth_m':H,'rho_kg_m3':rho,'cp_J_kg_K':cp,'k_W_m_K':50.,'dt_s':dt,'history_start_relative_s':-20000.,'late_end_relative_s':4000.,'bottom_boundary':'Actual fullyretained-source receiverfrontT(t) withconstantmeanbackpickup; hoodbottomloadsare smallandnotfed back into the existingreceiver solution.','side_boundary':'Insulated: no side radiation credit; lateralrecoil/chargedsecondary energy escape can reduce actualtipretention.','tip_source':'Full incidentmaterial kinetic+advectedinternal+delayedenergy assigned to projectedtiparea, pluspromptenergy. No1/tipfraction fluxamplification.','initial_condition':'Steady axialprofile with actual initialbaseT andinitialtipq,20000s beforetargetevent. Slowestfinmode~292s, so initialuncertainty is exponentiallyerased.'},'rows':rows,'profiles':profileout,'checks':checks,'limits':['This one-dimensionalfin model doesnotpredict thinwalllateraliontransport, anisotropicconductivity, joints, thermalwarping or temperaturedependentproperties.','All tips experience the same sourcefluences perprojectedarea as the floor. More side interception frommisalignment needsitsownheatdistribution.','The prompt microsecondtemperature uses theexistingfinite-depthTeffproxy andnot a newphoton/electron transportcalculation.','Temperature profiles are conditional calculations, not materialallowables or measuredhotwallperformance.']}
(P/'thermal.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(rows,indent=2))
