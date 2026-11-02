"""Combine separated slow and prompt thermal requirements without simultaneous peaks.
Uses actual source128 profiles and explicit finite-depth photon proxy alternatives.
"""
from pathlib import Path
import json,hashlib,math
P=Path(__file__).parent;sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();a=P/'source128-thin-thermal.json';b=P/'source128-prompt-proxy.json';t=json.loads(a.read_text());p=json.loads(b.read_text());cpath=P/'source128-prompt-radiance.json';rad=json.loads(cpath.read_text());rows=[]
for r in t['thermal_rows']:
 if r['k_W_m_K']!=50:continue
 for skin_k in[50.,4.]:
  for proxy in['outer_T','surface_flux_equivalent_Teff','outer_T_quarter']:
   stages={}
   for name,v in r['isolated_peaks'].items():
    m=p['stages'][name]['cases'][proxy]['materials']['carbon'];depth=m['resolved_surface_no_conduction_sensible_rise_K']/1.8;soft=m['below1keV_all_surface_conduction_peak_K']*math.sqrt(100000/(skin_k*1800));allface=m['all_photons_surface_conduction_peak_K']*math.sqrt(100000/(skin_k*1800));flash=min(depth+soft,allface)
    # Full added skin resistance on K+U plus maximum possible decay envelope.
    # This deliberately places even the more penetrating photons at the skin.
    rates={'n13':.001159113256783153,'o14':.009817978825996571,'o15':.005681556796002065};decay_q_upper=sum(E*rates[n]for n,E in t['pulse_fluences'][name]['decay_fluence'].items());skinextra=(v['maximum_K_plus_U_flux_W_m2']+decay_q_upper)*.0001/skin_k if skin_k==4 else 0.
    slow=v['slow_plus_prompt_tail_nonuniform_peak_K']+skinextra
    stages[name]={'prompt_depth_plus_soft_reference_K':depth+soft,'all_surface_prompt_reference_K':allface,'selected_prompt_reference_K':flash,'slow_reference_K':slow,'added100um_skin_resistance_K_on_slow_flux_upper':skinextra,'prompt_fraction_energy_requiring_W_gt1':rad['stages'][name]['cases'][proxy]['fraction_prompt_energy_in_bins_W_gt1'],'max_of_noncoincident_prompt_slow_K':max(flash,slow)}
   before=1000+r['common_background_upper_K']+max(v['max_of_noncoincident_prompt_slow_K']for n,v in stages.items()if n!='B1')+stages['B1']['max_of_noncoincident_prompt_slow_K'];pickup=r['bulk_floor_offset_to_keep_pickup_at_least1000K'];jitter=r['additional_bulk_swing_for_bounded1day_dispatch_lateness_K']
   rows.append({'sigma_kg_m2':r['sigma_kg_m2'],'bulk_k_W_m_K':50,'cp_J_kg_K':1800,'skin_k_W_m_K':skin_k,'proxy':proxy,'physical_scope':('Thermally permitted diluted/blackbody shape at the saved source area; actual spectrum unresolved.'if proxy!='outer_T_quarter'else'Unconstrained shape sensitivity: some bins exceed blackbody radiance at assigned colour. Not an attainable thermal photosphere at saved area or an imposed design requirement.'),'stages':stages,'periodic_all_phase_front_reference_K_bulk_floor1000':before,'with_pickup_floor1000_K':before+pickup,'with_pickup_floor1000_and_bounded1day_dispatch_lateness_K':before+pickup+jitter})
out={'scope':__doc__,'code_sha256':sha(__file__),'input_sha256':{str(a):sha(a),str(b):sha(b),str(cpath):sha(cpath)},'rows':rows,'limits':['Choose max(prompt,slow) within each individual pulse; independent big-pulse and B1 maxima plus storage/history form a conservative all-phase envelope, not a solved event-by-event surface history.','Photon shapes and local energy transfer are conditional proxies, not measured frequency-dependent transport. Their alternatives are not credible intervals.','The k4 skin adds the full100um resistance on a bound to K+U+delayed surface heat flux without changing bulk conductivity. Contact impedance remains a separate qualification.','The lower of the first-interaction depth proxy and all-surface reference is used within the respective homogeneous photon model; secondary transport and layered heat flow remain unsolved.','All source heat is retained in the backing for storage; no thin-receiver capture fraction is silently set to99%.']}
(P/'source128-thermal-combined.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps([{k:v for k,v in r.items()if k!='stages'}for r in rows],indent=2))
