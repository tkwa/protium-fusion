"""Primary commercial pump benchmarks plus explicitly conditional scale checks."""
from pathlib import Path
import json,math,hashlib
P=Path(__file__).resolve().parent
KB=1.380649e-23;MU=1.66053906660e-27
mass={'H2':2*1.00782503223*MU,'He':4.00260325413*MU}
bench=[dict(model='Pfeiffer HiPace700H',nominal_inlet_diameter_m=.160,speed_L_s={'H2':555.,'He':655.},K0={'H2':2e7,'He':1e9},rpm=49200.,mass_kg=11.8,source='https://www.pfeiffervacuum.com/global/assets/134975/1/Product%20Leaflet%20HiPace%20700%20H%20Global%20EN.pdf',manual='https://www.pfeiffervacuum.com/de/assets/123614/1/Instruction%20Manual%20HiPace%20700%20H%20with%20TC%20400_DE_de%2CEN_en.pdf'),
 dict(model='Edwards nEXT85H ISO63',nominal_inlet_diameter_m=.063,speed_L_s={'H2':54.,'He':78.},K0={'H2':5e5,'He':2e7},rpm=90000.,mass_kg=3.,source='https://www.edwardsvacuum.cn/content/dam/brands/edwards-vacuum/general-vacuum/downloads/turbomolecular-pumps/3601-0243-01-nEXT85-datasheet.pdf',rpm_source='https://digitalbrochure.edwardsvacuum.com/nEXT-TURBOMOLECULAR-PUMPS/en/technical-specifications.html')]
rows=[];checks=[]
for b in bench:
 area=math.pi*b['nominal_inlet_diameter_m']**2/4
 for sp in mass:
  S=b['speed_L_s'][sp]/1000;K=b['K0'][sp]
  for T in [300.,1000.]:
   v=math.sqrt(8*KB*T/(math.pi*mass[sp]));capture=4*S/(area*v)
   rows.append(dict(model=b['model'],species=sp,T_K=T,nominal_inlet_area_m2=area,
    mean_thermal_speed_m_s=v,nominal_speed_m3_s=S,effective_aperture_capture_if_same_speed=capture,
    aperture_fraction_for_global_g0p01_if_same_speed=.01/capture,
    scope='DN-derived nominal inlet area; actual bore/hub/baffle layout may differ.1000K row holds measured room-temperature speed fixed as a comparison, not a predicted hot-pump speed.'))
   checks.append(0<capture<1)
  hotK=math.exp(math.log(K)*math.sqrt(300/1000))
  rows.append(dict(model=b['model'],species=sp,all_hot_K0_if_lnK_scales_u_over_sqrtT=hotK,
   same_K0_required_speed_multiplier_300to1000=math.sqrt(1000/300),
   same_speed_required_log_compression_length_multiplier_300to1000=math.sqrt(1000/300),
   scope='Similar-geometry, fully molecular and all-stage1000K heuristic. A hot entrance followed by shielded cooler stages must be treated stagewise. No prediction for the commercial pump at1000K.'))
  fraction=1-1e-3/(K*1e-10)
  rows.append(dict(model=b['model'],species=sp,pressure_ratio_requested=1e7,
    zero_throughput_K0=K,signed_speed_fraction_at_Pin1e10minus_Pout1e3minus=fraction,
    feasible_forward_flow_at_those_partial_pressures=fraction>0,
    scope='Isothermal molecular Q=S0(Pin-Pout/K0), using separate species partial pressures. A negative result denotes backflow, not negative available hardware size.'))
required=[dict(retained_zero_backpressure_speed_fraction=f,K0_min_for_Pout_over_Pin1e7=1e7/(1-f)) for f in [.5,.9,.99]]
out={'scope':__doc__,'benchmarks':bench,'calculations':rows,'required_compression_margin':required,
 'tip_speed_primary_presentation':{'source':'https://indico.cern.ch/event/161398/contributions/235875/attachments/187008/262586/Turbopumps_from_theory_to_practice.pdf','author':'Florian Henss, Pfeiffer Vacuum R&D,2011-12-06,slide12','indexed_tip_speed_m_s':[300,400],'all_hot_similar_ratio_m_s':[300*math.sqrt(1000/300),400*math.sqrt(1000/300)],'access':'Search-indexed primary slide text read; direct PDF fetch returned403. This is an engineering reference, not a tested hot-rotor rating.'},
 'theory_sources':{'manufacturer':'https://www.pfeiffervacuum.com/no/en/knowledge/vacuum-technology/knowledge-book/4-vacuum-generation/4.9-turbomolecular-pumps.html','primary_DSMC':'https://arxiv.org/abs/physics/0401028'},
 'checks':{'capture_within0and1':all(checks),'row_count':len(rows)},'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
assert all(checks)
(P/'results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
