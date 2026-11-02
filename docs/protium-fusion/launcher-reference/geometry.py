"""Exact/diagnostic geometry for provisional launcher counts. No hardware optimization.
python geometry.py --root /path/to/worktree [--output geometry.json]
"""
import argparse,json,math,hashlib
from pathlib import Path
R=1.1e9;C=299792458.;SIGMA=5000.
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def calculate(root):
 p=root/'docs/protium-fusion/open-loop-redesign/station-inputs.json';d=json.loads(p.read_text());loss=d['cycle']['carbon_loss_per_cold_stage_output'];capture=d['cycle']['energy_capture_fraction']
 ip=root/'docs/protium-fusion/open-loop-redesign/collector-inputs.json';i=json.loads(ip.read_text());cp=i['stages']['c13']['record'];area=4*math.pi*R*R
 counts=[];apertures=[];budgets=[];fixed=[]
 for N in [100,1000,10000]:
  solid=4*math.pi/N;A=area/N;theta=math.acos(1-2/N)
  counts.append({'N':N,'solid_angle_cell_sr':solid,'cell_surface_area_m2':A,'equal_area_flat_disk_radius_m':math.sqrt(A/math.pi),'equal_area_cap_half_angle_rad':theta,'cap_geodesic_radius_m':R*theta,'characteristic_spacing_sqrt_area_m':math.sqrt(A),'spacing_scope':'sqrt(area per station); not a solved spherical packing or nearest-neighbor distance'})
  for akm in [1,10,100,1000]:
   a=akm*1000;A=N*math.pi*a*a;f=A/area
   # Exact solidangle of planar circular disks tangent toreceiver atdistanceR.
   u=math.sqrt(1+(a/R)**2);exact=N/2*(a/R)**2/(u*(u+1))
   apertures.append({'N':N,'aperture_radius_km':akm,'total_aperture_area_m2':A,'point_source_fraction_small_angle':f,'point_source_fraction_exact_planar_disks':exact,'fraction_of_energy_loss_budget':f/(1-capture),'fraction_of_fuel_loss_budget':f/loss,'local_C13_prompt_J_per_port':cp['radiation_escaped_J']/area*math.pi*a*a,'local_C13_kinetic_J_per_port':cp['final_kinetic_J']/area*math.pi*a*a,'equivalent_receiver_plug_mass_total_kg':SIGMA*A})
  for name,f in [('full_energy_loss_budget',1-capture),('full_fuel_loss_budget',loss)]:
   a=2*R*math.sqrt(f/N);exact=R*math.sqrt((1-2*f/N)**-2-1)
   budgets.append({'N':N,'criterion':name,'fraction':f,'aperture_radius_small_angle_m':a,'aperture_radius_exact_planar_disks_m':exact})
  # Hold the sumarea equal to1000 disks of1kmradius.
  Atot=1000*math.pi*1000**2;a=math.sqrt(Atot/(math.pi*N))
  fixed.append({'N':N,'total_aperture_area_m2':Atot,'each_radius_m':a,'loss_fraction':Atot/area,'total_circumference_m':2*math.pi*N*a,'equivalent_plug_mass_total_kg':SIGMA*Atot})
 timing=[]
 for k,it in i['stages'].items():
  rec=it['record'];vi=max(map(abs,it['initial']['velocity_m_s']));vo=max(rec['profile']['velocity_m_s'])
  timing.append({'stage':rec['stage'],'fastest_incoming_m_s':vi,'last_radial_ballistic_launch_before_focus_s':R/vi,'first_endpoint_ballistic_ion_after_focus_s':R/vo,'time_from_last_launch_to_prompt_arrival_s':R/vi+R/C,'time_from_last_launch_to_first_ion_s':R/vi+R/vo})
 gates=[]
 T=250.
 for akm in [1,10,100,1000]:
  a=akm*1000;m=SIGMA*math.pi*a*a;v=4*a/T;acc=8*a/(T*T)
  gates.append({'radius_km':akm,'closure_time_s':T,'stroke_m':2*a,'triangle_velocity_profile_peak_m_s':v,'acceleration_m_s2':acc,'plug_mass_kg':m,'peak_translation_kinetic_J':m*v*v/2,'scope':'Monolithic full-diameter rest-to-rest sliding plug illustration; segmentation/rotation can alter mechanics; kinetic energy may be recovered.'})
 atten=[{'photon_energy_MeV':E,'mu_m2_kg':mu,'sigma_kg_m2':SIGMA,'uncollided_fraction':math.exp(-mu*SIGMA),'sigma_for_1e_6_uncollided_kg_m2':-math.log(1e-6)/mu} for E,mu in [(10,.001959),(20,.001575)]]
 checks={'small_angle_agreement_max_relative':max(abs(x['point_source_fraction_exact_planar_disks']/x['point_source_fraction_small_angle']-1) for x in apertures),'fixed_total_area_fraction_spread':max(x['loss_fraction'] for x in fixed)-min(x['loss_fraction'] for x in fixed)}
 assert checks['small_angle_agreement_max_relative']<.0002;assert checks['fixed_total_area_fraction_spread']==0
 return {'scope':'Exact arithmetic for stated geometry and deterministic budget diagnostics. No inferred credible interval or validation of 1000 launchers. Each launcher may emit multiple packets.','radius_m':R,'receiver_area_m2':area,'energy_loss_budget':1-capture,'fuel_loss_budget':loss,'total_C_H_budget':d['cycle']['total_ongoing_C_H_mass_budget'],'fuel_C_H_budget':d['cycle']['fresh_C_H_mass_ratio'],'prompt_arrival_s':R/C,'C13_pulse_mass_kg':cp['mass_kg'],'C13_prompt_fluence_J_m2':cp['radiation_escaped_J']/area,'C13_kinetic_fluence_J_m2':cp['final_kinetic_J']/area,'counts':counts,'fixed_radius_apertures':apertures,'full_budget_diagnostic_radii':budgets,'fixed_total_aperture_area':fixed,'ballistic_closure_timing':timing,'monolithic_sliding_plug_examples':gates,'carbon_shield_uncollided_diagnostic':atten,'checks':checks,'input_hashes':{'station-inputs.json':sha(p),'collector-inputs.json':sha(ip)},'code_sha256':sha(Path(__file__))}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,default=Path(__file__).with_name('geometry.json'));a=p.parse_args();o=calculate(a.root);a.output.write_text(json.dumps(o,indent=2)+'\n');print(json.dumps({k:o[k] for k in ['counts','full_budget_diagnostic_radii','fixed_total_aperture_area','checks']},indent=2))
