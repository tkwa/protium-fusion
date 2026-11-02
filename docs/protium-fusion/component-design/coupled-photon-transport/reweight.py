"""Reweight a real Geant4 monochromatic response basis; standard library only.
All uncertainty intervals here are Monte Carlo sampling intervals, not model CIs.
"""
import math
Z80=1.2815515655446004

def _record(basis,sigma,particle,mode,E):
 rows=basis['columns'][str(int(sigma))]['responses']
 matches=[v for v in rows if v['particle']==particle and v['mode']==mode and abs(v['energy_MeV']-E)<1e-8]
 if len(matches)!=1:raise ValueError(f'No unique response for {sigma},{particle},{mode},{E}')
 return matches[0]

def evaluate(basis,spectrum,sigma,scenario='broad_D1p10',prompt_policy='EM20',cold_intercept_fraction=1.,hot_bypass_ppm=0.):
 if scenario not in ['nominal_lines_hardbeta','broad_D1p10']:raise ValueError('Explicit supported source scenario required')
 if prompt_policy not in ['EM20','all_to_cold']:raise ValueError('Explicit prompt policy required')
 if not 0<=cold_intercept_fraction<=1 or hot_bypass_ppm<0:raise ValueError('Invalid partition/reserve')
 H=spectrum['thermal_boundary_MeV_per_He'];lines=spectrum['primary_photon_lines'];beta=spectrum['main_beta_kinetic_MeV_per_He'];D=1.1 if scenario=='broad_D1p10' else 1.
 anni=sum(v['MeV_per_He'] for v in lines if abs(v['energy_MeV']-.511)<1e-8);main=sum(v['MeV_per_He'] for v in lines if abs(v['energy_MeV']-2.312593)<1e-8);minor=sum(v['MeV_per_He'] for v in lines if all(abs(v['energy_MeV']-e)>1e-8 for e in [.511,2.312593]))
 weights={2.5438523:D*main,5.657344:D*(anni+beta)} if D>1 else {.511:anni,2.312593:main,5.14304:beta}
 components=[];variance=0.;leak=0.
 for E,W in weights.items():
  r=_record(basis,sigma,'gamma','isotropic',E);q=W*r['leakage_fraction'];leak+=q;variance+=(W*r['sampling_standard_error'])**2;components.append({'source_energy_MeV_per_He':W,'primary_gamma_MeV':E,'leakage_MeV_per_He':q,'kernel':'isotropic25um'})
 prompt=spectrum['cycle']['escaped_prompt_MeV_per_He'];p=_record(basis,sigma,'gamma','radial',20.)
 prompt_leak=prompt if prompt_policy=='all_to_cold' else prompt*p['leakage_fraction']
 if prompt_policy=='EM20':variance+=(prompt*p['sampling_standard_error'])**2
 reserve=D*(minor+spectrum['other_unstable_whole_Q_upper_MeV_per_He'])+hot_bypass_ppm*1e-6*H;nonprompt=leak+reserve;leak=nonprompt+prompt_leak;f=leak/H;se=math.sqrt(variance)/H
 cold=(cold_intercept_fraction*nonprompt+(prompt_leak if prompt_policy=='all_to_cold' else cold_intercept_fraction*prompt_leak))/H
 return {'sigma_kg_m2':sigma,'scenario':scenario,'prompt_policy':prompt_policy,'thermal_boundary_MeV_per_He':H,'leakage_fraction':f,'hot_capture_fraction':1-f,'sampling_standard_error':se,'sampling_80pct_half_width':Z80*se,'hot_capture_sampling_80pct_interval':[max(0,1-f-Z80*se),min(1,1-f+Z80*se)],'components':components,'saved_prompt_source_MeV_per_He':prompt,'prompt_leakage_MeV_per_He':prompt_leak,'untransported_minor_lines_and_other_full_loss_MeV_per_He':D*(minor+spectrum['other_unstable_whole_Q_upper_MeV_per_He']),'carried_hot_bypass_MeV_per_He':hot_bypass_ppm*1e-6*H,'cold_equipment_fraction':cold,'external_escape_fraction':f-cold,'cold_intercept_fraction':cold_intercept_fraction,'source_uplift_drawn_from_existing_bulk_MeV_per_He':(D-1)*(sum(v['MeV_per_He'] for v in lines)+beta),'scope':'Conditional source-response estimate, not an energy-maximum theorem or full-device capture guarantee. Natural graphite; closed spherical cavity; two-photon in-flight annihilation approximation. Missing photonuclear physics is isolated by prompt-policy sensitivity. No extra nuclear energy is created by D uplift.'}

def positron_response_grid(basis,spectrum,sigma):
 """Energy-normalized comparisons, preserving the total beta+annihilation budget.
 Effective particle count varies with K; these are response brackets, not claimed
 physical beta spectra or actual nuclear parent-count replacements.
 """
 H=spectrum['thermal_boundary_MeV_per_He'];W=spectrum['main_beta_kinetic_MeV_per_He']+sum(v['MeV_per_He'] for v in spectrum['primary_photon_lines'] if abs(v['energy_MeV']-.511)<1e-8)
 base=evaluate(basis,spectrum,sigma,'nominal_lines_hardbeta');old=sum(v['leakage_MeV_per_He'] for v in base['components'] if v['primary_gamma_MeV'] in [.511,5.14304]);rows=[]
 for K in [.5,1.,1.8,4.12]:
  r=_record(basis,sigma,'e+','radial',K);q=W*r['leakage_fraction'];rows.append({'positron_kinetic_energy_MeV':K,'calorimetric_energy_per_positron_MeV':r['calorimetric_energy_MeV'],'beta_plus_annihilation_budget_MeV_per_He':W,'effective_positrons_per_He_for_energy_normalized_comparison':W/r['calorimetric_energy_MeV'],'response_leakage_fraction_of_positron_budget':r['leakage_fraction'],'response_sampling_80pct_half_width':Z80*r['sampling_standard_error'],'whole_thermal_leakage_fraction_with_this_response':base['leakage_fraction']+(q-old)/H})
 return {'sigma_kg_m2':sigma,'rows':rows,'scope':positron_response_grid.__doc__}
