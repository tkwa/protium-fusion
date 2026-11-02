"""Conditional retention and re-emission diagnostics; no reactor sticking posterior."""
import json,math,hashlib
from pathlib import Path
P=Path(__file__).resolve().parent;d=json.loads((P/'results.json').read_text());KB=8.617333262145e-5
rows=[]
for Edes in [1.4]:
 for T in [1000,1250,1500]:
  nu=1e13;Em=.47;a=.142e-9;tau=math.exp(Edes/(KB*T))/nu;hops=math.exp((Edes-Em)/(KB*T));rows.append({'T_K':T,'E_ads_assumed_desorption_barrier_eV':Edes,'migration_barrier_eV':Em,'attempt_rate_s':nu,'hop_length_m':a,'isolated_adatom_residence_s':tau,'expected_hops_before_desorption_equal_prefactors':hops,'RMS_distance_m_in_2D_random_walk':a*math.sqrt(hops),'scope':'Illustrative transition-state estimate using static clean-graphene DFT energies; no measured kinetic prefactor or damaged-surface trap network.'})
transport=[]
for name,m in d['model_diagnostics'].items():
 for E in [1.,3.]:
  initial=m['atomic_column_by_mean_energy_eV'][str(E)];thermal=m['atomic_column_by_mean_energy_eV'][str(2*1.380649e-23*1000/1.602176634e-19)]
  for s in [.1,.25,.5,1.]:transport.append({'source_diagnostic':name,'initial_mean_eV':E,'per_encounter_permanent_capture_probability':s,'steady_C_column_kg_m2':initial+(1-s)/s*thermal,'formula':'initial-flight column + (1-s)/s times1000K thermal single-flight column; cosine re-emission, no chemical conversion.'})
out={'scope':'Primary short-time chemisorption evidence plus explicitly assumed long-time capture probability. S is a scenario, not a calibrated irradiated-graphite property.','MD_amorphous_C_thermal_incidence':[{'T_K':1000,'P_chem':.247,'authors_sampling_error':.018},{'T_K':1250,'P_chem':.295,'authors_sampling_error':.016},{'T_K':1500,'P_chem':.299,'authors_sampling_error':.020}],'MD_source':'Bossion etal A&A692,A249(2024),TableA.1, https://arxiv.org/abs/2411.06125','basal_graphene_DFT_source':'Lehtinen etal PRL91,017202(2003), https://doi.org/10.1103/PhysRevLett.91.017202','isolated_adatom_kinetic_diagnostics':rows,'repeated_global_flight_scenarios':transport,'input_sha256':hashlib.sha256((P/'results.json').read_bytes()).hexdigest(),'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(P/'sticking-results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(rows,indent=2))
