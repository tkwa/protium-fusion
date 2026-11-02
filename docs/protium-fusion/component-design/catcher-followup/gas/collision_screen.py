"""Neutral-gas Knudsen sensitivity; no Boltzmann collision solver."""
from pathlib import Path
import json,math,hashlib
P=Path(__file__).resolve().parent;k=1.380649e-23;T=1000.
i=json.loads((P/'inputs.json').read_text());R=i['R_m'];m=[2*i['atomic_mass_kg']['p'],i['atomic_mass_kg']['he4']]
v=[math.sqrt(8*k*T/(math.pi*x)) for x in m]
# Primary NIST Table A1, DOI10.6028/jres.105.054: eta(4He,1000K)=46.378 microPa.s.
eta=46.378e-6;sig_eta=5*math.sqrt(math.pi*m[1]*k*T)/(16*eta)
cases=[];hashes={}
for name,fn,g in [('full_injection_g1','chamber-dt150.0-q128.json',1.),('cap100nm_g0p01','chamber-dt150.0-q128-cap0.1_eta100000.0.json',.01),('cap10nm_g0p01','chamber-dt150.0-q128-cap0.01_eta100000.0.json',.01)]:
 d=json.loads((P/fn).read_text());hashes[fn]=hashlib.sha256((P/fn).read_bytes()).hexdigest()
 s=next(x for x in d['stationary'] if x['T_K']==1000 and x['H_atoms_per_gas_particle']==2)
 p=[x/g for x in s['unit_capture_pressures_Pa']];n=[x/(k*T) for x in p];rows=[]
 for sigma in [1e-20,3e-20,1e-19,sig_eta,3e-19,1e-18]:
  nu=[sum(n[j]*sigma*math.sqrt(8*k*T/(math.pi*(m[a]*m[j]/(m[a]+m[j])))) for j in range(2)) for a in range(2)]
  lam=[v[a]/nu[a] for a in range(2)]
  rows.append({'common_pair_effective_cross_section_m2':sigma,'collision_frequency_H2_He_s1':nu,'mean_free_path_H2_He_m':lam,'Knudsen_lambda_over_R_H2_He':[x/R for x in lam],'collisions_during_unit_capture_residence_H2_He':[nu[a]*s['unit_capture_residence_s'][a] for a in range(2)],'collisions_during_actual_capture_residence_H2_He':[nu[a]*s['unit_capture_residence_s'][a]/g for a in range(2)],'collision_time_H2_He_s':[1/x for x in nu]})
 r=next(x for x in d['cases'] if x['model']=='well_mixed_volume_balance' and x['T_K']==1000 and x['H_atoms_per_gas_particle']==2 and x['effective_wall_capture_g']==g and x['material_release_time_s']==0)
 cases.append({'case':name,'pressure_H2_He_Pa':p,'number_density_H2_He_m3':n,'total_density_m3':sum(n),'diffuse_chronology_peak_over_mean_column':r['peak_radial_column_kg_m2']/r['mean_radial_column_kg_m2'],'rows':rows})
out={'primary_anchor':{'source':'https://doi.org/10.6028/jres.105.054','fulltext':'https://pmc.ncbi.nlm.nih.gov/articles/PMC4872693/','table':'A1,4He at1000K','viscosity_Pa_s':eta,'hard_sphere_viscosity_equivalent_area_m2':sig_eta,'conversion':'sigma_eta=5 sqrt(pi m kT)/(16 eta), obtained by inverting the hard-sphere viscosity law. This is an order-of-magnitude transport-area diagnostic, not the true total scattering cross section or the H2-He mixture kernel.'},'cases':cases,'scope':['Pair relative Maxwell speeds and explicit common effective areas provide a sensitivity, not inferred mixture collision data.','Uniform steady Maxwellian molecular wall flux remains n vbar/4 in the presence of gas-gas collisions. The corresponding globally uniform inventory balance therefore does not require a collisionless cavity.','The spherical first-flight impulse response requires collisionless chords and is not an upper or lower bound when Kn is near or below unity. Collisions can change source peaks, angular focusing, and release-front relaxation.','The routed g=.01 case has a nearly steady reservoir in the diffuse model; this does not bound localized fresh plumes or spatially segregated pumping.','A full calculation would use species/state-dependent cross sections and the actual source/sink angular geometry in a kinetic transport solver.'],'source_sha256':hashes,'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(P/'collision-screen.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'primary_anchor':out['primary_anchor'],'nominal100nm':[r for r in cases[1]['rows'] if r['common_pair_effective_cross_section_m2']==1e-19]},indent=2))
