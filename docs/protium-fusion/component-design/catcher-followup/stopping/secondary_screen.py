"""Conditional C13(alpha,n) secondary ceiling, not a fitted activation network."""
from pathlib import Path
import json,math,hashlib
P=Path(__file__).resolve().parent;MU=1.66053906660e-27;EV=1.602176634e-19
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
s=json.loads((P/'surface-spectra.json').read_text());d=json.loads((P/'surface-stopping.json').read_text());rows=[]
for x13 in [.0107,1.]:
 for barns in [1.,10.]:
  stages=[]
  for stage,a in s['stages'].items():
   N=a['species']['he4']['packet_nuclei_m2'];range_m=[x*1e-6 for x in d['stages'][stage]['species']['he4']['packet_CSDA_depth_um']];amu=12.011 if x13<1 else 13.;prob=[-math.expm1(-x*2000/(amu*MU)*x13*barns*1e-28) for x in range_m];events=sum(n*p for n,p in zip(N,prob));energy=sum(n*p*(4*eu+2.215)*1e6*EV for n,p,eu in zip(N,prob,a['packet_energy_MeV_u']));stages.append({'stage':stage,'reaction_events_m2_pulse':events,'probability_per_incident_He':events/sum(N),'all_reaction_kinetic_energy_J_m2_pulse_upper':energy,'neutron_mass_kg_m2_pulse_if_all_escape':events*1.00866491595*MU,'period_s':a['period_s']})
  rows.append({'assumed_C13_atomic_fraction':x13,'assumed_cross_section_ceiling_barn':barns,'stages':stages,'mean_all_reaction_kinetic_W_m2_upper':sum(x['all_reaction_kinetic_energy_J_m2_pulse_upper']/x['period_s'] for x in stages),'mean_neutron_mass_kg_m2_s_if_all_escape':sum(x['neutron_mass_kg_m2_pulse_if_all_escape']/x['period_s'] for x in stages)})
out={'scope':__doc__,'surface_spectra_sha256':sha(P/'surface-spectra.json'),'surface_stopping_sha256':sha(P/'surface-stopping.json'),'code_sha256':sha(__file__),'primary':'https://doi.org/10.1103/PhysRevC.108.L061601','rows':rows,'qualification':'The constant1/10barn curves are deliberately high conditional ceilings on sigma(E), not evaluated measured cross-section functions or universal proven bounds.2023direct data2.9–8MeV lie below1barn; dominantB1alpha spectrum is mostly below2MeV and requires the lower-energy measurements for precise yield. The all-C13 row tests possible enrichment; other nuclear channels,neutron multiplication/transport and local dose remain outside this screen. Charged-particle CSDA alone does not imply zero neutron escape. Kinetic upper assigns each reaction its incident alpha energy plusQ2.215MeV, even though the neutron receives only part.','checks':{'probabilities_bounded':all(0<=st['probability_per_incident_He']<=1 for r in rows for st in r['stages']),'all_energies_nonnegative':all(st['all_reaction_kinetic_energy_J_m2_pulse_upper']>=0 for r in rows for st in r['stages'])}}
(P/'secondary-screen.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(rows,indent=2))
