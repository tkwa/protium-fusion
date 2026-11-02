"""Stoichiometric collector capacity and regeneration work; no kinetic capture credit."""
from pathlib import Path
import json,hashlib
P=Path(__file__).resolve().parent;d=json.loads((P/'results.json').read_text());YEAR=365.25*86400;EV=1.602176634e-19;MU=1.66053906660e-27;MC=12*MU;AREA=d['constants']['area_m2'];NA=6.02214076e23
q1=sum(d['constant_yield_response_coefficients'][g]['C_kg_m2_year_per_Y1'] for g in ['C','N','O'])/YEAR
source={'Y_CNO_1':q1,**{k:v['C_kg_s']/AREA for k,v in d['model_diagnostics'].items()}}
rows=[]
for name,q in source.items():
 n=q/MC;r={'source':name,'C_kg_m2_year':q*YEAR,'C_atoms_m2_s':n,'electrical_work_W_m2_by_positive_process_eV_perC_at_eta0p5':{str(e):n*e*EV/.5 for e in [1,2,10,100]},'getters':[]}
 for metal,A,rho,G in [('Si',28.0855,2330.,{1000:-63.517,1500:-59.661}),('Ti',47.867,4500.,{1000:-172.653,1500:-165.821})]:
  qm=q*A/12
  for surfarea in [1.,40.]:
   r['getters'].append({'element':metal,'metal_atoms_per_C':1,'metal_kg_m2_receiver_year':qm*YEAR,'assumed_metal_density_kg_m3':rho,'collector_area_m2_per_receiver_m2':surfarea,'consumed_pure_metal_equivalent_thickness_nm_year':qm*YEAR/rho/surfarea*1e9,'100nm_initial_layer_capacity_years_if_all_atoms_accessible':100e-9*rho*surfarea/(qm*YEAR),'10nm_layer_inventory_kg_per_receiver_m2':10e-9*rho*surfarea,'10nm_layer_C_capacity_kg_per_receiver_m2':10e-9*rho*surfarea*12/A,'10nm_layer_capacity_years_if_all_atoms_accessible':10e-9*rho*surfarea/(qm*YEAR),'JANAF_carbide_regeneration_to_elements_eV_perC':{str(T):-g*1000/(NA*EV) for T,g in G.items()},'regeneration_W_m2_at50pct_efficiency':{str(T):n*(-g*1000/NA)/.5 for T,g in G.items()}})
 rows.append(r)
out={'scope':'Capacity and positive-work accounting for C capture in stoichiometric SiC orTiC. Protected collectors may be refreshed or regenerated; no sticking probability, access through passivating carbide, selective capture againstN/O/H, erosion lifetime, or nanosystem power overhead is inferred. Annual source mass is gross removal before local recapture; all values scale linearly with actual collected fraction and yield.','thermodynamic_sources':{'SiC':'https://janaf.nist.gov/tables/C-100.html','TiC':'https://janaf.nist.gov/tables/C-107.html'},'rows':rows,'input_sha256':hashlib.sha256((P/'results.json').read_bytes()).hexdigest(),'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(P/'getter-results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(rows[0],indent=2))
