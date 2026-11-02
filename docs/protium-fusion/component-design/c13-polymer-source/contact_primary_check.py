"""Primary CH2 Hugoniot comparison for the isotope-specific C13 polymer source."""
from pathlib import Path
import json,sys,hashlib
P=Path(__file__).resolve().parent
sys.path.insert(0,str(P.parent))
from launcher_screen import Table
src=Path('/home/tkwa/src/protium-fusion/output/c13-polymer-source/C13-H2-source128-carbon-median.inputs.json')
d=json.loads(src.read_text());s=d['contact'];t=Table('C14H28_EOS_09-18-20.txt',15)
f=t.shock(s['q_J_kg']);assert f['success']
out={'scope':'Independent primary CH2 table with C13H2 inertial mass15u. Same molecular volume as naturalCH2, matching chosen cold density to table precision. No isotope-dependent electronic EOS or3D trajectory.','FEOS':s,'FPEOS':f,'FEOS_over_FPEOS_T':s['T_K']/f['temperature_K'],'FEOS_over_FPEOS_compression':(s['rho_kg_m3']/d['cold_reference']['rho_kg_m3'])/f['compression'],'FPEOS_vs_chosen_initial_density_relative':f['rho_initial_kg_m3']/d['cold_reference']['rho_kg_m3']-1,'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'table_sha256':hashlib.sha256(t.path.read_bytes()).hexdigest()}
(P/'contact-primary-check.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
