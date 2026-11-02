"""All-charged-ion surface number flux on the exact saved quadrature interface."""
from pathlib import Path
import json,sys
import numpy as np
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P));import spectra as sp
rs,pm,st,decay,weak,hashes=sp.setup();s=json.loads((P/'surface-spectra.json').read_text());out={};checks=[]
for stage,row in s['stages'].items():
 x=sp.packets(Path(row['record_path']),Path(row['trace_path']),pm,decay,order=32);counts=x['n'][:,pm.ZZ>0].sum(axis=1);flux=counts/x['dvweight']*x['arrival_jacobian'];j=int(np.argmax(flux));N=float(counts.sum());expected=sum(row['whole_pulse_species_wall_nuclei_m2'][name] for name,z in zip(s['nuclear_species'],s['charges']) if z>0)
 fine=sp.packets(Path(row['record_path']),Path(row['trace_path']),pm,decay,order=64);nf=fine['n'][:,pm.ZZ>0].sum(axis=1);ff=nf/fine['dvweight']*fine['arrival_jacobian'];jf=int(np.argmax(ff));maxf=float(ff[jf]);checks.append({'stage':stage,'count_closure':bool(abs(N/expected-1)<1e-12),'ordered_interface':bool(np.array_equal(x['arrival'],row['packet_arrival_after_ignition_s']))})
 out[stage]={'incident_charged_nuclei_m2':N,'packet_charged_nuclei_m2':counts.tolist(),'packet_number_flux_m2_s':flux.tolist(),'peak_sampled_number_flux_m2_s':float(flux[j]),'peak_sampled_time_after_ignition_s':float(x['arrival'][j]),'quadrature64_peak_sampled_number_flux_m2_s':maxf,'quadrature64_peak_time_after_ignition_s':float(fine['arrival'][jf]),'quadrature32_to64_peak_relative_change':maxf/float(flux[j])-1}
for c in checks:assert all(c[k] for k in ['count_closure','ordered_interface'])
r={'scope':'Sum of all Z>0 modeled nuclei; charge state does not change nucleus number. Arrays match32-point surface-spectra packet order. Peaks are samples of the continuous within-cell arrival Jacobian, with64-point comparison; they are not a rigorous bound over unresolved hydrodynamic structure. No track-to-hole or leakage yield assigned.','surface_spectra_sha256':sp.sha(P/'surface-spectra.json'),'code_sha256':sp.sha(__file__),'stages':out,'checks':checks};sp.write(P/'surface-charged-flux.json',r);print(json.dumps({k:{a:b for a,b in v.items() if not a.startswith('packet_')} for k,v in out.items()},indent=2))
