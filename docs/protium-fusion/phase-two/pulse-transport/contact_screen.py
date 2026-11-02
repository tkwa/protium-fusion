"""Flag cold-contact states where an ideal fully ionized isentrope fails.

No ionization fractions or corrected isentrope are claimed. Hydrogenic
isolated-atom binding energies give energy scales, not dense-plasma thresholds.
"""
import argparse,json,hashlib,sys
from pathlib import Path
import numpy as np
from microphysics import *
ap=argparse.ArgumentParser()
ap.add_argument('--repo',default=str(Path(__file__).resolve().parents[2]))
ap.add_argument('--out',default=str(Path(__file__).parent))
a=ap.parse_args();ROOT=Path(a.repo);OUT=Path(a.out);OUT.mkdir(parents=True,exist_ok=True)
sys.path.insert(0,str(ROOT))
import radial_pulse as rp
import pulse_model as pm
p=ROOT/'phase-two/driver/delivery-screen.json'
data=json.loads(p.read_text());rows=[]
rydberg=13.605693123 # eV, infinite-nuclear-mass hydrogenic scale
for case in data['cases']:
    Y,_,_=rp.initial(case['stage'])
    assert np.all(pm.ZZ[Y>0]<=10)
    ionbinding=np.array([rydberg*z*z*(min(z,2)+max(z-2,0)/4) for z in pm.ZZ])
    # All species present are H/He/C/N. Higher shells are irrelevant here.
    upper=E/MU*np.sum(Y*ionbinding)
    dom=int(np.argmax(Y*pm.ZZ));zd=pm.ZZ[dom]
    for state in case['ideal_fully_ionized_contact_isentropes']:
        rho=state['density_kg_m3'];T=state['temperature_K']
        ne=rho/MU*np.sum(Y*pm.ZZ)
        ui,pi,gd,ge=coulomb(np.asarray(T),np.asarray(rho),Y,pm.ZZ)
        eta,nm,np_=fermions(T,ne)
        ae=(3/(4*np.pi*ne))**(1/3);ai=zd**(1/3)*ae
        # Hydrogenic last-electron binding and orbital radius are warning
        # scales only; screening and pressure ionization change both.
        shellbinding=rydberg*zd*zd
        rows.append({'stage':case['stage'],'zones':case['zones'],'tag':case['tag'],'rho_kg_m3':rho,'T_K':T,
           'driver_table_floor_reached':not state['reached_requested_density'],
           'ideal_specific_internal_J_kg':state['specific_internal_J_kg'],
           'fully_stripped_dominant_ion_Gamma':float(gd),
           'fully_stripped_ion_ion_pressure_fraction':float(-pi/state['pressure_Pa']),
           'fully_stripped_electron_eta':float(eta),
           'hydrogenic_K_binding_over_kT':shellbinding*E/(KB*T),
           'hydrogenic_K_orbit_over_ion_sphere':(5.29177210544e-11/zd)/ai,
           'isolated_atom_total_ionization_upper_J_kg':float(upper),
           'ionization_upper_over_ideal_contact_internal':float(upper/state['specific_internal_J_kg']),
           'ionization_upper_over_whole_prepared_input':float(upper*case['mass_kg']/case['input_energy_J'])})
out={'scope':'Flags and energy scales only: a chemically consistent nonideal EOS is required for contact entropy, pressure and sound speed.',
 'input_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),
 'hydrogenic_bound':'Nonrelativistic isolated atoms, shell occupancy n=1 then n=2 with electron-electron repulsion omitted. This overestimates binding and is not an actual ionization cost.',
 'rows':rows}
(OUT/'contact-regime.json').write_text(json.dumps(out,indent=2)+'\n')
for r in rows:
 if abs(r['rho_kg_m3']-2000)<1:
  print(r)
