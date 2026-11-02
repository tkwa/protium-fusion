"""Conservative residence-time screens for omitted two-proton captures.

Use the next hotter tabulated upper rate throughout a pulse, with its maximum
proton abundance, initial density and total modeled duration. These are upper
estimates within the published rate envelope, not rigorous nuclear bounds.
Reverse photodisintegration and reactant depletion would reduce leakage.
"""
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).parent
T=np.array([.3,.5,.6,.8,1.,1.5,2.,3.,5.])*1e9
# Grigorenko & Zhukov (2005), Tables I and IV. NA^2 <sigma v>,
# conventional cm^6 mol^-2 s^-1. Include the identical-proton factor 1/2.
RATES={
 'o15':np.array([2.9e-14,1.3e-12,2.8e-12,6.9e-12,1.3e-11,3.8e-11,6.8e-11,1.1e-10,2.3e-10]),
 'ne18':np.array([2.5e-25,1.8e-18,8.3e-17,9.5e-15,1.8e-13,1.1e-11,8e-11,4.4e-10,1.2e-9])}
cases=[('o15','radial-A-300m-1024-reduced.json',1/15),
       ('ne18','pulse-C0-5e+08-1e+10-30-fast.json',None)]
rows=[]
for isotope,file,yp in cases:
    d=json.loads((ROOT/file).read_text())
    temp=d.get('peak_temperature_K',d.get('peak_T_K'))
    j=int(np.searchsorted(T,temp));rate=float(np.max(RATES[isotope][:j+1]))
    if yp is None: yp=d['peak_mass_fractions']['p']
    rho=d['rho0_kg_m3']/1000
    exposure=.5*rho**2*yp**2*rate*d['elapsed_s']
    rows.append({'source_result':file,'reactant':isotope,
      'tabulated_temperature_ceiling_K':float(T[j]),
      'upper_rate_cm6_mol2_s':rate,'proton_molar_abundance_ceiling':yp,
      'constant_density_assumption_g_cm3':rho,'integrated_exposure_ceiling':exposure,
      'capture_probability_estimate_ceiling':float(-np.expm1(-exposure))})
out={'source':'https://doi.org/10.1103/PhysRevC.72.015803',
 'formula':'lambda = 0.5 * rho_cgs^2 * Y_p^2 * NA^2<sigma v>',
 'scope':'Published resonant-rate envelope and bounding pulse conditions; screening and newer rate evaluation not included.',
 'rows':rows}
(ROOT/'two-proton-screen.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
