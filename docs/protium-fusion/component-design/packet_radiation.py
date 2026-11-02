#!/usr/bin/env python3
"""Incoming-packet radiation acceptance screen. NumPy only; no network.
python packet_radiation.py --output packet-radiation.json
Optically thick, smooth unoxidized Al reference uses Rakic (1995) Table 2.
Room-temperature reference is not a qualified cryogenic grain coating.
"""
import argparse,hashlib,json,math
from pathlib import Path
import numpy as np
SIGMA=5.670374419e-8; KB_EV=8.617333262e-5
# E(eV), n, k. Primary table https://doi.org/10.1364/AO.34.004755
AL=np.array([list(map(float,r.split())) for r in '''
.0061993 423.96 483.70
.007 397.93 458.50
.008059 364.04 429.62
.009 339.62 408.92
.0099188 318.81 391.71
.012399 274.38 354.35
.015498 233.56 321.08
.018598 202.63 295.42
.021697 177.93 275.34
.024797 157.30 258.26
.027897 140.05 243.43
.030996 125.14 230.19
.037196 102.10 208.10
.04 91.955 199.99
.045 75.748 181.78
.05 68.535 164.81
.055 63.554 153.45
.06 58.580 144.23
.065 54.413 136.09
.07 50.951 129.49
.08 43.775 118.39
.09 38.461 108.96
.1 33.519 101.28
.12 26.216 88.197
.14 20.837 78.274
.16 16.755 69.857
.18 14.088 62.841
.2 12.195 57.156
.22 10.742 52.518
.24 9.5580 48.593
.26 8.5881 45.257
.28 7.7757 42.367
.3 7.0796 39.826
.32 6.4808 37.595
.34 5.9564 35.608
.36 5.4903 33.814
.38 5.0735 32.183
.4 4.7097 30.737
.45 3.9380 27.580
.5 3.3372 25.004
.6 2.4738 20.982
.7 1.9205 17.991
.8 1.5782 15.656
.9 1.3899 13.784
1 1.3157 12.245
1.1 1.3281 10.969
1.2 1.3998 9.8914
1.24 1.4359 9.4939
1.28 1.4867 9.0655
1.32 1.6784 8.5970
1.36 1.9739 8.3058
1.4 2.2802 8.1134
1.48 2.6945 8.1878
1.52 2.7668 8.2573
1.56 2.7675 8.3866
1.6 2.6154 8.4914
1.7 2.1606 8.3565
1.8 1.8301 8.0601
1.9 1.5724 7.7354
2 1.3660 7.4052
2.2 1.0728 6.7839
2.4 .87340 6.2418
2.6 .72780 5.7781
2.8 .60790 5.3676
3 .52135 5.0008
3.4 .39877 4.3957
3.8 .31474 3.9165
4 .28003 3.7081
5 .18137 2.9029
6 .12677 2.3563
'''.strip().splitlines()])
# Chang (1974), Table 2; molar convention CH2=14.027g. H(T)-H(0) J/mol.
POLY_T=np.array([20,50,100,150,200,250,300.])
POLY_H={'linear':np.array([4.85,83.55,448.9,1016,1746,2653,3778]),'branched':np.array([7.32,95.63,468.7,1045,1808,2825,4241])}
def planck(e,T):
 x=e/(KB_EV*T)
 return e**3/np.expm1(np.minimum(x,700))/(math.pi**4/15*(KB_EV*T)**4)
def optics(e,angular=160):
 n=np.exp(np.interp(np.log(e),np.log(AL[:,0]),np.log(AL[:,1])))+1j*np.exp(np.interp(np.log(e),np.log(AL[:,0]),np.log(AL[:,2])))
 x,w=np.polynomial.legendre.leggauss(angular);mu=(x+1)/2;w=w/2
 ct=np.sqrt(n[:,None]**2-1+mu[None,:]**2)
 rs=(mu[None,:]-ct)/(mu[None,:]+ct)
 rp=(n[:,None]**2*mu[None,:]-ct)/(n[:,None]**2*mu[None,:]+ct)
 absorption=1-(abs(rs)**2+abs(rp)**2)/2
 return np.sum(absorption*2*mu*w,axis=1), 1-abs((1-n)/(1+n))**2,n

def calc(points=8000):
 e=np.geomspace(AL[0,0],AL[-1,0],points);a,an,n=optics(e)
 baseT=1000.;r=1.1e9;v=2.4e6;rho=1000.;flight=r/v
 def alpha(T):return float(np.trapz(planck(e,T)*a,e))
 hot_alpha=alpha(baseT);normal=float(np.trapz(planck(e,baseT)*an,e))
 low=np.geomspace(1e-8,e[0],1600);high=np.geomspace(e[-1],100,800)
 tail=float(np.trapz(planck(low,baseT),low)+np.trapz(planck(high,baseT),high))
 allowances={k:float((h[-1]-h[0])/.014027) for k,h in POLY_H.items()}
 rows=[];checks=[]
 for radius in [1e-4,1e-3,1e-2]:
  dose=3*SIGMA*baseT**4*flight/(rho*radius)
  row={'radius_m':radius,'black_absorber_dose_J_kg':dose,'optically_thick_Al_reference_dose_before_self_emission_J_kg':dose*hot_alpha,
       'allowed_Planck_absorptivity_polyethylene20to300K':{k:h/dose for k,h in allowances.items()},
       'allowed_absorptivity_for_explicit_30kJkg_budget':30000/dose,
       'Al100nm_coating_mass_per_core_mass':3*2700*1e-7/(rho*radius),
       'polyethylene_reference_transit_limits':{}}
  # Integrate actual enthalpy increments with temperature-dependent emission using Kirchhoff.
  for kind,H in POLY_H.items():
   tt=np.linspace(20,300,4001);hh=np.interp(tt,POLY_T,H)/.014027
   emiss=np.array([alpha(T) for T in np.linspace(20,300,100)])
   ep=np.interp(tt,np.linspace(20,300,100),emiss)
   rate=3*SIGMA/(rho*radius)*(hot_alpha*baseT**4-ep*tt**4)
   time=float(np.sum(np.diff(hh)/(.5*(rate[:-1]+rate[1:]))))
   row['polyethylene_reference_transit_limits'][kind]={'time_to300K_s':time,'maximum_unshielded_path_m':time*v,'full_transit_at_or_below300K':time>=flight}
  rows.append(row);checks.append(all(0<x<1 for x in row['allowed_Planck_absorptivity_polyethylene20to300K'].values()))
 # Maximum normal single-pass material leakage through 100nm unoxidized bulk Al.
 lam=1.239841984e-6/e
 attenuation=np.exp(-4*math.pi*n.imag*1e-7/lam)
 # This omits incident reflection and is a loose straight normal-path material leakage screen.
 leak=float(np.trapz(planck(e,baseT)*attenuation,e))
 checks.extend([0<normal<hot_alpha<1,abs(float(np.trapz(planck(e,baseT),e))+tail-1)<2e-6,leak<1e-3])
 return {'scope':'Geometric-optics sphere inside uniform 1000K cavity. Radiation only; room-temperature smooth unoxidized thick Al optical reference, no coating qualification.','parameters':{'radius_receiver_m':r,'speed_m_s':v,'flight_s':flight,'core_density_kg_m3':rho,'wall_temperature_K':baseT},'primary_sources':{'Al_table':'https://doi.org/10.1364/AO.34.004755','Al_public_table':'https://opg.optica.org/ao/abstract.cfm?uri=ao-34-22-4755','polyethylene':'https://pmc.ncbi.nlm.nih.gov/articles/PMC6565414/'},'optics':{'Al_reference_hemispherical_Planck_absorptivity':hot_alpha,'Al_reference_normal_Planck_absorptivity':normal,'unresolved_Planck_energy_fraction_zero_full_bound':tail,'Al100nm_normal_material_transmission_upper_resolved':leak,'table_E_n_k':AL.tolist()},'polyethylene_enthalpy20to300_J_kg':allowances,'cases':rows,'emission_assumption':'Same temperature-independent optical spectrum used for Kirchhoff emissivity. This is a conditional room-temperature optical reference; no cryogenic improvement granted.','unknowns':['Cryogenic film absorption, oxidation, roughness and size-scattering','Coating mechanical/electrical integrity and contamination by high-Z addition','Actual cryogenic nitrogen compound enthalpy, gas retention and phase stability','Nonuniform internal grain temperatures and thermal stress; energy acceptance is necessary but not sufficient'],'checks':{'passed':sum(checks),'count':len(checks),'all_passed':all(checks)},'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',default=str(Path(__file__).with_name('packet-radiation.json')));a=p.parse_args();o=calc();Path(a.output).write_text(json.dumps(o,indent=2)+'\n');print(json.dumps({k:o[k] for k in ['optics','polyethylene_enthalpy20to300_J_kg','cases','checks']},indent=2))
