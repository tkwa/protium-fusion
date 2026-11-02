"""Conditional helium first-passage, retention and extraction requirements.
Run with numpy/scipy: python helium_release.py
The Jung coefficient is an abstract-reported mobile/orientation reference;
no trapped fraction or effective diffusivity of the selected graphite is fitted.
"""
from pathlib import Path
import hashlib,json,math
import numpy as np
from scipy.special import erfc,erfcinv,spence
from scipy.optimize import brentq
P=Path(__file__).resolve().parent
KB=1.380649e-23;KBEV=8.617333262145e-5;MU=1.66053906660e-27
SIGMA=5.670374419e-8;YEAR=365.25*86400
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()

def diffusivity(T,material='graphite',policy='extrapolate'):
    T=np.asarray(T,dtype=float)
    if policy=='hold_1050':T=np.minimum(T,1050.)
    D0,E={'graphite':(5e-5,1.17),'SiC':(1.1e-6,1.14)}[material]
    return D0*np.exp(-E/(KBEV*T))

def released(depth,integrated_D):
    x,A=np.broadcast_arrays(np.asarray(depth,dtype=float),np.asarray(integrated_D,dtype=float))
    out=np.zeros_like(x);ok=A>0
    out[ok]=erfc(x[ok]/(2*np.sqrt(A[ok])))
    out[x<=0]=1.
    return out

def t_fraction(depth,D,f):
    return np.asarray(depth,dtype=float)**2/(4*D*erfcinv(f)**2)

def strip_front_probability(depth,pitch,eta):
    # Uniform lateral implantation in a strip with absorbing front and side walls.
    # Sum of odd n terms written via Li_2(z)=spence(1-z).
    z=np.exp(-np.pi*np.asarray(depth,dtype=float)*math.sqrt(eta)/pitch)
    return (8/np.pi**2)*(spence(1-z)-.25*spence(1-z*z))

def weighted_quantile(x,w,q):
    k=np.argsort(x);c=np.cumsum(w[k]);return float(x[k[min(np.searchsorted(c,q*c[-1]),len(k)-1)]])

def main():
    sf=P.parent/'stopping/surface-spectra.json';s=json.loads(sf.read_text())
    b=s['stages']['B1'];period=b['period_s'];he=b['species']['he4'];fluence=he['nuclei_m2']
    rows=[]
    for T in [750,1000,1050,1500,1668,1700]:
      for material in ['graphite','SiC']:
       for x in np.array([.1,1,10,100])*1e-6:
        D=float(diffusivity(T,material));Dc=float(diffusivity(T,material,'hold_1050'))
        rows.append(dict(T_K=T,material=material,depth_m=float(x),D_m2_s=D,
          coefficient_within_reported_T_domain=(750<=T<=1050 if material=='graphite' else 800<=T<=1050),
          median_first_passage_s=float(t_fraction(x,D,.5)),release95_time_s=float(t_fraction(x,D,.95)),
          mobile_released_after_B1_period=float(released(x,D*period)),
          mobile_released_after_B1_period_hold1050=float(released(x,Dc*period)),
          max_D_for_95pct_retention_m2_s=float(x*x/(4*period*erfcinv(.05)**2)),
          max_D_for_99pct_retention_m2_s=float(x*x/(4*period*erfcinv(.01)**2)),
          eventual_back_fraction_for_isotropic_175mm_slab=float(x/.175)))
    retention=[]
    for T in [1000,1500,1668,1700]:
      for hold in [100.,period,YEAR]:
       for nu in [1e12,1e13]:
        retention.append(dict(T_K=T,hold_s=hold,attempt_frequency_s=nu,
          E_for_95pct_survival_eV=KBEV*T*math.log(nu*hold/(-math.log(.95))),
          scope='single first-order trap, no retrapping/irradiation detrapping; requirement, not measured He barrier'))
    geometry=[]
    for p in [20e-6,100e-6]:
     for x in [.1e-6,1e-6,10e-6]:
      vals={str(eta):float(strip_front_probability(x,p,eta)) for eta in [1,10,100,1000,10000]}
      req={str(eps):float(10**brentq(lambda lg:float(strip_front_probability(x,p,10**lg))-eps,-8,16)) for eps in [.05,.01,.001]}
      geometry.append(dict(pitch_m=p,depth_m=x,front_probability_by_anisotropy=vals,
        required_Dparallel_over_Dnormal_for_front_fraction=req,
        mean_side_exit_time_ignoring_front_s=p*p/(12*float(diffusivity(1000)))))
    # Equivalent room-temperature throughput is a convention, not gas temperature in the panel.
    qshot=fluence*KB*273.15
    cryo=dict(He_B1_mass_kg_m2=he['mass_kg_m2_baryon'],He_mean_kg_m2_s=he['mass_kg_m2_baryon']/period,
      B1_equivalent_gas_load_273K_Pa_m3_per_receiver_m2=qshot,
      cryopanel_area_fraction_cases=[])
    cp=2.5*KB/(4*MU)
    for f in [1e-6,1e-4,.01]:
      cryo['cryopanel_area_fraction_cases'].append(dict(area_fraction=f,
       load_per_B1_shot_Pa_m3_per_panel_m2=qshot/f,
       shots_to_1000_Pa_m3_per_panel_m2=1000*f/qshot,
       black_200K_view_heat_W_per_receiver_m2=f*SIGMA*(200**4-5**4),
       reversible_electric_for_that_heat_to_200K_W_per_receiver_m2=f*SIGMA*(200**4-5**4)*(200/5-1),
       black_1000K_view_heat_W_per_receiver_m2=f*SIGMA*(1000**4-5**4)))
    cryo.update(dict(ideal_constant_pressure_He_sensible_cooling_200_to_5K_J_kg=cp*(200-5),
      reversible_work_for_that_cooling_with_200K_sink_J_kg=cp*(200*math.log(200/5)-(200-5)),
      scope='Sensible-gas thermodynamic bounds only. Adsorption heat, compression, panel regeneration heat, finite efficiency and installed mass excluded. 1000Pa m3/m2 is illustrative graph-scale working load, not universal capacity or guaranteed chamber pressure.'))
    checks=[]
    def check(name,ok):checks.append(dict(name=name,passed=bool(ok)));assert ok,name
    check('median identity',abs(float(released(1e-6,float(diffusivity(1000))*t_fraction(1e-6,float(diffusivity(1000)),.5)))-.5)<1e-14)
    check('zero front depth',abs(float(strip_front_probability(0,1.,1.))-1)<1e-14)
    check('strip series',abs(float(strip_front_probability(.1,1.,100.))-(8/math.pi**2)*sum(math.exp(-n*math.pi)/(n*n) for n in range(1,1000,2)))<1e-14)
    for row in rows:
      check('range '+str((row['T_K'],row['material'],row['depth_m'])),0<=row['mobile_released_after_B1_period']<=1)
    for row in geometry:
      vv=list(row['front_probability_by_anisotropy'].values());check('monotone anisotropy '+str((row['pitch_m'],row['depth_m'])),all(a>=b for a,b in zip(vv,vv[1:])))
    result=dict(scope='Conditional material/geometry kernels; no measured reactor retention fraction or selected anisotropy. Temperatures above1050K explicitly extrapolated.',
      source_spectra=str(sf),source_spectra_sha256=sha(sf),code_sha256=sha(__file__),
      B1_period_s=period,B1_He_nuclei_m2=fluence,B1_He_mass_ug_m2=he['mass_ug_m2_baryon'],
      diffusion_rows=rows,trap_hold_requirements=retention,lamellar_geometry=geometry,cryo_bounds=cryo,checks=checks)
    (P/'helium-release.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(checks=len(checks),all_pass=all(x['passed'] for x in checks),D1000=float(diffusivity(1000)),B1period=period,geometry=geometry[1],cryo=cryo),indent=2))
if __name__=='__main__':main()
