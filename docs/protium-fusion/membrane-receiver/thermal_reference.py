#!/usr/bin/env python3
"""Bounded membrane heat/opacity screen; no hydro or material-lifetime model.
Run: python thermal_reference.py --root /path/to/website [--output /path/to/thermal.json]
Needs numpy. Inputs are the sibling pulse-loading.json and optional stopping-results.json.
NIST coefficients below are transcribed primary table rows, not fitted/extrapolated.
"""
from __future__ import annotations
import argparse, hashlib, json, math
from pathlib import Path
import numpy as np
SIGMA=5.670374419e-8
C=299792458.0
ENERGY=[1,2,3,5,10,30,100]
# NIST X-ray mass coefficients, (mu/rho, mu_en/rho), cm^2/g.
ELEMENTS={
 'C': [(2211,2209),(302.6,301.6),(90.33,89.63),(19.12,18.66),(2.373,2.078),(.2562,.06614),(.1514,.02147)],
 'O': [(4590,4576),(694.9,692.6),(217.1,215.8),(47.90,47.15),(5.952,5.565),(.3779,.1729),(.1551,.02355)],
 'Si':[(1570,1567),(2777,2669),(978.4,951.6),(245,240),(33.89,32.89),(1.436,1.164),(.1835,.04513)],
 'B': [(1229,1228),(159.7,159.1),(46.67,46.17),(9.683,9.332),(1.255,1.006),(.2063,.03506),(.1391,.01940)],
 'N': [(3311,3306),(476.9,475.5),(145.6,144.7),(31.44,30.86),(3.879,3.545),(.3066,.1099),(.1529,.02231)]}
AM={'C':12.011,'O':15.999,'Si':28.0855,'B':10.81,'N':14.007}
# cp=1000 is a controlled sensitivity, not a fitted temperature-dependent EOS.
# Carbon/BN cross-plane and in-plane k deliberately distinct; these are scenarios.
MATERIALS={
 'silica':{'rho':2200.,'cp':1000.,'k_lateral':1.3,'k_through':1.3,'parts':{'Si':1,'O':2}},
 'carbon':{'rho':2000.,'cp':1000.,'k_lateral':100.,'k_through':2.,'parts':{'C':1}},
 'BN':{'rho':2100.,'cp':1000.,'k_lateral':100.,'k_through':2.,'parts':{'B':1,'N':1}},
 'SiC':{'rho':3210.,'cp':1000.,'k_lateral':50.,'k_through':50.,'parts':{'Si':1,'C':1}}}
THICKNESS_NM=[10,30,100,300,1000,3000,10000]
URLS={e:f'https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z{z:02}.html' for e,z in [('C',6),('O',8),('Si',14),('B',5),('N',7)]}
URLS.update({'NIST_definitions':'https://physics.nist.gov/PhysRefData/XrayMassCoef/chap3.html',
 'silica_properties':'https://www.corning.com/media/worldwide/csm/documents/1e624b3034474193b5b00ea6f558dd3d2.pdf',
 'BN_room_temperature_microbridge':'https://pmc.ncbi.nlm.nih.gov/articles/PMC4855177/',
 'SiC_bulk_1073K':'https://doi.org/10.1016/j.jallcom.2018.05.321',
 'NASA_optical_cooled_wall':'https://ntrs.nasa.gov/citations/19710000905'})
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def coeff(m,e):
    p=MATERIALS[m]['parts']; total=sum(AM[s]*n for s,n in p.items()); j=ENERGY.index(e)
    return tuple(.1*sum(AM[s]*n/total*ELEMENTS[s][j][k] for s,n in p.items()) for k in [0,1])
def flash(m,dnm,e,t,Fcum,modes=128,substeps=4):
    """Insulated-face Fourier solution; exact piecewise-constant input, sampled maxima.
    Retains all mu_en deposition, ignoring radiation/phase change/electron escape.
    Returns sensible-heat-equivalent deltaT, not a material state above breakdown.
    """
    p=MATERIALS[m];rho=p['rho'];cp=p['cp'];d=dnm*1e-9;alpha=p['k_through']/(rho*cp)
    mu,me=coeff(m,e);a=mu*rho;b=me*rho;n=np.arange(modes+1,dtype=float);kn=n*np.pi/d
    an=2/d*b*a*(1-(-1)**n*np.exp(-a*d))/(a*a+kn*kn)
    an[0]=(b/a)*(-np.expm1(-a*d))/d
    lam=alpha*kn*kn; val=np.zeros_like(n); peak=0.;at=0.
    for lo,hi,flo,fhi in zip(t[:-1],t[1:],Fcum[:-1],Fcum[1:]):
        dt=(hi-lo)/substeps; q=(fhi-flo)/(hi-lo);dec=np.exp(-lam*dt)
        gain=np.empty_like(n);gain[0]=dt;gain[1:]=-np.expm1(-lam[1:]*dt)/lam[1:]
        for j in range(substeps):
            val=val*dec+q*an/(rho*cp)*gain
            if val.sum()>peak:peak=float(val.sum());at=float(lo+(j+1)*dt)
    av=float(val[0]); local=float(Fcum[-1]*me/cp)
    return {'surface_peak_sensible_delta_K':peak,'surface_peak_time_s':at,
            'final_average_sensible_delta_K':av,'local_no_conduction_delta_K':local,
            'through_diffusion_time_s':d*d/alpha,'modes':modes,'substeps_per_input_bin':substeps}
def pitch(k,d,q,delta=100.):return math.sqrt(8*k*d*delta/q) if q>0 else None
def grey_cavity(a,t):
    # Symmetric gray membrane, back side facing a black bath; r=1-a-t.
    # Inner J=(a*sigma*Tm^4+t*sigma*Tb^4)/(a+t).
    return a*(a+2*t)/(a+t) if a+t>0 else 0.
def main(root,out):
    here=root/'docs/protium-fusion/membrane-receiver'; pp=here/'pulse-loading.json'
    pl=json.loads(pp.read_text());s=pl['stages']['C13'];R=pl['radius_m'];area=4*math.pi*R*R
    prompt=s['prompt'];t=np.array(prompt['radial_ray_arrival_time_s']);t=t-t[0]
    F=np.array(prompt['cumulative_fluence_J_m2']);peakph=float(np.max(np.diff(F)/np.diff(t)))
    peakion=s['smooth']['peak_kinetic_power_W_m2']
    old=json.loads((root/'docs/protium-fusion/open-loop-redesign/collector-inputs.json').read_text())['stages']['c13']['record']
    K=old['final_kinetic_J']/area;U=old['final_internal_J']/area
    rows=[];checks=[]
    for mat,p in MATERIALS.items():
        for dnm in THICKNESS_NM:
            d=dnm*1e-9;arealm=p['rho']*d
            for e in ENERGY:
                mu,me=coeff(mat,e);inter=-math.expm1(-mu*arealm);dep=me/mu*inter
                z=flash(mat,dnm,e,t,F)
                checks.append(z['surface_peak_sensible_delta_K']+1e-6>=z['final_average_sensible_delta_K'])
                checks.append(z['surface_peak_sensible_delta_K']<=z['local_no_conduction_delta_K']*(1+1e-9))
                avg=F[-1]*dep/(arealm*p['cp'])
                checks.append(abs(z['final_average_sensible_delta_K']/avg-1)<1e-12)
                rows.append({'material':mat,'thickness_nm':dnm,'photon_keV':e,'areal_mass_kg_m2':arealm,
                    'whole_shell_film_mass_kg':arealm*area,'mu_m2_kg':mu,'mu_en_m2_kg':me,
                    'uncollided_transmission':math.exp(-mu*arealm),'local_retained_energy_fraction_proxy':dep,
                    'retained_fluence_J_m2_at_eta1':float(F[-1]*dep),
                    'retained_peak_flux_W_m2_at_eta1':peakph*dep,
                    'steady_cooling_pitch_for_100K_peak_at_eta1_m':pitch(p['k_lateral'],d,peakph*dep),
                    **z})
    refinement=[]
    for mat,d,e in [('silica',10000,1),('silica',10000,3),('carbon',10000,1),('BN',1000,3),('SiC',10000,3)]:
        lo=flash(mat,d,e,t,F);hi=flash(mat,d,e,t,F,512,16)
        refinement.append({'material':mat,'thickness_nm':d,'photon_keV':e,'relative_peak_change':hi['surface_peak_sensible_delta_K']/lo['surface_peak_sensible_delta_K']-1})
        checks.append(abs(refinement[-1]['relative_peak_change'])<.003)
    cool=[]
    for mat,p in MATERIALS.items():
        for dnm in THICKNESS_NM:
            for fi in [0.001,.01,.1,1.]:
                q=fi*peakion;L=pitch(p['k_lateral'],dnm*1e-9,q)
                cool.append({'material':mat,'thickness_nm':dnm,'net_heat_fraction_of_kinetic':fi,
                  'conditional_peak_heat_W_m2':q,'cooling_pitch_for_100K_m':L,
                  'lateral_relaxation_time_s':L*L/(np.pi**2*p['k_lateral']/(p['rho']*p['cp'])),
                  'buffer_kg_m2_at_cp1000_delta100K':fi*K/1e5,
                  'back_bath1000K_equilibrium_K_at_effective_emissivity1':(1000.**4+q/SIGMA)**.25})
    radiation=[]
    for Tm in [1000,1200,1500,2000]:
        for Tb in [450,1000,1500,2000,3000,5000]:
            for eff in [.001,.01,.1,1.]:
                radiation.append({'membrane_K':Tm,'back_bath_K':Tb,'effective_emissivity_viewfactor':eff,
                    'outward_cooling_W_m2':eff*SIGMA*(Tm**4-Tb**4)})
    ion=[];stop=here/'stopping-results.json'
    if stop.exists():
        st=json.loads(stop.read_text())
        for r in st['direct_rows']:
            if r['material'] not in MATERIALS:continue
            p=MATERIALS[r['material']];d=r['thickness_nm']*1e-9;frac=r['energy_loss_fraction'];q=frac*peakion
            ion.append({'material':r['material'],'thickness_nm':r['thickness_nm'],
                'projectile_removal_fraction':frac,'projectile_removal_J_m2':r['projectile_energy_loss_J_m2'],
                'no_cooling_sensible_equivalent_delta_K':r['projectile_energy_loss_J_m2']/(p['rho']*d*p['cp']),
                'kinetic_waveform_scaled_peak_W_m2_DIAGNOSTIC':q,
                'cooling_pitch100K_from_scaled_waveform_m_DIAGNOSTIC':pitch(p['k_lateral'],d,q),
                'buffer_kg_m2_at_cp1000_delta100K':r['projectile_energy_loss_J_m2']/1e5})
    mixture={m:{'mass_fraction':{s:AM[s]*n/sum(AM[j]*nj for j,nj in p['parts'].items()) for s,n in p['parts'].items()},
      'rows':[{'photon_keV':e,'mu_m2_kg':coeff(m,e)[0],'mu_en_m2_kg_approximate':coeff(m,e)[1],
       'thin_limit_sensible_delta_K_eta1':float(F[-1]*coeff(m,e)[1]/p['cp'])} for e in ENERGY]} for m,p in MATERIALS.items()}
    cavity=[]
    for a,tr in [(1,0),(.9,0),(.1,0),(.01,0),(.5,.5),(.1,.9),(.01,.99)]:
        k=grey_cavity(a,tr);cavity.append({'absorptivity_emissivity':a,'transmission':tr,
             'reflectivity':1-a-tr,'whole_uniform_shell_effective_outward_coefficient':k})
        checks.append(0<=k<=1)
    data={'scope':'Parameterized heat/opacity screen for a uniform supported membrane facing a renewable absorber. Not material lifetime qualification.',
      'source_sha256':{'pulse-loading.json':sha(pp),'stopping-results.json':sha(stop) if stop.exists() else None},
      'code_sha256':sha(__file__),'source_urls':URLS,
      'pulse':{'stage':'C13','R_m':R,'area_m2':area,'prompt_fluence_J_m2':float(F[-1]),
       'prompt_radial_ray_duration_s':float(t[-1]),'prompt_piecewise_peak_W_m2':peakph,
       'prompt_emission_duration_s':float(prompt['emission_time_s'][-1]),
       'kinetic_fluence_J_m2':K,'residual_internal_fluence_J_m2_separate':U,
       'kinetic_smooth_peak_W_m2':peakion,'arrival_1_50_99_percent_energy_s':s['smooth']['arrival_1_50_99_percent_energy_s'],
       'electron_bookkeeping':s['electron_bookkeeping']},
      'assumptions':{'materials':MATERIALS,'cp_scope':'Constant 1000J/kg/K sensitivity. Large deltaT values are sensible-heat equivalents; phase changes and decomposition occur first.',
       'k_scope':'Silica k1.3 is room-temperature bulk data. Carbon and BN lateral100/through2 and SiC50 are controlled material scenarios, not hot irradiated film qualifications.',
       'photon_scope':'NIST bulk coefficients. Mass-fraction mu additive; approximate mu_en additive. Fluence partition neglects secondary transport. No extrapolation below1keV.',
       'electron_retention':'Rows show eta_retained=1. Linear scale all photon temperature and energy results by eta=.1,.01 as sensitivities. Secondary electrons arriving from the absorber are extra; electrostatic self-charging can suppress net escape.',
       'flash':'Insulated faces, exponential primary deposition, Fourier conduction, fixed properties. No phase change, lateral supports or radiation. Actual radial-ray pulse bins; finite-angle broadening unresolved.',
       'cooling':'Conductive sink lines assumed isothermal. Lateral L uses strip solution qL²/(8kd); ignores contact resistance, support heat and coolant limitations.',
       'ion_scope':'Generic fractions are sensitivities. Historical scaled-peak diagnostics are preserved here. Exact timed stopping peaks are now available and used by mechanics.py. Projectile energy removal requires an independent local heat-retention model.',
       'cavity':'Uniform closed sphere; no steady radiative heat removal to central vacuum. Opaque inward radiation is returned. Semitransparent gray examples face a black exterior bath.'},
      'mixture_coefficients':mixture,'photon_flash_rows':rows,'modal_refinement':refinement,
      'cooling_fraction_rows':cool,'back_bath_exchange_rows':radiation,'coupled_stopping_diagnostics':ion,
      'cavity':{'max_return_delay_s':2*R/C,'mean_return_delay_s':4*R/(3*C),
       'Lambertian_return_pdf':'K(t)=2*t/tmax² for0<t<tmax; integral1, mean2tmax/3.',
       'closed_opaque_uniform_steady_front_cooling_W_m2':0.,
       'symmetric_gray_film_black_back_bath_formula':'q_net=a*(a+2*t)/(a+t)*sigma*(Tm^4-Tb^4), r=1-a-t.',
       'gray_rows':cavity},
      'flash_absolute_radiation_upper_J_m2':{str(T):2*SIGMA*T**4*t[-1] for T in [1000,1500,2000]},
      'radiation_upper_scope':'Two black faces during the microsecond flash, no incoming photons; only an absolute local upper bound. Whole-shell steady rejection is limited by the actual outward sink.',
      'checks':{'passed':all(checks),'count':len(checks)}}
    out.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'checks':data['checks'],'pulse':data['pulse'],'refinement':refinement,'ion_rows':len(ion)},indent=2))
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[3]);a.add_argument('--output',type=Path)
    x=a.parse_args();main(x.root,x.output or Path(__file__).with_name('thermal-reference.json'))
