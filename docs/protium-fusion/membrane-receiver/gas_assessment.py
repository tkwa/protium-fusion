#!/usr/bin/env python3
"""Independent neutral-column and incoming-packet screen from saved plume outputs.
No source/hydrodynamic rerun. Run python gas_assessment.py --root /path/to/website.
Dependencies numpy/scipy; output gas-assessment.json next to this script by default.
Subjective intervals below are explicit conditional judgments, not measured posteriors.
"""
from pathlib import Path
import argparse,hashlib,json,math
import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
from scipy.special import ndtri,expit
K=1.380649e-23;MU=1.66053906660e-27;M_C=12.0107*MU;YEAR=365.25*86400

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def mono_column(a):
    # Radius-normalized radial column for mono-speed cosine wall emission, first hit sticks.
    if a<=0:return 1.
    if a>=2:return 0.
    if a<1e-3:return 1+sum(a**k*(1/k+1/(k+2)) for k in range(1,9))
    return -1/a+.5-(1+a*a)/(a*a)*math.log(abs(1-a))
def column(w):
    # w=t*sqrt(2kT/m)/R; x=mv²/(2kT), flux-Maxwell density x exp(-x).
    top=min(65.,(2/w)**2);sing=1/w**2
    return quad(lambda x:x*math.exp(-x)*mono_column(w*math.sqrt(x)),0,top,
        points=[sing] if sing<top else None,epsabs=2e-10,limit=200)[0]
def column_phase_bound(delta,vt,R,gmax=3.26):
    # Gamma(w)<=gmax checked numerically; exact large-time inequality Gamma<=8/w^4.
    # Bound every older shot independently, then integral-bound the residual series.
    n=2000;k=np.arange(1,n+1,dtype=float);fac=delta*vt/R
    return gmax+float(np.minimum(gmax,8/(fac*k)**4).sum())+8/fac**4/(3*n**3)
def quant(a):return {str(q):float(np.quantile(a,q)) for q in [.1,.5,.9]}
def main(root,out):
    here=root/'docs/protium-fusion/membrane-receiver';p=here/'results.json';s=json.loads(p.read_text())
    R=s['radius_m'];area=s['area_m2'];V=area*R/3;T=1700.
    mass=np.arange(1,4)*M_C;vbar=np.sqrt(8*K*T/(math.pi*mass));vt=np.sqrt(2*K*T/mass)
    q=np.zeros(3);phase=np.zeros(3);shot=[]
    for name,r in s['pulse_cases'].items():
        sigma=np.array(r['total_mass_kg_m2_by_Cn']);dt=r['mean_interval_s'];q+=sigma*area/dt
        factors=np.array([column_phase_bound(dt,v,R) for v in vt]);phase+=sigma*factors
        shot.append({'stage':name,'single_shot_mass_kg':r['total_mass_kg'],
            'single_shot_atom_column_max_m2_at_factor326':float(3.26*sigma.sum()/M_C),
            'independent_phase_bound_mass_column_kg_m2':float((sigma*factors).sum()),
            'phase_bound_factors_C1_C2_C3':factors.tolist()})
    tau=4*R/(3*vbar);inventory=q*tau;column_mass=3*inventory/area
    Ngas=float(np.sum(column_mass/mass));Natom=float(np.sum(column_mass)/M_C)
    density=float(np.sum(inventory/mass)/V);pressure=density*K*T
    # A deliberately conservative species-level encounter screen counts C atoms independently.
    base={'gas_T_K':T,'source_kg_s':float(q.sum()),'source_C1_C2_C3_kg_s':q.tolist(),
      'unit_sticking_inventory_kg':float(inventory.sum()),'mean_residence_by_species_s':tau.tolist(),
      'mass_density_kg_m3':float(inventory.sum()/V),'molecular_number_density_m3':density,
      'equivalent_kinetic_pressure_Pa':pressure,'mean_radial_mass_column_kg_m2':float(column_mass.sum()),
      'mean_radial_molecular_column_m2':Ngas,'mean_radial_carbon_atom_column_m2':Natom,
      'phase_independent_mass_column_bound_kg_m2':float(phase.sum()),
      'phase_independent_carbon_atom_column_bound_m2':float(phase.sum()/M_C)}
    gfit=minimize_scalar(lambda l:-column(math.exp(l)),bounds=(-4,2),method='bounded')
    grid=np.geomspace(.002,20,240); vals=np.array([column(w) for w in grid]);gmax=-gfit.fun
    ct=quad(lambda a:2*a**3*mono_column(a),0,2,points=[1.],epsabs=1e-10)[0]
    checks=[abs(q.sum()/s['average_gross_evaporation_kg_s']-1)<1e-12,gmax<3.26,vals.max()<3.26,abs(ct-8)<1e-10]
    scenarios=[]
    # m/s intentionally variable: these are assessment points, not selected launcher designs.
    packet=[];dust=[]
    for G,stick,label in [(1,1,'saved_unit_sticking'),(1,.5,'subjective_central_inputs'),(400,1,'400x_surface_stress'),(400,.1,'400x_and_low_capture_stress')]:
        cm=base['phase_independent_mass_column_bound_kg_m2']*G/stick
        na=cm/M_C
        scenarios.append({'label':label,'source_multiplier':G,'effective_condensing_probability':stick,
          'thermal_source_kg_s':q.sum()*G,'phase_mass_column_kg_m2':cm,'phase_atom_column_m2':na,
          'encounter_probability_at_sigma1e_minus19':-math.expm1(-na*1e-19),
          'encounter_probability_at_sigma1e_minus18':-math.expm1(-na*1e-18),
          'sigma_max_for_encounter_fraction1e_minus3_m2':-math.log1p(-.001)/na,
          'fraction_structural_C_allowance_if_all_escapes':q.sum()*G/s['material_C_allowance_kg_s']})
        for rho in [100.,1000.,2000.]:
            for speed in [1e6,1e7,2e7]:
                for radius in [1e-6,1e-5,1e-4,1e-3,1e-2]:
                    mp_area=4*rho*radius/3;heat=.5*cm*speed**2/mp_area
                    packet.append({'gas_case':label,'packet_radius_m':radius,'packet_density_kg_m3':rho,
                      'speed_m_s':speed,'all_intercepted_relative_kinetic_energy_retained_J_kg':heat,
                      'intercepted_mass_fraction':cm/mp_area,
                      'radius_min_for_heat1e4_J_kg_m':3*cm*speed**2/(8*rho*1e4)})
    for ap in [1e-4,1e-3,1e-2]:
        rho_p=1000.;rho_d=2000.;speed=1e7;eta=.5;budget=1e4;mp=4*math.pi/3*rho_p*ap**3
        md=2*budget*mp/(eta*speed**2);ad=(3*md/(4*math.pi*rho_d))**(1/3)
        sig=1e-3*md/(math.pi*(ap+ad)**2)
        # Firstflight dust mean residence4R/(3v), assuming isotropic directions and fixed speed.
        # Sigma=4R Q/(A v) for uniform wall emission and absorption.
        qd=sig*area*100/(4*R)
        dust.append({'packet_radius_m':ap,'rho_packet_kg_m3':rho_p,'relative_speed_m_s':speed,
         'energy_retention_fraction':eta,'per_packet_specific_heat_budget_J_kg':budget,
         'critical_dust_radius_m_rho2000':ad,'critical_grain_mass_kg':md,
         'dust_mass_column_for_1e_minus3_critical_impacts_per_packet_kg_m2':sig,
         'source_kg_s_at_dust_speed100_m_s_unit_sticking':qd,
         'source_kg_s_if_any_of1000_packets_must_fail_less_than1e_minus3':qd/1000,
         'scope':'Conservative single-grain energy criterion; not a crater model or proof that one affected packet spoils a whole pulse. Smaller grains require aggregate heat test.'})
    # Subjective80% inputs: split-log source distribution; logit capture/retention;
    # lognormal effective atomic encounter crosssection. All independent for this diagnostic.
    rng=np.random.default_rng(3102026);n=200000;z=rng.standard_normal((4,n));z90=ndtri(.9)
    G=np.exp(np.where(z[0]<0,z[0]*math.log(10)/z90,z[0]*math.log(400)/z90))
    ss=expit(np.where(z[1]<0,z[1]*math.log(9)/z90,z[1]*math.log(19)/z90))
    xs=1e-19*np.exp(z[2]*math.log(10)/z90);eta=expit(z[3]*math.log(9)/z90)
    cm=base['phase_independent_mass_column_bound_kg_m2']*G/ss
    optical=cm/M_C*xs
    rcmin=3*eta*cm*1e14/(8*1000*1e4)
    uncertainty={'classification':'Subjective conditional uncertainty, not simulation-derived probability; distributional shape and independence are analyst judgments.',
      'acceptance':{'dispersed_projectiles_encounter_fraction':.001,'condensed_packet_extra_heat_J_kg':1e4,
       'scope':'Engineering screening allocations. Actual implosion/packet tolerances are not specified or validated.'},
      'input_central_and80pct':{
       'thermal_source_multiplier':{'central':1.,'p10':.1,'p90':400.,'rationale':'Broad conditional source range covers accommodation belowunity and strong rough/nanoparticle enhancement. Assumes saved thermal profile remains representative; excludes impact-driven material ejection.'},
       'effective_condensing_probability_per_wall_encounter':{'central':.5,'p10':.1,'p90':.95,'rationale':'Purpose-built cool collecting surfaces; unmeasured sticking/reemission/port behavior. Neutral capture must be designed, not assumed instantaneous.'},
       'atomic_interaction_cross_section_m2':{'central':1e-19,'p10':1e-20,'p90':1e-18,'rationale':'Order-of-magnitude atomic collision screen only. Replace with charge/species/energy-resolved data; not fitted to this beam.'},
       'retained_gas_collision_energy_fraction':{'central':.5,'p10':.1,'p90':.9,'rationale':'Macroscopic packet heat partition not simulated. Complete retention1 is also reported as a bounding screen.'}},
      'distribution':'Source multiplier: median1, split normal inlog with10/90percentiles. Capture: split logit-normal median.5. Crosssection: lognormal. Energyretention:logit-normal. Independent draws, seed3102026,200000samples.',
      'derived_quantiles':{'thermal_source_kg_s':quant(q.sum()*G),'phase_mass_column_kg_m2':quant(cm),
       'optical_depth':quant(optical),'radius_min_rho1000_v1e7_heat1e4_m':quant(rcmin)},
      'explicit_limits':'These80% intervals do not include nonthermal sputter/dust/plasma loads or unknown spatialhotspots. No credible interval for total contamination is defensible from currentdata.'}
    output={'scope':__doc__,'inputs':{'results.json':sha(p),'plume.py':sha(here/'plume.py')},'code_sha256':sha(__file__),
     'basis':base,'pulse_bounds':shot,'column_kernel':{'maximum':gmax,'dimensionless_time_at_maximum':math.exp(gfit.x),
       'used_peak_bound':3.26,'exact_tail_inequality':'Gamma(w)<=8/w^4','tail_integral_constant_check':ct,
       'formula':'Gamma(w)=integral_0^min(65,4/w²) x exp(-x) B(w sqrt(x)) dx, B(a)=.5−1/a−(1+a²)/a² log|1−a|,0<a<2.',
       'interpretation':'Instantaneous uniform spherical wall mass source; diffuse cosine directions,flux-Maxwell speeds,no collisions,first wallhitcaptured. Numerical peak guardedby0.3%; pastshotbound usesmin(3.26,8/w^4).',
       'g_grid_w':grid.tolist(),'g_grid':vals.tolist()},
     'gas_scenarios':scenarios,'packet_rows':packet,'dust_acceptance_rows':dust,'uncertainty':uncertainty,
     'acceptance_source_scales':{'uniform_neutral_C_source_kg_s_for_encounter1e_minus3_sigma1e_minus18_unit_capture':(-math.log1p(-.001))*area*vbar[0]*M_C/(4*R*1e-18),'uniform_neutral_C_source_kg_s_for_packet_heat1e4_r1mm_rho1000_v1e7_all_retained':(2*1e4*(4*1000*.001/3)/1e14)*area*vbar[0]/(4*R),'scope':'Steady uniform source only. Source limits scale with capture probability and gas speed; burst and local-port columns need direct integration.'},
     'checks':{'passed':all(checks),'count':len(checks)},
     'limitations':['Thermal data are saved64-grid receiver temperatures; no rerun or totalerosion claim.',
      '1/s scaling of source and periodic bound is a molecular-flow renewal screen. Actual returned gas can cool, change species and stick differently.',
      'Molecular and atom columns aredistinct: atom-count opticaldepth is a conservative additiveinteraction proxy, not a C2/C3 crosssection model.',
      'Surface collisions deposit heat in a stopping layer; bulk specificenergy is necessary but insufficient to ruleout surfaceablation or beam halo.',
      'Crosssection-based dilute-neutral tests do not describe chargedplasma collective fields or micron dust.',
      'No manufacturing energy or cost coefficient is assumed.1000launchers do notdivide theambientcolumn perpacket.']}
    out.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'basis':base,'scenarios':scenarios,'dust':dust,'uncertainty':uncertainty['derived_quantiles'],'checks':output['checks']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[3]);p.add_argument('--output',type=Path);a=p.parse_args()
 main(a.root,a.output or Path(__file__).with_name('gas-assessment.json'))
