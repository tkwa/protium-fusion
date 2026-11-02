#!/usr/bin/env python3
"""Bounded sealed-capsule thermal/pressure/composition screen; no pulse rerun.

Requires CoolProp 8.0.0. Example isolated runtime:
PYTHONPATH=/home/tkwa/src/protium-fusion/output/capsule-runtime \
 /home/tkwa/miniforge3/bin/python screen.py --output screen.json
"""
from pathlib import Path
import argparse, hashlib, json, math
import CoolProp
from CoolProp.CoolProp import PropsSI

HERE=Path(__file__).resolve().parent
COMP=HERE.parent
SIGMA=5.670374419e-8
EPS0=8.8541878128e-12
A=.01
R=1.1e9
V=2.4e6
TW=1000.
ALPHA=.02182
E=1e9
PE=.5*EPS0*E**2
RHO_S=2000.
RHO_AL=2700.
D_AL=100e-9

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p): return json.loads(Path(p).read_text())
def dose(rho,fill=1.,path=1.,outer_ratio=1.):
    # Fixed 1 cm INNER cavity. All radiative/electric energy assigned to core;
    # no shell heat capacity or reradiation credited. Outer radius may be larger.
    ir=3*ALPHA*SIGMA*TW**4*R*path/(rho*fill*A*V)*outer_ratio**2
    electric=1.5*EPS0*E**2/(rho*fill)*outer_ratio**3
    return dict(IR_J_kg_core=ir,electric_J_kg_core=electric,total_J_kg_core=ir+electric)

def shell(pi,allow,rho,fill=1.):
    # Exact thick-sphere linear-elastic maximum hoop stress, with equivalent
    # outward electrostatic surface traction PE. pi+PE is a conservative
    # simultaneous load; the actual maxima generally occur at different times.
    den=allow-.5*pi-1.5*PE
    if den<=0: return dict(valid=False)
    x=(allow+pi)/den
    beta=x**(1/3)
    s=RHO_S/(rho*fill)*(x-1)
    return dict(valid=True,internal_pressure_Pa=pi,allowable_Pa=allow,
      outer_inner_radius_ratio=beta,thickness_m=A*(beta-1),
      carbon_skin_kg_per_core_kg=s,carbon_mass_fraction=s/(1+s),
      aluminum_kg_per_core_kg=3*RHO_AL*D_AL*beta**2/(rho*fill*A),
      reconstructed_max_hoop_Pa=(pi*(1+.5*x)+1.5*PE*x)/(x-1),
      thin_shell_mass_ratio=3*RHO_S*(pi+PE)/(2*rho*fill*allow))

def methane(fill,path=1.,outer_ratio=1.):
    dn=439.*fill
    di=dn*17/16
    u0=PropsSI('U','T',102.4,'D',dn,'Methane')
    q=dose(466.4375,fill,path,outer_ratio)
    un=u0+q['total_J_kg_core']*17/16
    t=PropsSI('T','U',un,'D',dn,'Methane')
    p=PropsSI('P','T',t,'D',dn,'Methane')
    return dict(fill=fill,free_path_fraction=path,**q,temperature_K=t,
      pressure_Pa=p,initial_pressure_Pa=PropsSI('P','T',102.4,'D',dn,'Methane'),
      natural_density_kg_m3=dn,isotope_density_kg_m3=di,
      EOS_within_published_T_P_range=bool(t<=625 and p<=1e9),
      energy_closure_relative=(PropsSI('U','T',t,'D',dn,'Methane')-un)/max(abs(un),1.))

def geometry_corrected_methane(fill,path,allow=1e10):
    # Iterate small outer-area/field-energy corrections; neglecting shell heat
    # capacity remains conservative for a well-equilibrated core.
    beta=1.
    for n in range(100):
        m=methane(fill,path,beta)
        s=shell(m['pressure_Pa'],allow,466.4375,fill)
        if not s['valid']: raise ValueError('No spherical stress solution')
        new=s['outer_inner_radius_ratio']
        if abs(new-beta)<1e-12: break
        beta=.5*(beta+new)
    else: raise ValueError('Geometry iteration failed')
    s['effective_H_per_C_if_C13_skin']=4/(1+s['carbon_skin_kg_per_core_kg']*17/13)
    return dict(**m,shell=s,iterations=n+1,geometry_fixed_point_error=abs(new-beta),
      scope='10 GPa only; fixed inner volume, exact outer area, no shell heat capacity; no thermoelastic strain solution')

def pressure_limited_path(fill,p_limit=30e6):
    dn=439*fill
    t0=102.4
    def p(t): return PropsSI('P','T',t,'D',dn,'Methane')
    lo,hi=t0,625.
    for _ in range(80):
        mid=(lo+hi)/2
        if p(mid)<p_limit: lo=mid
        else: hi=mid
    t=(lo+hi)/2
    du=(PropsSI('U','T',t,'D',dn,'Methane')-PropsSI('U','T',t0,'D',dn,'Methane'))*16/17
    q=dose(466.4375,fill)
    fraction=max(0.,min(1.,(du-q['electric_J_kg_core'])/q['IR_J_kg_core']))
    return dict(fill=fill,pressure_cap_Pa=p_limit,temperature_at_cap_K=t,
      allowed_core_heat_J_kg=du,electric_retention_fraction=1.,
      maximum_free_path_fraction=fraction,maximum_free_path_m=fraction*R,
      scope='Thin geometry reference; pressure cap chosen, not a survival ceiling')

def ethane(fill,allow=1e10):
    # New comparator: saturated natural liquid at100K, isotope volume32/30.
    rho0=PropsSI('D','T',100.,'Q',0,'Ethane')
    dn=rho0*fill;di=rho0*32/30
    u0=PropsSI('U','T',100.,'D',dn,'Ethane')
    beta=1.
    for n in range(100):
        q=dose(di,fill,1.,beta)
        un=u0+q['total_J_kg_core']*32/30
        t=PropsSI('T','U',un,'D',dn,'Ethane')
        p=PropsSI('P','T',t,'D',dn,'Ethane')
        sh=shell(p,allow,di,fill)
        new=sh['outer_inner_radius_ratio']
        if abs(beta-new)<1e-12:break
        beta=.5*(beta+new)
    else:raise ValueError('ethane geometry')
    sh['effective_H_per_C_if_C13_skin']=3/(1+sh['carbon_skin_kg_per_core_kg']*16/13)
    return dict(fill=fill,initial_T_K=100.,natural_liquid_density_kg_m3=rho0,
       isotope_liquid_density_kg_m3=di,**q,temperature_K=t,pressure_Pa=p,shell=sh,
       EOS_within_published_T_P_range=bool(t<=675 and p<=900e6),
       energy_closure_relative=(PropsSI('U','T',t,'D',dn,'Ethane')-un)/max(abs(un),1.),
       isotope_scaling='32/30 density and30/32 specificenergy; no isotope-specific EOS',
       scope='New 100K liquid recipe; no matching FEOS entropy map or hydro run.')

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,default=HERE/'screen.json');args=ap.parse_args()
    checks=[]
    def check(name,ok,value=None):
        checks.append(dict(name=name,passed=bool(ok),value=value))
        if not ok: raise AssertionError(name)
    states={
     'A_NH3':dict(rho=850.,T=150.,M_per_catalyst=17.,fill=.8,pressure=6060.,
       latent_J_mol=5650.,molecular_mass_kg_mol=.017,Tm=194.95),
     'B1_N15hydrazine':dict(rho=1062.5,T=250.,M_per_catalyst=17.,fill=.8,pressure=319.,
       latent_J_mol=12657.,molecular_mass_kg_mol=.034,Tm=274.69),
     'B2_CH2':dict(rho=944.2,T=300.,M_per_catalyst=14.,fill=.8,pressure=1e6),
     'C13_CH2':dict(rho=944.2*15/14,T=300.,M_per_catalyst=15.,fill=.8,pressure=1e6),
    }
    out={}
    for name,s in states.items():
        q=dose(s['rho'],s['fill'])
        row=dict(cold_reference=s,**q,shells={})
        if 'latent_J_mol' in s:
            latent=s['latent_J_mol']/s['molecular_mass_kg_mol']
            row['fusion_energy_J_kg']=latent
            row['melt_fraction_upper_ignoring_sensible_heat']=q['total_J_kg_core']/latent
            row['bulk_equilibrium_T_upper_K']=s['Tm']
            check(name+'_total_dose_below_fusion',q['total_J_kg_core']<latent)
        for strength in [1e9,1e10]:
            sh=shell(s['pressure'],strength,s['rho'],s['fill'])
            sh['C12_atoms_per_native_catalyst']=sh['carbon_skin_kg_per_core_kg']*s['M_per_catalyst']/12
            sh['charge_energy_J_per_kg_carbon_skin']=q['electric_J_kg_core']/sh['carbon_skin_kg_per_core_kg']
            if name in ['B2_CH2','C13_CH2']:
                massC=12 if name=='B2_CH2' else 13
                sh['effective_H_per_C_same_isotope_skin']=2/(1+sh['carbon_skin_kg_per_core_kg']*s['M_per_catalyst']/massC)
            row['shells'][str(int(strength))]=sh
            check(name+'_stress_'+str(strength),abs(sh['reconstructed_max_hoop_Pa']/strength-1)<1e-10)
        out[name]=row
    # Original Overstreet/Giauque Cp values (cal/mol/K); interpolate the final
    # 4.95 K linearly from the two preceding measured smooth-table entries.
    t=[150.,160.,170.,180.,190.,194.95]
    cp=[9.272,9.846,10.42,11.03,11.71,11.71+.068*4.95]
    sensible=sum((t[i+1]-t[i])*(cp[i+1]+cp[i])/2 for i in range(len(t)-1))*4.184/.017
    out['A_NH3']['sensible_150K_to_melt_J_kg_primary_Cp']=sensible
    out['A_NH3']['estimated_melt_fraction_with_primary_Cp']=max(0.,(out['A_NH3']['total_J_kg_core']-sensible)/out['A_NH3']['fusion_energy_J_kg'])
    rhol=PropsSI('D','T',195.5,'Q',0,'Ammonia')
    out['A_NH3']['liquid_density_reference_195p5K']=rhol
    out['A_NH3']['cavity_fraction_if_all_liquid']=.8*850/rhol
    check('NH3_complete_melt_volume_fits_20pct_headspace',.8*850/rhol<1)
    # PE phase model: explicit Cp/Tm/crystallinity screen, no new fitted EOS.
    pe=[]
    for chi in [.25,.5,1.]:
        cp=1600.;tm=400.;lf=281e3
        q=out['B2_CH2']['total_J_kg_core'];qm=cp*(tm-300)
        tf=300+q/cp if q<qm else tm+max(0,q-qm-chi*lf)/cp
        pe.append(dict(initial_crystallinity=chi,Cp_J_kg_K=cp,Tm_K=tm,
           pure_crystal_latent_J_kg=lf,temperature_K=tf,
           melted_mass_fraction=min(chi,max(0,(q-qm)/lf)),
           scope='Piecewise screening model, phase-dependent measured capsule enthalpy and expansion still needed'))
    out['B2_CH2']['polymer_phase_scenarios']=pe
    out['C13_CH2']['polymer_phase_scope']='Same molecular thermodynamics approximately scaled14/15 perkg; use B2 scenarios as adjacent isotope reference, not isotope-specific calorimetry.'
    out['C13_CH2']['precooling_20_to300K_enthalpy_bounds_J_kg']=[269e3*14/15,302e3*14/15]
    out['C13_CH2']['precooling_20_to300K_dose_fraction_of_lower_bound']=out['C13_CH2']['total_J_kg_core']/(269e3*14/15)
    check('precooled_C13polymer_reference_bulk_below300K',out['C13_CH2']['total_J_kg_core']<269e3*14/15)
    meth=[]
    for path in [.2,1.]:
        for fill in [1.,.8,.5]:
            row=methane(fill,path);row['mechanical_only_shells']={}
            for strength in [1e9,1e10]:
                sh=shell(row['pressure_Pa'],strength,466.4375,fill)
                sh['effective_H_per_C_if_C13_skin']=4/(1+sh['carbon_skin_kg_per_core_kg']*17/13)
                sh['natural_C12_atoms_per_core_C13_if_C12_skin']=sh['carbon_skin_kg_per_core_kg']*17/12
                row['mechanical_only_shells'][str(int(strength))]=sh
                check('methane_stress_'+str((path,fill,strength)),abs(sh['reconstructed_max_hoop_Pa']/strength-1)<1e-10)
            row['geometry_corrected_10GPa']=geometry_corrected_methane(fill,path)
            check('methane_EOS_energy_'+str((path,fill)),abs(row['energy_closure_relative'])<1e-8)
            check('methane_domain_'+str((path,fill)),row['EOS_within_published_T_P_range'])
            meth.append(row)
    cpath=COMP/'source-family-carbon-median128/integration.json';cycle=load(cpath)['families']['source128']['cycle']
    gross={}
    for strength in [1e9,1e10]:
        bystage={}
        for name,cat in [('A_NH3','n14'),('B1_N15hydrazine','n15'),('B2_CH2','c12')]:
            bystage[cat]=out[name]['shells'][str(int(strength))]['carbon_skin_kg_per_core_kg']*cycle['stage_pulse_data'][cat]['gross_fuel_kg_per_He']/cycle['fresh_hydrogen_kg_per_He']
        total=sum(bystage.values());reserve=.007-cycle['fresh_C_H_mass_ratio']
        gross[str(int(strength))]=dict(gross_shell_carbon_kg_per_fresh_H_kg=bystage,
          subtotal_excluding_methane=total,
          maximum_terminal_shell_loss_fraction_if_entire_remaining_budget_allocated=min(1,reserve/total),
          scope='Inventory scale only against unchanged source128 fuel flows. Skins change actual inputs; this is not a re-solved isotope ledger or full material allowance.')
    tracerpath=COMP/'isotope-tracer-final8.json';tr=load(tracerpath)
    tracer={s:dict(c12=tr['stages'][s]['cold_tracer']['c12'],native_reconstruction_L1=tr['stages'][s]['native_hot_reconstruction_L1_per_initial_catalyst']) for s in ['A','B1']}
    mix={}
    for h in [3,4]:
        p=COMP.parent/f'mixture-design/rate-results/carbon128/C13-H{h}-fixed_energy-z128-k0p001-t1-rates-updated-median.json'
        j=load(p);mix[str(h)]=dict(path=str(p),sha256=sha(p),gain=j['nuclear_energy_J']/j['initial_energy_J'],success=j['success'],initial_energy_J=j['initial_energy_J'])
    qf=dose(1000.)
    check('reference_IR',abs(qf['IR_J_kg_core']/170125.4085412591-1)<1e-9)
    check('field_pressure_energy_identity',abs(qf['electric_J_kg_core']/(3*PE/1000)-1)<1e-12)
    # Verify nonzero-inner-pressure + outward-traction mechanics by two limits.
    sh=shell(0.,1e10,1000.)
    check('thin_shell_limit',abs(sh['carbon_skin_kg_per_core_kg']/sh['thin_shell_mass_ratio']-1)<.001)
    for f in [.5,.8,1.]:
        check('dose_inverse_filling_'+str(f),abs(dose(1000,f)['total_J_kg_core']*f/qf['total_J_kg_core']-1)<1e-12)
    neutral=[]
    for strength in [1e9,1e10]:
        u=2*strength/RHO_S
        neutral.append(dict(allowable_Pa=strength,field_energy_per_skin_mass_thin_charge_limited_J_kg=u,
          shell_sensible_allowance_J_kg=1e6,maximum_skin_deposition_fraction_for_allowance=min(1,1e6/u),
          scope='Illustrative cp1000 J/kg/K and deltaT1000 K; no inferred electron energy-deposition fraction'))
    eth=[ethane(f) for f in [1.,.8,.5]]
    for z in eth:
        check('ethane_energy_'+str(z['fill']),abs(z['energy_closure_relative'])<1e-8)
        check('ethane_domain_'+str(z['fill']),z['EOS_within_published_T_P_range'])
    rc=load(cpath)['families']['source128']['stage_inputs']['c13']['selected']['receiver']
    ec=rc['energy_channels_J'];tubes=[]
    for g in [.2,.5]:
        rt=R*g;nt=1e5;b=250.;area=2*math.pi*b*(R-rt)*nt;af=nt*b*b/(4*rt*rt)
        tubes.append(dict(free_path_fraction=g,tip_radius_m=rt,outlet_count=nt,bore_radius_m=b,
          projected_intercepted_pulse_fraction=af,total_outer_cylindrical_area_m2=area,
          tube_area_per_receiver_area=area/(4*math.pi*R*R),
          C13_kinetic_fluence_J_m2=ec['outgoing_kinetic_nodal']/(4*math.pi*rt*rt),
          C13_advected_internal_upper_fluence_J_m2=ec['endpoint_internal']/(4*math.pi*rt*rt),
          C13_prompt_fluence_J_m2=ec['escaped_prompt_radiation']/(4*math.pi*rt*rt),
          C13_peak_kinetic_flux_W_m2=rc['ballistic']['peak_kinetic_flux_W_m2']*(rc['radius_m']/rt)**3,
          C13_kinetic_intercepted_energy_J=ec['outgoing_kinetic_nodal']*af,
          mass_equivalent_per_receiver_m2_at10_and100kg_per_tube_m2=[area/(4*math.pi*R*R)*s for s in [10.,100.]],
          insulation_heat_per_receiver_m2_at1_and10W_per_tube_m2=[area/(4*math.pi*R*R)*q for q in [1.,10.]],
          scope='Radially aligned tip/disks only, no claim that sidewalls, permanent holes, shutter timing, support or local stopping are qualified. Hot exterior/cold interior with separately budgeted insulation.'))
    sourcepaths=[cpath,tracerpath,COMP/'launcher-mixture-entropy.json',COMP/'packet-radiation.json']
    result=dict(schema='capsule-survival-v1',scope='Bounded sealed capsule thermodynamics and elastic stress, no material qualification or changed pulse simulation',
      constants=dict(inner_radius_m=A,wall_T_K=TW,R_m=R,speed_m_s=V,flight_time_s=R/V,
      absorptivity_reference=ALPHA,surface_field_V_m=E,electrostatic_pressure_Pa=PE,
      carbon_shell_density_kg_m3=RHO_S,Al_optical_coat_m=D_AL),
      nonmethane=out,methane=meth,ethane_100K_comparator=eth,protected_tube_examples=tubes,pressure_limited_methane_paths=[pressure_limited_path(f) for f in [1.,.8,.5]],
      shell_inventory_scale=gross,legacy_passive_tracers=tracer,neighboring_nominal128_C13_mixtures=mix,
      neutralization_local_heat=neutral,
      thermal_gradient_reference=dict(absorbed_surface_flux_W_m2=ALPHA*SIGMA*TW**4,
        core_diffusivity_assumed_m2_s=1e-7,diffusion_length_flight_m=math.sqrt(1e-7*R/V),
        diffusion_length_0p33ms_m=math.sqrt(1e-7*.00033),
        scope='Diffusivity is an example; headspace/contact thermal resistance unresolved'),
      provenance=dict(CoolProp_version=CoolProp.__version__,CoolProp_git_revision=CoolProp.__gitrevision__,
        isotope_scaling='Natural-methane EOS with rho13/rhonatural=17/16 and u13/unatural=16/17, matching stored source convention; no isotope-specific EOS',
        methane_EOS='Setzmann/Wagner 1991, doi:10.1063/1.555898',
        sources={str(p):sha(p) for p in sourcepaths},code_sha256=sha(__file__)),checks=checks,all_checks_pass=all(c['passed'] for c in checks))
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(output=str(args.output),checks=len(checks),all_pass=result['all_checks_pass'],methane_full=[dict(fill=m['fill'],T=m['temperature_K'],P_MPa=m['pressure_Pa']/1e6,geometry_corrected_skin=m['geometry_corrected_10GPa']['shell']['carbon_skin_kg_per_core_kg']) for m in meth if m['free_path_fraction']==1],nitrogen_shell_inventory=gross),indent=2))
if __name__=='__main__':main()
