#!/usr/bin/env python3
"""Local membrane, pressure closure and cooling requirements, not a qualified design.
Run python3 mechanics.py; uses exact timed outgoing ion stopping plus existing station cadence.
"""
from pathlib import Path
import hashlib,json,math
import numpy as np
P=Path(__file__).resolve().parent
SB=5.670374419e-8;C=299792458.;YEAR=365.25*86400
MATS={
'silica':{'rho':2200.,'cp':1000.,'k':1.3,'E':72.7e9,'nu':.16,'alpha':.5e-6,'strength_requirement_Pa':10e6,'property_scope':'Bulk room-temperature Corning k/E/nu; 10MPa is a required allowable for actual film, not measured hot irradiated strength.'},
'carbon':{'rho':2000.,'cp':1000.,'k':100.,'E':50e9,'nu':.2,'alpha':1e-6,'strength_requirement_Pa':10e6,'property_scope':'Controlled in-plane graphite k100 scenario. E50GPa bulk Pyroid reference; nu.2 and alpha1ppm are scenario values, actual anisotropy and irradiation require qualification.'},
'BN':{'rho':2100.,'cp':1000.,'k':60.,'E':30e9,'nu':.2,'alpha':2e-6,'strength_requirement_Pa':10e6,'property_scope':'Momentive bulk in-plane k60. E30GPa and nu.2 scenario only; no hot-film strength implied.'},
'SiC':{'rho':3210.,'cp':1000.,'k':34.8,'E':420e9,'nu':.17,'alpha':4e-6,'strength_requirement_Pa':10e6,'property_scope':'Saint-Gobain bulk k at1000C and E atRT, CTE to1100C; nu.17 controlled scenario. 10MPa remains film requirement.'}}
SOURCES={
'silica':'https://www.corning.com/media/worldwide/csm/documents/1e624b3034474193b5b00ea6f558dd3d2.pdf',
'SiC':'https://www.ceramicsrefractories.saint-gobain.com/files/5681/download',
'BN':'https://dev.momentivetech.com/products/ceramics/ceramic-coatings/pyrolytic-boron-nitride-pbn-coating-and-shapes',
'BN_current_use_limit':'https://www.momentivetech.com/materials/ceramic-materials/pyrolytic-boron-nitride-pbn',
'carbon':'https://www.mineralstech.com/docs/default-source/refractories-documents/pyrogenics/technical-data-sheets/physical-properties-for-pyroid-ht.pdf'}
def sha(x):return hashlib.sha256(Path(x).read_bytes()).hexdigest()
def cap_stress(p,a,t,d):return p*a*(1+d*d)/(4*t*d)
def flatcap(p,a,t,m):
    # Spherical-cap meridional extension, zero pretension; controlled approximation.
    lo=1e-8;hi=1.
    for _ in range(100):
        d=(lo+hi)/2;eps=(1+d*d)/d*math.atan(d)-1;sig=MATS[m]['E']/(1-MATS[m]['nu'])*eps
        pressure=4*sig*t*d/(a*(1+d*d))
        if pressure<p:lo=d
        else:hi=d
    return {'sag_ratio':d,'sag_m':a*d,'tension_stress_Pa':sig,'meridional_strain':eps,'pressure_relative_residual':pressure/p-1}
def coolspan(k,t,q,dt):return math.sqrt(8*k*t*dt/q) if q>0 else None

def main():
    pl=json.loads((P/'pulse-loading.json').read_text());st=json.loads((P/'stopping-results.json').read_text());op=P.parent/'open-loop-redesign/operating-cases.json';base=json.loads(op.read_text())['cases']['high_source_rating_surface']
    R=pl['radius_m'];area=4*math.pi*R*R
    frequency={n:1/base['pulse_intervals_s'][k] for n,k in [('C13','c13'),('A','n14'),('B1','n15'),('B2','c12')]}
    ionrows=[dict(x,stage='C13') for x in st['direct_rows'] if x['material'] in MATS]+[x for x in st['other_stages_direct_rows'] if x['material'] in MATS]
    rows=[];checks=[]
    for r in ionrows:
        m=MATS[r['material']];t=r['thickness_nm']*1e-9;q=r['peak_projectile_loss_flux_W_m2'];stage=r['stage'];pk=pl['stages'][stage]['smooth']['peak_kinetic_power_W_m2']
        for eta in [1.,.1,.01]:
            heat=eta*q;L=coolspan(m['k'],t,heat,100.)
            # Illustrative 1%-obscuration mechanical grid, reusable seals/actuators separately.
            width_ratio=.005;f=2*width_ratio-width_ratio**2;rho_grid=3000.;sigma_grid=1e8
            pressure=pl['stages'][stage]['smooth']['peak_momentum_flux_Pa']
            depth_stress=L*math.sqrt(3*pressure/(4*width_ratio*sigma_grid));depth=max(.001,depth_stress)
            # Sum separate maxima is deliberately conservative: they need not be simultaneous.
            total_q=(1-f)*heat+f*pk
            absorbedF=(1-f)*eta*r['projectile_energy_loss_J_m2']+f*pl['stages'][stage]['totals']['kinetic_fluence_smooth_J_m2']
            rows.append({'stage':stage,'material':r['material'],'thickness_nm':r['thickness_nm'],'retained_fraction_scenario':eta,'projectile_removed_peak_W_m2':q,'local_retained_peak_W_m2':heat,'strip_span_for_100K_m':L,'strip_span_for_1000K_m':L*math.sqrt(10),'lateral_diffusion_time_s':L*L*m['rho']*m['cp']/(math.pi**2*m['k']),'through_thickness_temperature_drop_K_if_backed':heat*t/m['k'],'minimum_contact_conductance_W_m2_K_for_100K':heat/100,'no_cooling_sensible_equivalent_K':eta*r['projectile_energy_loss_J_m2']/(m['rho']*t*m['cp']),'film_areal_mass_kg_m2':m['rho']*t,'grid_width_at_1pct_obscuration_m':width_ratio*L,'grid_obscuration_fraction':f,'grid_depth_stress_floor_m':depth_stress,'grid_assumed_depth_m':depth,'grid_material_floor_kg_m2':rho_grid*f*depth,'grid_plus_film_peak_to_cool_upper_W_m2':total_q,'coolant_inventory_for_10m_loop_10m_s_cp1000_delta100_kg_m2':2*total_q*10/(10*1000*100),'cooling_hardware_requirement_at_1kW_kg_peak_kg_m2':total_q/1000,'per_stage_mean_captured_film_grid_heat_W_m2':absorbedF*frequency[stage],'conditional_pumping_W_m2_if1pct_heat_work':.01*absorbedF*frequency[stage],'scope':'All absorbed heat remains in station heat/rejection balance. Frame interception and local heat routing added; cooling specific power, loop, width/depth and pumping are explicit requirements, not validated hardware.'})
            checks.append(abs(heat*L*L/(8*m['k']*t)/100-1)<1e-12)
    pres=[]
    for nm in [10,1000]:
        for p in [.1,1.,100.,1000.]:
            for a in [.00005,.005,.5]:
                t=nm*1e-9;flat=flatcap(p,a,t,'silica')
                pres.append({'pressure_Pa':p,'thickness_nm':nm,'tile_radius_m':a,'preformed_cap_sag_ratio':.05,'preformed_cap_stress_Pa':cap_stress(p,a,t,.05),'flat_initial_film':flat})
                checks.append(abs(flat['pressure_relative_residual'])<1e-4)
    closures=[]
    for p in [.1,1.,100.,1000.]:
        for gap in [.01,1.,10.]:
            closures.append({'pressure_Pa':p,'gap_m':gap,'global_unbalanced_shell_tension_mass_floor_kg_m2':3000*p*R/(2e8),'local_closed_cassette_tie_mass_floor_kg_m2':3000*p*gap/1e8,'support_density_kg_m3':3000,'support_allowable_Pa':1e8,'scope':'Opposing pressure forces must close locally. Tile edges alone do not remove global shell load; buckling, seals, attachments and thermal distortion omitted from mass floor.'})
    ir=[]
    for tb in [1000,1500,2000]:
        for eff in [.001,.01,.1,1.]:
            q=eff*SB*(tb**4-1000**4)
            ir.append({'film_temperature_K':1000,'back_bath_K':tb,'net_coupling_coefficient':eff,'inward_heat_to_film_W_m2':q,'scope':'Uniform closed front cavity has zero steady net sink. Effective coefficient must come from actual two-sided optical/geometry solution; it is not automatically emissivity.'})
    allocations=[]
    for sigma in [1.,10.,100.]:
        for fp in [1e-6,1e-4,.01]:
            allocations.append({'reusable_frame_seal_cooling_allowance_kg_m2':sigma,'port_fraction_scenario':fp,'shutter_100kg_per_open_area_kg_m2':100*fp,'whole_station_added_mass_kg':area*(sigma+100*fp),'scope':'Allocation ceiling only. Actual 1000-port geometry and accelerator are supplied separately; film not placed across prepared inflow.'})
    duty=[]
    for nm in [10.,1000.]:
        rr=[x for x in rows if x['material']=='silica' and x['thickness_nm']==nm and x['retained_fraction_scenario']==1]
        duty.append({'material':'silica','thickness_nm':nm,'heat_all_stages_mean_W_m2':sum(x['per_stage_mean_captured_film_grid_heat_W_m2'] for x in rr),'pump_if1pct_mean_W_m2':sum(x['conditional_pumping_W_m2_if1pct_heat_work'] for x in rr),'max_cooling_hardware_at1kW_kg_peak_kg_m2':max(x['cooling_hardware_requirement_at_1kW_kg_peak_kg_m2'] for x in rr),'scope':'Outgoing projectile-loss and opaque frame interception only. No IR back-bath burden, incident electrons/sheath, prompt or replacement cost included; not total plant auxiliary.'})
    summary={'scope':'Conditional mechanical and cooling requirements for tiled membrane/cassette alternatives. No radiation-damaged film lifetime or sealed pressure-vessel qualification.','input_sha256':{x:sha(P/x) for x in ['pulse-loading.json','stopping-results.json']},'station_cadence_source_sha256':sha(op),'code_sha256':sha(__file__),'radius_m':R,'area_m2':area,'frequency_per_s':frequency,'material_assumptions':MATS,'material_sources':SOURCES,'pressure_rows':pres,'pressure_closure_rows':closures,'cooling_rows':rows,'duty_summary':duty,'closed_cavity':{'max_return_delay_s':2*R/C,'mean_return_delay_s':4*R/(3*C),'opaque_uniform_front_steady_cooling_W_m2':0.,'return_kernel':'K(t)=2t/tmax², 0<t<tmax, unit integral for diffuse emission on sphere.','back_bath_rows':ir},'allocation_scenarios':allocations,'clamp_thermal_stress_Pa_per500K':{m:p['E']*p['alpha']*500/(1-p['nu']) for m,p in MATS.items()},'checks':{'count':len(checks),'passed':int(sum(checks)),'all_passed':all(checks)}}
    (P/'mechanics.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps({'checks':summary['checks'],'duty_summary':duty,'silica_C13_rows':[x for x in rows if x['material']=='silica' and x['stage']=='C13' and x['thickness_nm'] in [10,1000] and x['retained_fraction_scenario']==1]},indent=2))
if __name__=='__main__':main()
