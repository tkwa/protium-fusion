"""Conditional complete-mass allocation with explicit isotope processors.

Uses the already derived winding/finite-field model and a finite equilibrium
family. Nuclear side reactions replace the previous generic restoration cost.
Conductor capability, driven extraction, electron clearing, 3-D mechanics,
shield transport and service life remain engineering assumptions.
"""
from pathlib import Path
import json,gzip
from math import pi,sqrt
from copy import deepcopy
from bundle_solenoid import winding,finite_axis_factor,YEAR,SB
from nonneutral_screen import MEV,C,EPS0,Q,DD_CHARGED_PER_NEW_N,DD_THERMAL_PER_NEW_N
from catalyst_catalog import atomicmass,nucmass
from capture_bank_screen import screen
from allocate_network_channels import choose
ROOT=Path(__file__).parent


def prepare(energy,family,B,heat=.7,heater=.95,direct=.95,potential_in=.999,potential_out=.999,
            end_barrier_kT=25,end_in=.999,end_out=.999,extraction=.2,atomic=.1):
    e=deepcopy(energy);fam={'configurations':{}}
    g=direct+(1-direct)*heat
    e['bare_electricity_MeV_per_primary']=g*e['charged_MeV_per_primary']+heat*e['thermal_MeV_per_primary']-e['input_thermal_MeV_per_primary']/heater-e['network']['discard_electricity_reserve_1000MeV_per_event_MeV']
    for key,rows in family['configurations'].items():
        values=[]
        s=next(s for s in e['stages'] if str(1000*s['ZA'][0]+s['ZA'][1])==key);nex=s['extraction_events_per_pn']
        for orig in rows:
            r=orig.copy();scale=B/r['B_T'];r['B_T']=B
            r['radius_m']/=scale;r['reaction_rate_per_m_s']*=scale**2;r['area_per_pn_rate_m3_s']/=scale**4
            r['wall_E_V_m']*=scale;r['wall_circulation_MeV_per_pn']/=scale
            rot=r['rotation_input_MeV_per_pn'];nr=r['neutron_rotation_MeV_per_pn']
            r['rotation_loss_MeV_per_pn']=rot/heater-g*(rot-nr)-heat*nr
            r['wall_loss_MeV_per_pn']=r['wall_circulation_MeV_per_pn']*(1/heater-g)
            r['potential_loss_MeV_per_pn']=r['potential_circulation_MeV_per_pn']*(1/potential_in-potential_out-(1-potential_out)*heat)
            r['radiation_loss_MeV_per_pn']=(r['Larmor_MeV_per_pn']+r['ion_brems_MeV_per_pn'])*(1/heater-heat)
            r['end_circulation_MeV_per_pn']=nex*(r['ZA'][0]+1)*r['T_MeV']*end_barrier_kT
            r['end_loss_MeV_per_pn']=r['end_circulation_MeV_per_pn']*(1/end_in-end_out-(1-end_out)*heat)
            r['extraction_loss_MeV_per_pn']=nex*extraction;r['atomic_preparation_reserve_MeV_per_pn']=nex*atomic
            r['total_loss_MeV_per_pn']=sum(r[k] for k in ['rotation_loss_MeV_per_pn','wall_loss_MeV_per_pn','potential_loss_MeV_per_pn','radiation_loss_MeV_per_pn','end_loss_MeV_per_pn','extraction_loss_MeV_per_pn','atomic_preparation_reserve_MeV_per_pn'])
            values.append(r)
        fam['configurations'][key]=values
    return e,choose(e,fam)


def plant(energy,family,B=140,heat=.7,core_radius=1000,shield=5,packing=.7,J=1e8,
          sigma_hoop=16e9,rho_hoop=1600,sigma_comp=8e9,rho_comp=3500,Ecomp=1e12,
          comp_mass_factor=1.3,end_guard_radii=3,end_system_fraction=.15,
          euler_safety=2,global_margin=1.2,stress_poisson=.2,
          heater=.95,direct=.95,potential_in=.999,potential_out=.999,
          end_barrier_kT=25,end_in=.999,end_out=.999,extraction=.2,atomic=.1,
          neutron_capture=.995,neutron_leak=1e-6,gamma_leak=1e-5,
          coil_cold_K=4.2,radiator_K=30,hot_K=300,ambient_K=3,
          radiator_areal_mass=1,radiator_emissivity=.9,cryo_Carnot_fraction=.3,
          warm_to_cold_W_m2=.05,outer_to_cold_W_m2=.001,
          J_fraction_of_critical=.1,power_law_n=30,Ec_V_m=1e-4,
          warm_water_density=1000,wall_E=1e12,wall_density=3500,water_pressure=1e4,wall_buckling_safety=3,
          conversion_kg_W=.1,Dplant_kg_W=.2,processing_kg_W=.1,other_fraction=.1,
          capture_stock_turnover_years=1000,capture_processing_J_kg=1e8,
          extra_processing_MeV_per_primary=0.,
          annual_replacement_fraction=1e-4,embodied_J_kg=1e8,
          manufacturing_kg_W=1,replacement_inventory_years=10):
    assumptions={k:v for k,v in locals().items() if k not in ['energy','family']}
    a=core_radius+shield;w=winding(B,a,J,sigma_hoop,rho_hoop)
    F=w['field_energy_J_per_m'];area_comp=F/sigma_comp;rout=w['outer_radius_m']
    rcomp=sqrt(rout*rout+area_comp/pi);IoverA=(rcomp*rcomp+rout*rout)/4
    L=pi*sqrt(Ecomp*IoverA/(euler_safety*sigma_comp*global_margin))
    activeL=max(0,L-2*end_guard_radii*rout);finite=finite_axis_factor(w,L,activeL/2)
    e,allocation=prepare(energy,family,B*finite,heat,heater,direct,potential_in,potential_out,end_barrier_kT,end_in,end_out,extraction,atomic)
    rows=allocation['selected'];flows=e['network']['processor_pn_rates'];network=e['network']
    def total(key):return sum(flows[z]*r[key] for z,r in rows.items())
    rate=packing*pi*core_radius**2*activeL/allocation['channel_volume_per_primary_rate_m3_s']
    base=rate*allocation['net_source_MeV_per_primary']*MEV
    g=direct+(1-direct)*heat;usefuln=network['net_neutrons_available_for_light_fuels']-network['light_net_before_D_completion'].get('2003',0)
    lostn=(1-neutron_capture)*usefuln
    capture_penalty=rate*lostn*(g*DD_CHARGED_PER_NEW_N+heat*(DD_THERMAL_PER_NEW_N+Q['np_capture']))*MEV
    qn=e['neutron_MeV_per_primary']+total('neutron_rotation_MeV_per_pn')
    qg=e['prompt_photon_MeV_per_primary']+network['cold_capture_Q_MeV']+usefuln*neutron_capture*Q['np_capture']
    cold_radiation=rate*(neutron_leak*qn+gamma_leak*qg)*MEV
    inner_area=2*pi*a*L+2*pi*a*a;outer_area=2*pi*rout*L+2*pi*rout*rout
    inner_leak=warm_to_cold_W_m2*inner_area;outer_leak=outer_to_cold_W_m2*outer_area
    creep=Ec_V_m*J_fraction_of_critical**power_law_n*J*w['conductor_mass_kg_per_m']/8000*L
    cold_load=cold_radiation+inner_leak+outer_leak+creep
    fridge=cold_load*(radiator_K/coil_cold_K-1)/cryo_Carnot_fraction
    lost_heat=heat*(cold_radiation+inner_leak)
    wall_volume=0.;channel_rows=[];active_ions=0.
    for z,r in rows.items():
        rr=r['radius_m'];channels=rate*flows[z]/(activeL*r['reaction_rate_per_m_s'])
        pressure=water_pressure+EPS0*r['wall_E_V_m']**2/2
        thickness=rr*(4*(1-stress_poisson**2)*pressure*wall_buckling_safety/wall_E)**(1/3)
        volume=channels*pi*((rr+thickness)**2-rr*rr)*activeL;wall_volume+=volume
        za=r['ZA'];active_ions+=rate*flows[z]/r['reaction_rate_per_m_s']*(r['proton_inventory_per_m']*nucmass(1,1)+r['catalyst_inventory_per_m']*nucmass(*za))*MEV/C**2
        channel_rows.append({'ZA':za,'pn_rate_s':rate*flows[z],'channels':channels,'radius_m':rr,'wall_m':thickness,'wall_pressure_Pa':pressure,'wall_volume_m3':volume,'volume_fraction_of_channels':allocation['area_fractions'][z]})
    core_vol=pi*core_radius**2*activeL
    water_vol=max(0,core_vol*(1-packing)-wall_volume)+pi*(a*a-core_radius**2)*L+2*pi*a*a*shield
    water_mass=warm_water_density*water_vol;wall_mass=wall_volume*wall_density
    nh=2*warm_water_density/((atomicmass(8,16)+2*atomicmass(1,1))*MEV/C**2)
    sigma_h=nh*.3326e-28;flux=rate*usefuln*neutron_capture/(sigma_h*water_vol)
    bank=screen(network,rate,flux);stock=bank['total_listed_inventory_lower_bound_kg']
    stock_processing=stock*capture_processing_J_kg/(capture_stock_turnover_years*YEAR)
    extra_processing=rate*extra_processing_MeV_per_primary*MEV
    preaux=base-capture_penalty-lost_heat-fridge-creep-stock_processing-extra_processing
    other=other_fraction*max(0,preaux);bus=preaux-other
    comp=area_comp*rho_comp*L*comp_mass_factor;hoop=w['hoop_mass_kg_per_m']*L;conductor=w['conductor_mass_kg_per_m']*L
    sigma_shell=Ecomp*(rcomp-rout)/(sqrt(3*(1-stress_poisson**2))*rcomp)
    decay=rate*sum(r['kg_per_primary_pn_s'] for r in network['decay_inventories'])
    qnuclear=e['nuclear_heat_MeV_per_primary'];heat_source=rate*qnuclear*MEV
    nuclear_after_capture=rate*(qnuclear-lostn*(DD_CHARGED_PER_NEW_N+DD_THERMAL_PER_NEW_N+Q['np_capture']))*MEV
    grosscirculation=e['input_thermal_MeV_per_primary']+qnuclear+extra_processing_MeV_per_primary+total('rotation_input_MeV_per_pn')+total('wall_circulation_MeV_per_pn')+total('potential_circulation_MeV_per_pn')+total('end_circulation_MeV_per_pn')
    converter=conversion_kg_W*rate*grosscirculation*MEV
    # All downstream light-fuel heat excluding its initial H-capture photons.
    Dgross=rate*(network['downstream_light_Q_MeV']-usefuln*Q['np_capture']-lostn*Q['net_D']/2)*MEV
    masses={'circumferential_structure':hoop,'conductor':conductor,'axial_compression_support':comp,
      'water_moderator_and_shield':water_mass,'channel_walls':wall_mass,'cold_decay_inventory':decay,
      'capture_isotope_stock_lower_bound':stock,'active_plasma_ions':active_ions,
      'conversion_and_potential_handling':converter,'closed_D_plant':Dplant_kg_W*Dgross,
      'isotope_processing':processing_kg_W*heat_source,'capture_bank_processing_equipment':manufacturing_kg_W*stock_processing,
      'end_system_allowance':end_system_fraction*(hoop+conductor+comp)}
    fixed=sum(masses.values());radflux=radiator_emissivity*SB*(radiator_K**4-ambient_K**4);radkgW=radiator_areal_mass/radflux
    k=annual_replacement_fraction*embodied_J_kg/YEAR;spare=annual_replacement_fraction*replacement_inventory_years
    denominator=1-k*(radkgW+manufacturing_kg_W)-spare
    if denominator<=0:raise ValueError('No finite maintenance closure')
    mass=(fixed+(nuclear_after_capture-bus)*radkgW)/denominator
    maintenance=k*mass;net=bus-maintenance;waste=nuclear_after_capture-net
    masses.update({'radiators':waste*radkgW,'replacement_manufacturing_equipment':manufacturing_kg_W*maintenance,'replacement_inventory':spare*mass})
    return {'scope':__doc__,'assumptions':assumptions,'net_W':net,'total_mass_kg':mass,'net_W_kg':net/mass,
      'plant_bus_net_W':bus,'plant_bus_mass_kg':fixed+(nuclear_after_capture-bus)*radkgW,
      'primary_pn_rate_s':rate,'source_power_W':base,'source_net_MeV_per_primary':allocation['net_source_MeV_per_primary'],
      'source_field_T':B*finite,'length_m':L,'active_length_m':activeL,'outer_winding_radius_m':rout,
      'primary_volume_per_rate_m3_s':allocation['channel_volume_per_primary_rate_m3_s'],
      'channel_allocations':channel_rows,'selected_source_profiles':rows,
      'source_loss_MeV_per_primary':{key:total(key) for key in ['rotation_loss_MeV_per_pn','wall_loss_MeV_per_pn','potential_loss_MeV_per_pn','radiation_loss_MeV_per_pn','end_loss_MeV_per_pn','extraction_loss_MeV_per_pn','atomic_preparation_reserve_MeV_per_pn']},
      'capture_penalty_W':capture_penalty,'cold_radiation_W':cold_radiation,'refrigeration_W':fridge,'total_cold_load_W':cold_load,
      'neutron_source_energy_MeV_per_primary':qn,'prompt_and_capture_gamma_MeV_per_primary':qg,
      'capture_stock_processing_W':stock_processing,'extra_processing_W':extra_processing,'other_auxiliary_W':other,'replacement_manufacturing_W':maintenance,
      'nuclear_after_capture_W':nuclear_after_capture,'waste_heat_W':waste,'radiator_area_m2':waste/radflux,
      'conversion_gross_circulation_MeV_per_primary':grosscirculation,'mass_breakdown_kg':masses,
      'capture_bank':bank,'water_volume_m3':water_vol,'H_macroscopic_capture_m_inv':sigma_h,
      'global_buckling_margin':global_margin,'local_buckling_margin':sigma_shell/(euler_safety*sigma_comp),
      'stored_field_energy_J':F*L,'jupiter_masses':mass/1.89813e27,
      'embodied_energy_per_bus_power_years':mass*embodied_J_kg/bus/YEAR if bus>0 else None,
      'mass_ledger_error_kg':sum(masses.values())-mass,'energy_ledger_error_W':nuclear_after_capture-net-waste}

if __name__=='__main__':
    energy=json.loads((ROOT/'network-ion-budget.json').read_text());family=json.loads(gzip.decompress((ROOT/'network-channel-family.json.gz').read_bytes()))
    rows=[]
    for B in [100,120,140,160,180]:
        for heat in [.6,.7]:
            for eta in [.999,.9999]:
                r=plant(energy,family,B=B,heat=heat,end_in=eta,end_out=eta);rows.append(r)
                print('B',B,'heat',heat,'eta',eta,'mW/kg',r['net_W_kg']*1e3,'TW',r['net_W']/1e12,'mass',r['total_mass_kg'],'sourceMeV',r['source_net_MeV_per_primary'],flush=True)
    (ROOT/'network-plant.json').write_text(json.dumps({'scope':__doc__,'cases':rows},indent=2,allow_nan=False)+'\n')
