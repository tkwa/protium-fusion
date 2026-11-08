#!/usr/bin/env python3
"""Flow, thermal-transport and maintenance acceptance scales for the 16x case."""
import csv
import json
import math
from pathlib import Path
import plant
P=Path(__file__).resolve().parent

def loop(L,q,T,dT,u,p=2e6,stress=1e8):
    cp=2.5*8.314462618/.004002602;rs=8.314462618/.004002602;rho=p/(rs*T)
    flow=q*L*L/(cp*dT);area=flow/(rho*u);diam=math.sqrt(4*area/math.pi)
    length=2*L;wall=p*diam/(2*stress)
    mass=math.pi*((diam/2+wall)**2-(diam/2)**2)*length*8000+area*length*rho
    dp=(.015*length/diam+6)*rho*u*u/2
    return dict(territory_m=L,header_kg_m2=mass/(L*L),helium_kg_s_m2=flow/(L*L),
        bore_m=diam,wall_m=wall,pump_W_m2=flow/rho*dp/.8/(L*L),
        mean_K=T,delta_K=dT,allowable_stress_Pa=stress,pressure_Pa=p)

def main():
    x=json.loads((P/'plant-results.json').read_text());c=x['cases'][0];q=c['energy_W_m2']['engine_hot']
    qc=c['energy_W_m2']['waste']+11.0254
    loops=[]
    for size in [30,100,300,1000]:
        h=loop(size,q,1570.,20.,10.);cold=loop(size,qc,420.,20.,5.)
        loops.append(dict(territory_m=size,hot=h,cold=cold,
            combined_kg_m2=h['header_kg_m2']+cold['header_kg_m2'],
            external_pump_W_m2=h['pump_W_m2']+cold['pump_W_m2']))
    # Keep macro rework, selective damage repair and source preparation distinct.
    rate=c['flows_kg_s']['gross_fuel_kg_s'];cm=c['flows_kg_s']['recovered_CNO_kg_s']
    processing=[]
    for work in [1e6,1e8,1e10]:
        for capacity in [1e5,1e8]:
            processing.append(dict(J_kg=work,equipment_kg_per_kg_s=capacity,
                circulating_fuel_W=rate*work,processing_mass_kg=rate*capacity,
                processing_kg_m2=rate*capacity/plant.AREA,
                fraction_export=rate*work/c['net_W']))
    events=list(csv.DictReader((P/'large-n15-16x-events.csv').open()))
    counts={n:sum(e['stage']==n for e in events) for n in ['n14','n15','c12','c13']}
    # Baseline NRT proxies weighted by source mass, with unchanged per-ion damage law.
    dpa={'n14':.0295,'n15':.000319*100,'c12':.0316,'c13':.0197}
    annual=sum(counts[n]*dpa[n] for n in counts)*365.25/370
    layer=[]
    shots=len(events)*365.25/370
    for depth,events_per_year,label in [(1e-7,shots,'100nm face per shot'),(1.5e-4,12,'150um layer monthly')]:
        flow=2000*depth*events_per_year/plant.YEAR
        layer.append(dict(operation=label,flow_kg_s_m2=flow,electric_W_m2=2e6*flow,
                          equipment_kg_m2=1e8*flow))
    out=dict(scope=__doc__,loops=loops,processing=processing,maintenance=dict(
        event_counts_370days=counts,NRT_proxy_dpa_year=annual,
        selected_cap_proxy_10dpa_interval_days=365.25*10/annual,
        operations=layer,qualification='NRT events measure opportunities for damage; creep, cracking, irradiation conductivity, connected pores and annealing require a material model.'),
        fuel=dict(gross_kg_s=rate,catalyst_return_kg_s=cm,hydrogen_return_kg_s=c['flows_kg_s']['recovered_H_kg_s'],
                  ready_CNO_kg_by_hold_days={str(days):cm*days*86400 for days in [1,30,182.625,365.25]},
                  catalyst_loss_kg_s_by_eventual_fraction={str(e):cm*e for e in [1e-6,1e-4,1e-3]}),
        capsule=dict(black_cavity_dose_ratio_1600_to1000=(1600/1000)**4,
            reference_radius_m=.01,equal_dose_radius_m=.01*(1600/1000)**4,
            reference_path_m=1.1e9,equal_dose_exposed_path_m=1.1e9/(1600/1000)**4,
            absorptivity_ceiling_relative_to_old=.022/(1600/1000)**4,
            qualification='Equal-spectrum gray dose scaling only. The 1600 K optical response, phase enthalpy and coating gradients need a new spectrum/material calculation.'),
        insulation=dict(black_flux_W_m2=plant.SB*(1600**4-400**4),
            effective_emissivity_for6p55_W_m2=6.55/(plant.SB*(1600**4-400**4))),
        lifetime=dict(availability90_maintenance_days_year=36.525,
            specific_energy_repay_years_per_MJ_kg=1e6/(c['net_W_kg']*plant.YEAR),
            twenty_year_export_MJ_kg=c['net_W_kg']*20*plant.YEAR/1e6))
    grains=[]
    for radius in [.01,.065536]:
        mass=4*math.pi/3*1000*radius**3;qm=3*8.8541878128e-12*1e9/(1000*radius)
        grains.append(dict(radius_m=radius,grain_kg=mass,charge_C_kg=qm,
            average_grains_s_per_100k_outlet=rate/mass/1e5,
            mean_grains_s_per_outlet_in700s_large_shot=1.1309733552923254e13/700/mass/1e5,
            mean_current_A_per_outlet_in700s_large_shot=qm*1.1309733552923254e13/700/1e5,
            acceleration_length_m_at2p339Mm_s_10MV_m=2.339e6**2/(2*qm*1e7)))
    out['grain_preparation_and_launch_requirements']=grains
    out['cold_equipment_scenarios']=plant.cold_requirements()
    out['guide_buffer_dose_scenarios']=[dict(chosen_service_dose_Gy=limit,
        dose_Gy_year=43413.63018799364*16,
        reworked_polymer_kg_m2_year=20.76929422482458*43413.63018799364*16/limit,
        processing_kg_m2_at1e8_kg_per_kg_s=20.76929422482458*43413.63018799364*16/limit/plant.YEAR*1e8,
        processing_W_m2_at2MJ_kg=20.76929422482458*43413.63018799364*16/limit/plant.YEAR*2e6,
        qualification='Selected dose-to-service trigger, not a measured polyethylene lifetime. Working heat capacity, radiation chemistry, repair and outgassing require qualification.') for limit in [1e5,1e6,1e7]]
    (P/'equipment-results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
