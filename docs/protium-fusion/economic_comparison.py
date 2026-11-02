"""Energy-capital and fuel alternatives for the conditional protium plants.

Engineering efficiencies, solar areal masses, and economic discount rates are
explicit scenarios. This is an energy/mass comparison, not a price forecast.
The opaque allocation is extended to include replacement manufacture and stock.
"""
from pathlib import Path
import argparse,json,math
import numpy as np
from scipy.optimize import brentq
from nonneutral_screen import Q,MEV,C,MU,mass

YEAR=365.25*86400;KB=1.380649e-23;G=6.67430e-11
SOLAR_FLUX=1361.;SOLAR_L=3.828e26;GMJ=1.2668653e17;RJ=7.1492e7


def opaque_replacement(base,fraction=1e-4,embodied=1e8,factory_kg_W=1.,spare_years=10.):
    """Reclose the existing opaque requirements allocation with replacement.

    The same pulse, collector, driver, recovery and power-component assumptions
    are retained. Higher throughput must still be physically qualified.
    """
    p=base['parameters'];R=p['collector_radius_m'];v=p['support_stream_speed_m_s']
    fixed_names=['absorber','insulation','processing_and_inventory_allowance',
        'support_guides_controls_and_redundancy_allowance','driver_facility_allowance','pulse_electrical_energy_storage']
    fixed={k:base['mass_budget_kg'][k] for k in fixed_names};M0=sum(fixed.values())
    q=base['cycle_heat_MeV_per_helium']*MEV;u=base['prepared_energy_MeV_per_helium']*MEV
    generated=p['thermal_conversion_efficiency']*(q+u);drive=u/p['driver_efficiency'];cycle_net=generated-drive
    radiator_flux=p['radiator_emissivity']*5.670374419e-8*p['radiator_T_K']**4
    h_per_he=base['operations']['hydrogen_consumption_kg_s']/MU/base['operations']['helium_atoms_per_s']
    def state(M):
        stored=G*M*M/(4*R);support=p['support_fractional_energy_loss_per_orbit']*stored/(2*np.pi*R/v)
        useful=.001*M;replacement=fraction*M*embodied/YEAR
        rate=(useful+support+replacement)/((1-p['other_auxiliary_fraction_of_preaux_net'])*cycle_net)
        gross=generated*rate;driver=drive*rate;aux=p['other_auxiliary_fraction_of_preaux_net']*cycle_net*rate
        heat=q*rate-useful
        variable={'circulating_support_streams':2*stored/v**2,
            'thermal_power_conversion':p['gross_generator_kg_W']*gross,
            'driver_average_power_conditioning':p['driver_conditioning_kg_W']*driver,
            'radiators':p['radiator_kg_m2']*heat/radiator_flux,
            'replacement_manufacturing_equipment':factory_kg_W*replacement,
            'replacement_inventory':fraction*M*spare_years}
        return sum(fixed.values())+sum(variable.values()),{
            'total_mass_kg':M,'net_W':useful,'net_W_kg':.001,'helium_atoms_s':rate,
            'gross_electric_W':gross,'driver_electric_W':driver,'support_loss_W':support,'other_auxiliary_W':aux,
            'replacement_manufacturing_W':replacement,'waste_heat_W':heat,'nuclear_captured_W':q*rate,
            'support_stored_energy_J':stored,'hydrogen_atoms_s':rate*h_per_he,
            'radiator_area_m2':heat/radiator_flux,'radiator_area_over_collector':heat/radiator_flux/(4*np.pi*R**2),
            'mass_budget_kg':{**fixed,**variable},'electric_balance_W':gross-driver-support-aux-replacement-useful,
            'thermal_balance_W':q*rate-useful-heat}
    # Solve in units of the old mass, avoiding an absolute kg tolerance.
    scale=base['total_mass_kg'];grid=np.geomspace(M0/scale,1.89813e27/scale,1000)
    brackets=[(a,b) for a,b in zip(grid,grid[1:]) if (state(a*scale)[0]/scale-a)*(state(b*scale)[0]/scale-b)<0]
    if not brackets:return {'closes':False,'replacement_fraction_year':fraction,'embodied_J_kg':embodied}
    x=brentq(lambda x:state(x*scale)[0]/scale-x,*brackets[0],xtol=1e-13)
    total,d=state(x*scale)
    # An independent quadratic root uses the fact that support is quadratic
    # in M and every other newly scaled component is linear in M or support.
    y0=state(0)[0]/scale;y1=state(scale)[0]/scale-1;y2=state(2*scale)[0]/scale-2
    aa=(y2-2*y1+y0)/2;bb=y1-y0-aa;cc=y0
    analytic=2*cc/(-bb+math.sqrt(bb*bb-4*aa*cc))
    d.update(closes=True,replacement_fraction_year=fraction,embodied_J_kg=embodied,
        factory_kg_W=factory_kg_W,spare_years=spare_years,mass_relative_residual=total/d['total_mass_kg']-1,
        independent_quadratic_mass_relative_error=analytic/x-1,
        mass_fixed_point_derivative=2*aa*x+bb+1,
        throughput_over_original=d['helium_atoms_s']/base['operations']['helium_atoms_per_s'],
        original_mean_pulse_interval_s=base['operations']['mean_interval_between_all_pulses_s'],
        revised_mean_pulse_interval_s=base['operations']['mean_interval_between_all_pulses_s']*base['operations']['helium_atoms_per_s']/d['helium_atoms_s'])
    d['radiator_area_exceeds_collector_area']=d['radiator_area_over_collector']>1
    d['allocation_warnings']=(['Radiator area exceeds the entire collector surface; additional geometry and support must be modeled.']
        if d['radiator_area_exceeds_collector_area'] else
        ['Radiator area exceeds half the collector surface; detailed view factors and placement are required.']
        if d['radiator_area_over_collector']>.5 else [])
    return d


def discounted_payback(s,capital_J_kg,rate_year):
    tau=capital_J_kg/(s*YEAR)
    if rate_year==0:return tau
    if rate_year*tau>=1:return None
    return -math.log1p(-rate_year*tau)/rate_year


def run(root):
    plants=json.loads((root/'network-plant.json').read_text())['cases']
    p=next(x for x in plants if x['assumptions']['B']==140 and x['assumptions']['heat']==.7 and x['assumptions']['end_in']==.999)
    network=json.loads((root/'network-ion-budget.json').read_text())['network']
    opaque=json.loads((root/'plant-budget.json').read_text())['cases'][0]
    mH=mass('p')*MEV/C**2;mD=mass('d')*MEV/C**2;mHe=mass('he4')*MEV/C**2
    H_rate=p['primary_pn_rate_s']*network['net_H_per_primary_pn'];Hgas=H_rate*mH
    hwork=p['net_W']/Hgas
    # Catalyzed DD plus burning the T and He3 products: 6D -> 2He+2p+2n.
    # Neutron recapture on the two product protons closes 4D -> 2He.
    qdd=Q['DD_Tp']+Q['DD_He3n']+Q['DT']+Q['DHe3']
    rawD=qdd*MEV/(6*mD);closedD=2*Q['net_D']*MEV/(4*mD)
    variants=[]
    for f in [0.,1e-5,1e-4,3e-4,1e-3]:variants.append(opaque_replacement(opaque,f))
    s=p['net_W_kg'];capital=[]
    for e in [1e6,1e7,1e8,1e9]:
        for charge_field in [False,True]:
            k=e+(p['stored_field_energy_J']/p['total_mass_kg']/.95 if charge_field else 0.)
            capital.append({'construction_J_kg':e,'field_charge_included':charge_field,'startup_charge_efficiency':.95,
                'total_capital_energy_J_kg':k,'simple_payback_year':k/(s*YEAR),
                'maximum_continuous_energy_growth_rate_year':s*YEAR/k,
                'ideal_continuous_energy_doubling_year':math.log(2)*k/(s*YEAR),
                'discounted_payback':[{'rate_year':r,'payback_year':discounted_payback(s,k,r)} for r in [0.,1e-4,3e-4,5e-4,.001,.01]]})
    solar=[]
    for eta in [.1,.3,.6]:
        for sigma in [.1,1.,10.,100.]:
            solar.append({'net_solar_efficiency':eta,'complete_effective_areal_mass_kg_m2':sigma,
                'mass_crossover_AU_for_1mW_kg':math.sqrt(eta*SOLAR_FLUX/(sigma*.001)),
                'mass_crossover_AU_for_nonneutral_reference':math.sqrt(eta*SOLAR_FLUX/(sigma*s)),
                'net_W_kg_at_100AU':eta*SOLAR_FLUX/(sigma*1e4)})
    deuterium=[]
    for ratio in [1.9e-5,2.6e-5,3.3e-5,1e-7]:
        x=ratio/(1+ratio);average=(1-x)*mH+x*mD;yf=x*mD/average
        entropy=-x*math.log(x)-(1-x)*math.log1p(-x)
        entropy_work=KB*300*entropy/average
        for eta in [.3,.6]:
            for cycle,qkg in [('DD_secondary_burn',rawD),('DD_with_neutron_recapture',closedD)]:
                workD=eta*qkg;feedwork=yf*workD
                deuterium.append({'D_over_H_number_ratio':ratio,'D_atom_fraction':x,'D_mass_fraction_in_hydrogen':yf,
                    'net_D_energy_conversion_efficiency':eta,'D_cycle':cycle,'nuclear_J_per_kg_D':qkg,'net_J_per_kg_D':workD,
                    'net_J_per_kg_hydrogen_feed':feedwork,'ideal_atomic_mixture_separation_J_per_kg_feed':entropy_work,
                    'ideal_atomic_mixture_separation_J_per_kg_D':entropy_work/yf,
                    'ideal_separation_fraction_of_D_work':entropy_work/feedwork,
                    'raw_hydrogen_feed_kg_s_for_reference_power':p['net_W']/feedwork,
                    'D_kg_s_for_reference_power':p['net_W']/workD,
                    'protium_work_over_D_work_in_same_feed':hwork/feedwork,
                    'lift_all_feed_from_Jupiter_energy_fraction':GMJ/RJ/feedwork,
                    'lift_only_D_from_Jupiter_energy_fraction':GMJ/RJ/workD,
                    'separator_mass_capacity_limit_kg_per_kg_s_if_D_reactor_1W_kg':(1/s-1)*feedwork,
                    'years_from_1e27kg_hydrogen_at_reference_power':1e27*feedwork/p['net_W']/YEAR,
                    'years_from_1e27kg_hydrogen_at_solar_luminosity':1e27*feedwork/SOLAR_L/YEAR})
    accretion=[]
    for efficiency in [.0570,.302]:
        for recovery in [.3,.5]:
            work=efficiency*recovery*C*C
            accretion.append({'thin_disk_outgoing_radiation_fraction_of_rest_mass':efficiency,
                'assumed_capture_conversion_efficiency':recovery,'useful_J_per_kg_feed':work,
                'fuel_work_over_nonneutral_protium_reference':work/hwork,
                'feed_kg_s_for_reference_power':p['net_W']/work})
    return {'scope':__doc__,'constants':{'seconds_per_year':YEAR,'solar_irradiance_W_m2':SOLAR_FLUX,'solar_luminosity_W':SOLAR_L,'Jupiter_GM_m3_s2':GMJ,'Jupiter_equatorial_radius_m':RJ,'hydrogen_atomic_mass_kg':mH,'deuterium_atomic_mass_kg':mD},
        'nonneutral_reference':{'net_W':p['net_W'],'bus_W':p['plant_bus_net_W'],'mass_kg':p['total_mass_kg'],'net_W_kg':s,
            'neutral_hydrogen_feed_kg_s':Hgas,'helium_kg_s':p['primary_pn_rate_s']*network['net_alpha_per_primary_pn']*mHe,
            'net_J_per_kg_neutral_hydrogen':hwork,'useful_rest_mass_fraction':hwork/C**2,
            'stored_field_J':p['stored_field_energy_J'],'field_charge_output_year_at_95percent':p['stored_field_energy_J']/.95/p['net_W']/YEAR,
            'replacement_W':p['replacement_manufacturing_W'],'replacement_fraction_of_bus':p['replacement_manufacturing_W']/p['plant_bus_net_W'],
            'replacement_fraction_year':p['assumptions']['annual_replacement_fraction'],'embodied_J_kg':p['assumptions']['embodied_J_kg'],
            'already_includes_replacement':True},
        'opaque_replacement_cases':variants,'capital_scenarios':capital,'solar_scenarios':solar,'deuterium_scenarios':deuterium,'accretion_scenarios':accretion,
        'nuclear_inputs':{'atomic_masses_MeV':{n:mass(n) for n in ['p','d','t','he3','he4','n']},'reaction_Q_MeV':{k:Q[k] for k in ['DD_Tp','DD_He3n','DT','DHe3','np_capture','net_D']},'library':'pynucastro 3.1.0; same retained mass data as the nuclear ledgers'},
        'DD_nuclear_ledger':{'six_D_secondary_burn_Q_MeV':qdd,'four_D_closed_Q_MeV':2*Q['net_D'],'two_np_capture_Q_MeV':2*Q['np_capture'],'closure_residual_MeV':qdd+2*Q['np_capture']-2*Q['net_D']},
        'limits':['Whole-plant isotope separation and conventional D reactor equipment are not designed here.',
            'Atomic ideal-mixture entropy at 300 K is a reversible benchmark; molecular chemistry, throughput, pumping, cooling and finite recovery require separate costs.',
            'Solar areal mass includes all hardware assigned to local useful electricity; transport and remote transmission must be added when needed.',
            'Thin-disk energy fractions are published model results, not a demonstrated black-hole power plant. Seed mass, acquisition, feeding and heat rejection remain explicit requirements.',
            'Opaque replacement cases preserve unproved pulse-driver, collection and structural assumptions; their increased cadence needs engineering qualification.',
            'Discount and embodied-energy values are scenarios, not forecasts or market prices.']}


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=Path(__file__).parent);ap.add_argument('--out',type=Path);a=ap.parse_args();d=run(a.root)
    (a.out or a.root/'economic-comparison.json').write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
    print(json.dumps(d['nonneutral_reference'],indent=2));print('opaque replacement',[(x['replacement_fraction_year'],x.get('total_mass_kg'),x.get('throughput_over_original')) for x in d['opaque_replacement_cases']])
