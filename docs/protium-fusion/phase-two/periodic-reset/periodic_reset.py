"""Whole-inventory periodic reset of saved electron-poor ion processors.

The confining magnetic field remains energized. This is a prescribed inventory,
loss/capacity and maintenance screen, not a demonstrated drain/refill solution.
The main cases preserve the 95% random-ion preparation efficiency. Uniform
near-perfect efficiencies are separately marked counterfactual capabilities.
"""
from pathlib import Path
import hashlib,json,sys
import numpy as np
from scipy.integrate import simpson
from scipy.optimize import minimize_scalar,brentq
from scipy.special import erf
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT))
from electron_wall_orbits import Cylinder
from catalyst_penning import MP
from nonneutral_screen import MEV,E_CHARGE,EPS0
YEAR=365.25*86400.


def inventories(p):
    primary=p['primary_pn_rate_s'];rows=[]
    for alloc in p['channel_allocations']:
        key=str(alloc['ZA'][0]*1000+alloc['ZA'][1]);src=p['selected_source_profiles'][key]
        c=Cylinder(src);eq=c.eq;r=np.linspace(0,c.R,12001);x=r/c.scale
        nh=c.ns*np.interp(x,eq['x'],eq['h']);nc=c.ns*np.interp(x,eq['x'],eq['c'])
        z,a=src['ZA'];omega=eq['kappa']*eq['xi']*E_CHARGE*c.B/MP
        Nh=simpson(2*np.pi*r*nh,x=r);Nc=simpson(2*np.pi*r*nc,x=r)
        thermal=1.5*c.T*(Nh+Nc)
        rotation=simpson(np.pi*MP*(nh+eq['M']*nc)*omega**2*r**3/MEV,x=r)
        electric=.5*simpson(2*np.pi*r*(nh+z*nc)*c.potential_MV(r),x=r)
        field=-c.T*MEV/(E_CHARGE*c.scale)*c.spline.derivative()(np.maximum(x,eq['x'][0]))
        electric_field=simpson(np.pi*EPS0*r*field**2/MEV,x=r)
        transfer=(Nh+z*Nc)*(src['T_MeV']*p['assumptions']['end_barrier_kT'])
        total_length=alloc['pn_rate_s']/src['reaction_rate_per_m_s'];scale=total_length/primary
        row=dict(processor=key,thermal_MeV_s=thermal*scale,rotation_MeV_s=rotation*scale,
          electric_MeV_s=electric*scale,barrier_transfer_MeV_s=transfer*scale,
          ion_inventory_s=(Nh+Nc)*scale,proton_inventory_s=Nh*scale,catalyst_inventory_s=Nc*scale,
          charge_inventory_C_s=(Nh+z*Nc)*E_CHARGE*scale,
          charge_per_channel_C=(Nh+z*Nc)*E_CHARGE*p['active_length_m'],
          channel_plasma_energy_J=(thermal+rotation+electric)*MEV*p['active_length_m'],
          electric_integral_relative_difference=electric_field/electric-1,
          inventory_relative_difference=max(abs(Nh/src['proton_inventory_per_m']-1),abs(Nc/src['catalyst_inventory_per_m']-1)))
        rows.append(row)
    keys=['thermal_MeV_s','rotation_MeV_s','electric_MeV_s','barrier_transfer_MeV_s','ion_inventory_s','proton_inventory_s','catalyst_inventory_s','charge_inventory_C_s']
    sums={k:sum(row[k] for row in rows) for k in keys}
    return rows,sums


def parameters(p,total):
    a=p['assumptions'];m=p['mass_breakdown_kg'];R=p['primary_pn_rate_s']
    k=a['annual_replacement_fraction']*a['embodied_J_kg']/YEAR
    radkgW=m['radiators']/p['waste_heat_W']
    fixed=sum(v for name,v in m.items() if name not in ['radiators','replacement_manufacturing_equipment','replacement_inventory'])
    ri=a['core_radius']+a['shield'];L=p['length_m'];ro=p['outer_winding_radius_m']
    inner=a['warm_to_cold_W_m2']*(2*np.pi*ri*L+2*np.pi*ri*ri)
    outer=a['outer_to_cold_W_m2']*(2*np.pi*ro*L+2*np.pi*ro*ro)
    # The saved conductor creep is negligible at the prescribed J/Jcrit and n;
    # recover the tiny residual from the saved total cold load where possible.
    creep=max(0.,p['total_cold_load_W']-p['cold_radiation_W']-inner-outer)
    fridge_factor=(a['radiator_K']/a['coil_cold_K']-1)/a['cryo_Carnot_fraction']
    static=fridge_factor*(inner+outer+creep)+a['heat']*inner+creep+p['capture_stock_processing_W']
    return dict(R=R,fixed=fixed,radkgW=radkgW,k=k,spare=a['annual_replacement_fraction']*a['replacement_inventory_years'],
      manufacturing=a['manufacturing_kg_W'],converter_kg_W=a['conversion_kg_W'],static_preaux_W=static,
      preaux_W=p['plant_bus_net_W']/(1-a['other_fraction']),aux=a['other_fraction'],heat=a['heat'],
      nuclear_W=p['nuclear_after_capture_W'],heater=a['heater'])


def recovery_cases():
    cases=[]
    for thermal_out in [.95,.99]:
        for coherent in [.999,.9999,.99999]:
            cases.append(dict(name=f'thermal_in95_out{thermal_out}_coherent{coherent}',thermal_in=.95,thermal_out=thermal_out,
              coherent_in=coherent,coherent_out=coherent,rotation_coherent=True,category='retains_95pct_random_ion_preparation'))
    cases.append(dict(name='reference_rotation_also95_coherent999',thermal_in=.95,thermal_out=.95,
              coherent_in=.999,coherent_out=.999,rotation_coherent=False,category='retains_95pct_random_and_rotation_preparation'))
    for eta in [.95,.99,.999,.9999]:
        cases.append(dict(name=f'uniform_{eta}',thermal_in=eta,thermal_out=eta,
              coherent_in=eta,coherent_out=eta,rotation_coherent=True,category='uniform_recovery_capability_counterfactual'))
    return cases


def reset_coefficients(total,case,barrier,heat):
    # rho phi/2 is the stored self-energy. Using rho phi would double-count it.
    thermal=total['thermal_MeV_s'];coherent=total['electric_MeV_s']
    if case['rotation_coherent']:coherent+=total['rotation_MeV_s']
    else:thermal+=total['rotation_MeV_s']
    if barrier:coherent+=total['barrier_transfer_MeV_s']
    loss=lambda ei,eo:1/ei-eo-heat*(1-eo)
    A=thermal*loss(case['thermal_in'],case['thermal_out'])+coherent*loss(case['coherent_in'],case['coherent_out'])
    # Add extraction input energy and reinjection electrical energy. Two legs
    # are processed in separate parts of a channel's reset interval.
    throughput=thermal*(1+1/case['thermal_in'])+coherent*(1+1/case['coherent_in'])
    return A,throughput


def ledger(tau,delta,para,total,A,throughput,B,layout,pair_supply):
    period=tau+delta;duty=tau/period;conversion=para['R']*MEV
    reset_W=conversion*A/period
    accumulation_W=conversion*B*tau*tau/(2*period)
    # A small stated all-in per-prompt-electron allowance, kept separate from
    # the trapped-positron population; it is not derived from the escape model.
    electron_W=conversion*.1*pair_supply*duty
    loss=reset_W+accumulation_W+electron_W
    preaux=duty*(para['preaux_W']+para['static_preaux_W'])-para['static_preaux_W']-loss
    bus=preaux-para['aux']*max(0.,preaux)
    peak=conversion*throughput/(period if layout=='shared_staggered' else delta)
    converter=peak*para['converter_kg_W']
    nuclear=duty*para['nuclear_W']
    denom=1-para['k']*(para['radkgW']+para['manufacturing'])-para['spare']
    mass=(para['fixed']+converter+(nuclear-bus)*para['radkgW'])/denom
    maintenance=para['k']*mass;net=bus-maintenance;waste=nuclear-net
    return dict(burn_s=float(tau),reset_downtime_s=delta,duty=float(duty),reset_loss_MeV_per_active_primary=A/tau,
      positron_loss_MeV_per_active_primary=B*tau/2,electron_allowance_MeV_per_active_primary=.1*pair_supply,
      average_reset_loss_W=reset_W,average_positron_loss_W=accumulation_W,
      extra_converter_rated_W=peak,extra_converter_mass_kg=converter,total_mass_kg=mass,
      plant_bus_W=bus,replacement_W=maintenance,net_W=net,specific_power_mW_kg=net/mass*1000,
      nuclear_heat_W=nuclear,waste_heat_W=waste,
      end_of_burn_positrons=pair_supply*para['R']*tau,
      end_of_burn_positron_charge_fraction=pair_supply*tau/(total['charge_inventory_C_s']/E_CHARGE),
      average_inventory_reset_current_A=total['charge_inventory_C_s']*para['R']/period,
      synchronous_dump_average_current_A=2*total['charge_inventory_C_s']*para['R']/delta if delta>0 else None,
      energy_closure_W=net+waste-nuclear)


def optimum(delta,para,total,A,throughput,B,layout,pair_supply):
    def objective(logtau):return -ledger(np.exp(logtau),delta,para,total,A,throughput,B,layout,pair_supply)['specific_power_mW_kg']
    grid=np.linspace(np.log(1e-3),np.log(1e7),150);values=[objective(x) for x in grid];i=int(np.argmin(values))
    lo=grid[max(0,i-1)];hi=grid[min(len(grid)-1,i+1)]
    sol=minimize_scalar(objective,bounds=(lo,hi),method='bounded',options={'xatol':1e-10})
    out=ledger(np.exp(sol.x),delta,para,total,A,throughput,B,layout,pair_supply)
    out['optimization_success']=bool(sol.success)
    out['loss_only_optimum_burn_s']=float(np.sqrt(2*A/B))
    out['loss_only_minimum_MeV_per_active_primary']=float(np.sqrt(2*A*B))
    return out


def free_stream_tail(L,vth,t):
    x=L/(np.sqrt(2)*vth*t)
    return float(erf(x)-(-np.expm1(-x*x))/(np.sqrt(np.pi)*x))


def main():
    source=ROOT/'network-plant.json';allp=json.loads(source.read_text())
    p=next(p for p in allp['cases'] if p['assumptions']['B']==140 and p['assumptions']['heat']==.7 and p['assumptions']['end_in']==.999)
    rows,total=inventories(p);para=parameters(p,total);cases=recovery_cases();scans=[]
    supplied=[.000253,.000906];powers=[10.,30.,100.]
    cost_models={'radiated_ion_energy':1/para['heater']-para['heat'],'already_electrical':1.}
    for case in cases:
        for barrier in [False,True]:
            A,throughput=reset_coefficients(total,case,barrier,para['heat'])
            for b in supplied:
                for power in powers:
                    for cost_model,factor in cost_models.items():
                        for delta in [.03,.1,1.,10.]:
                            for layout in ['shared_staggered','synchronous']:
                                out=optimum(delta,para,total,A,throughput,b*power*factor,layout,b)
                                scans.append(dict(recovery=case['name'],category=case['category'],barrier_transfer=barrier,
                                   pair_supply_per_primary=b,positron_power_MeV_s=power,cost_model=cost_model,
                                   electrical_cost_factor=factor,converter_layout=layout,A_MeV_s=A,
                                   reset_two_leg_throughput_MeV_s=throughput,**out))
    thresholds=[]
    for case in cases:
        A,throughput=reset_coefficients(total,case,False,para['heat'])
        for b in supplied:
            def f(logp):
                return optimum(.1,para,total,A,throughput,b*np.exp(logp)*cost_models['radiated_ion_energy'],'shared_staggered',b)['specific_power_mW_kg']-1.
            try:
                lim=np.exp(brentq(f,np.log(1e-7),np.log(1e5),xtol=1e-9))
            except ValueError:lim=None
            thresholds.append(dict(recovery=case['name'],pair_supply_per_primary=b,
                maximum_radiated_positron_power_MeV_s_for_1mWkg=lim,reset_downtime_s=.1,converter_layout='shared_staggered',barrier_transfer=False))
    # Reproduce the original plant with no reset converter and no cycling loss.
    baseline=ledger(1.,0.,para,total,0.,0.,0.,'shared_staggered',0.)
    L=p['active_length_m'];species=[]
    for key,m in [('H',MP),('Ti49',MP*Cylinder(p['selected_source_profiles']['22049']).eq['M'])]:
        s=np.sqrt(1.6*MEV/m)
        species.append(dict(species=key,thermal_parallel_sigma_m_s=s,characteristic_transit_s=L/s,
          residual_fraction_at_time_s={str(t):free_stream_tail(L,s,t) for t in [.03,.1,1.,10.]}))
    short_burn=[]
    ti_lifetime=rows[0]['catalyst_inventory_s']
    for case in cases:
        A,throughput=reset_coefficients(total,case,False,para['heat'])
        for tau in [1.,10.,30.,100.]:
            entry=ledger(tau,.1,para,total,A,throughput,supplied[0]*10.*cost_models['radiated_ion_energy'],'shared_staggered',supplied[0])
            short_burn.append(dict(recovery=case['name'],pair_supply_per_primary=supplied[0],positron_power_MeV_s=10.,
              Ti_depletion_fraction_if_no_replacement_and_fixed_proton_profile=float(-np.expm1(-tau/ti_lifetime)),
              constant_profile_power_is_optimistic_without_selective_processing=True,**entry))
    checks={
       'baseline_net_relative_error':abs(baseline['net_W']/p['net_W']-1),
       'baseline_mass_relative_error':abs(baseline['total_mass_kg']/p['total_mass_kg']-1),
       'max_electric_integral_relative_difference':max(abs(r['electric_integral_relative_difference']) for r in rows),
       'max_inventory_integral_relative_difference':max(abs(r['inventory_relative_difference']) for r in rows),
       'all_optimizations_successful':all(row['optimization_success'] for row in scans),
       'largest_energy_closure_relative_error':max(abs(row['energy_closure_W'])/max(1.,row['nuclear_heat_W']) for row in scans),
       'saved_thermal_matches_three_halves_NkT_relative_error':abs(total['thermal_MeV_s']/(1.5*1.6*total['ion_inventory_s'])-1),
    }
    assert checks['max_electric_integral_relative_difference']<1e-6
    assert checks['max_inventory_integral_relative_difference']<1e-6
    assert checks['baseline_net_relative_error']<1e-10
    assert checks['baseline_mass_relative_error']<1e-10
    assert checks['all_optimizations_successful']
    assert checks['largest_energy_closure_relative_error']<1e-9
    output=dict(scope=__doc__,input_sha256={path.name:hashlib.sha256(path.read_bytes()).hexdigest() for path in [source,Path(__file__),ROOT/'electron_wall_orbits.py',ROOT/'catalyst_penning.py']},
      baseline={k:p[k] for k in ['net_W','total_mass_kg','net_W_kg','primary_pn_rate_s','stored_field_energy_J','active_length_m']},
      inventory_per_primary_rate=total,inventory_energy_J={k:total[k]*para['R']*MEV for k in total if k.endswith('MeV_s')},processor_inventories=rows,
      baseline_reproduction_relative_error={k:baseline[k]/p[k]-1 for k in ['net_W','total_mass_kg']},
      static_preauxiliary_power_W=para['static_preaux_W'],recovery_cases=cases,
      free_streaming_tail_diagnostic=species,thresholds=thresholds,cases=scans,short_burn_sensitivity=short_burn,checks=checks,
      qualifications=[
       'Uniform magnetic energy is excluded because the coil remains energized; gate capacitance and time-dependent plasma self-fields need a switching design.',
       'Lowered-gate refill avoids assigning 40MV times total ion charge automatically. Barrier-transfer cases show the added circulation if that voltage is traversed.',
       'All hot ion inventory is processed each cycle. Selective V extraction and catalyst recycling during burn remain baseline assumptions; resetting alone does not prove their control.',
       'Shared staggered conversion assumes channels can route their high-power discharge/refill to common converters. Local switch, cable and pulse-current capability are not proved. Synchronous cases rate all channels at once.',
       'The main cases retain 95% random-ion heating and 95% or99% random-energy extraction; treating rotation as coherent is optimistic. Uniform high efficiencies would require additional random-energy recovery or adiabatic-expansion apparatus and entropy export.',
       'Positrons accumulate linearly from the stated prompt-pair supply until a complete reset; no assumed steady abundance or annihilation equilibrium is used. Constant10--100MeV/s per trapped positron is a diagnostic parameter, not a measured rate.',
       'The radiated-ion-energy cost factor is1/0.95−0.7, matching heater replacement minus thermal-conversion credit. Already-electrical scenarios take factor1.',
       'Prompt electrons are given a separate assumed0.1MeV electrical cost each. Their earlier4.7ms escape and77.8keV radiation benchmark does not establish the positron lifetime or full gas-load cost.',
       'Free-streaming residuals are diagnostics, not bounds with evolving axial fields. A specified low residual and post-dump axial sweep need demonstration.',
       'Extra converter mass uses the existing0.1kg/Wgross assumption and is added before re-closing radiator, replacement factory and spare inventory. Unspecified buffer storage and drain/refill hardware would add mass.',
       'No reactor success, future construction date or empirical performance is established by these conditional allocations.'
      ])
    (HERE/'periodic-reset.json').write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    print('baseline',output['baseline_reproduction_relative_error'])
    for case in cases:
        chosen=[x for x in scans if x['recovery']==case['name'] and not x['barrier_transfer'] and x['reset_downtime_s']==.1 and x['converter_layout']=='shared_staggered' and x['cost_model']=='radiated_ion_energy']
        for x in [chosen[0],chosen[-1]]:
            print(case['name'],x['pair_supply_per_primary'],x['positron_power_MeV_s'],'burn',x['burn_s'],'mWkg',x['specific_power_mW_kg'],'netTW',x['net_W']/1e12,'convPW',x['extra_converter_rated_W']/1e15)
    print('thresholds',json.dumps(thresholds,indent=2))
if __name__=='__main__':main()
