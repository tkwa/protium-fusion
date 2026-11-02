"""Independent equal-share launcher accounting and overhead headroom.

Existing active-source hardware is partitioned, not replicated, as N changes.
New per-complex overhead is charged explicitly. A launcher can emit many packets.
"""
import hashlib
import json
import math
from pathlib import Path

HERE=Path(__file__).resolve().parent
PARENT=HERE.parent
YEAR=365.25*86400
G=6.67430e-11
SB=5.670374419e-8
TARGET=1e-3
COUNTS=[100,1000,10000]


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def positive_limit(a,b,c):
    assert a>=0 and b>0 and c<0
    return -2*c/(b+math.sqrt(b*b-4*a*c)) if a else -c/b


def main():
    cp=PARENT/'open-loop-redesign/operating-cases.json'
    ip=PARENT/'open-loop-redesign/station-inputs.json'
    mp=PARENT/'manufacturing-efficiency/sensitivity.json'
    case=json.loads(cp.read_text())['cases']['high_source_rating_surface']
    inp=json.loads(ip.read_text())
    manufacturing=json.loads(mp.read_text())
    M0=case['mass_kg'];R=case['radius_m'];A=4*math.pi*R*R
    shares=[]
    for n in COUNTS:
        mass=case['mass_breakdown_kg']
        allocated={k:mass[k]/n for k in ['peak_driver_hardware','average_driver_conditioning','two_pulse_storage']}
        shares.append({'N_complexes':n,'electric_peak_W_per_complex':case['source_peak_W']/n,
            'electric_mean_W_per_complex':case['driver_W']/n,
            'assigned_source_mass_components_kg_per_complex':allocated,
            'assigned_source_mass_sum_kg_per_complex':sum(allocated.values()),
            'stored_electric_J_per_complex':case['stored_electric_J']/n,
            'gross_flow_kg_s_per_complex':case['all_stage_flow_kg_s']/n,
            'largest_stage_electric_J_per_complex':case['stored_electric_J']/2/n,
            'stage_fuel_kg_per_complex_per_event':{k:v['mass_kg']/n for k,v in inp['cycle']['stage_pulse_data'].items()},
            'stage_electric_J_per_complex_per_event':{k:v['initial_energy_J']/.6/n for k,v in inp['cycle']['stage_pulse_data'].items()},
            'sphere_equal_share_area_m2':A/n,
            'sqrt_equal_share_area_m':math.sqrt(A/n),
            'equal_area_spherical_cap_halfangle_deg':math.degrees(math.acos(1-2/n)),
            'scope':'All N complexes share every event equally. Patch geometry is area-equivalent; neither one complex per packet nor one pencil beam per complex is assumed.'})
        assert abs(n*sum(allocated.values())-sum(mass[k] for k in allocated))<1e-12*sum(mass[k] for k in allocated)
    rows=[]
    speed=1e4
    ksupport=.001*G*speed/(8*math.pi*R**2)
    kstream=G/(2*R*speed**2)
    Mt=case['tidal_equivalent_mass_kg']
    flux=.9*SB*case['radiator_K']**4
    # Preserve baseline manufacturing machinery per unit refurbished throughput;
    # lower energy need not shrink or expand that specific throughput allocation.
    factory_mass_per_installed_kg=case['mass_breakdown_kg']['replacement_manufacturing']/M0
    for cost in [2.,3.1,10.]:
        mfg=next(x for x in manufacturing['rows'] if x['effective_total_MJ_kg']==cost)
        P0=mfg['net_W']
        rep=cost*1e6*case['replacement_fraction_per_year']/YEAR
        for standby_fraction in [0.,.01,.1]:
            standby=standby_fraction*case['driver_W']
            x=positive_limit(ksupport,TARGET+rep+ksupport*(2*M0+Mt),TARGET*M0-P0+standby)
            d_support=ksupport*x*(2*M0+Mt+x)
            d_stream=kstream*x*(2*M0+Mt+x)
            d_repair=rep*x
            # Bare additional hardware is less than total installed headroom.
            # Preserve 20% reserve, ten-year spares, factory throughput ratio,
            # support streams, and extra heat-rejection panel mass.
            additions={'reserve':.2*x,'ten_year_spares':.001*x,
                       'repair_throughput_hardware':factory_mass_per_installed_kg*x,
                       'support_streams':d_stream,
                       'extra_radiator_panels':20*(d_support+d_repair+standby)/flux}
            bare=x-sum(additions.values())
            P=P0-standby-d_support-d_repair
            assert abs(P/(M0+x)/TARGET-1)<1e-12
            assert bare>0
            rows.append({'effective_refurbishment_MJ_kg':cost,
                 'additional_standby_fraction_of_original_driver_mean':standby_fraction,
                 'additional_standby_W':standby,
                 'optimistic_extra_total_mass_kg_ignoring_added_loads':(P0-standby)/TARGET-M0,
                 'extra_total_installed_mass_at_threshold_kg':x,
                 'final_total_mass_kg':M0+x,
                 'threshold_net_W':P,
                 'additional_support_W':d_support,'additional_refurbishment_W':d_repair,
                 'bare_overhead_available_kg':bare,'necessary_mass_additions_kg':additions,
                 'maximum_bare_fixed_overhead_per_complex_kg':{str(n):bare/n for n in COUNTS},
                 'scope':'All additional unallocated station functions compete for this headroom. Added standby includes source-control/handling consumption. Factory throughput ratio is preserved from baseline, not inferred from efficiency.'})
    # Transparent fixed overhead scenarios. Totals grow with N. The mass-only
    # version is an optimistic diagnostic, not a replacement for the above bound.
    pref=next(x['net_W'] for x in manufacturing['rows'] if x['effective_total_MJ_kg']==2.)
    simple=[]
    for h in [1e17,1e18,1e19,1e20]:
        for n in COUNTS:
            simple.append({'N_complexes':n,'new_fixed_kg_per_complex':h,'new_total_kg':n*h,
                           'mass_only_net_W_kg':pref/(M0+n*h)})
    result={'scope':'Conditional 100/1000/10000-complex allocation, not component-resolved launchers or an optimization. Energy-cost sensitivities respect the >=50% premise; none assumes 100 MJ/kg is necessary.',
            'baseline_station_mass_kg':M0,'original_source_peak_W':case['source_peak_W'],
            'original_source_mean_W':case['driver_W'],
            'mean_to_peak_electric_ratio':case['driver_W']/case['source_peak_W'],
            'largest_stage_energy_over_peak_s':case['stored_electric_J']/2/case['source_peak_W'],
            'equal_share_rows':shares,'overhead_limits':rows,
            'fixed_per_complex_overhead_mass_only_scenarios':simple,
            'unallocated_overhead_equation':'M_source(N) = P_peak/r_peak + 0.1 P_mean + E_store/100000 + H_shared + N*h_fixed. Existing service/reserve may cover documented items, but the same mass cannot be assigned twice.',
            'threshold_equation':'P0 - S - r*x - k*x*(2*M0 + Mt + x) = 0.001*(M0+x), where x includes all new installed mass.',
            'count_uncertainty':'No quantitative probability distribution for N is supported. 100, 1000, 10000 are design scenarios, not 80% credible bounds. Count should follow independent spatial channels, aperture and fidelity, dose/reliability and per-complex overhead.',
            'hashes':{str(p.relative_to(PARENT)):sha(p) for p in [cp,ip,mp]},'code_sha256':sha(Path(__file__)),
            'checks':{'share_sums_match_existing_allocation':True,'threshold_power_mass_equalities_relative_1e12_inverse':True}}
    (HERE/'mass-count.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'shares':shares,'limits_2MJkg':[r for r in rows if r['effective_refurbishment_MJ_kg']==2.],'selected_fixed_overhead_scenarios':[s for s in simple if s['new_fixed_kg_per_complex']==1e19]},indent=2))

if __name__=='__main__':main()
