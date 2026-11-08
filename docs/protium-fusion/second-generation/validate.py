#!/usr/bin/env python3
"""Independent conservation, chronology, forcing, source and provenance checks."""
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import numpy as np
P=Path(__file__).resolve().parent;ROOT=P.parents[2]
checks=[]
def ck(name,value):
    checks.append(dict(name=name,passed=bool(value)))
    if not value:print('FAIL',name)
def read(path):return json.loads(Path(path).read_text(),parse_constant=lambda x:(_ for _ in ()).throw(ValueError(x)))
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def main():
    rates=read(P/'rates.json');A=np.array(rates['A']);Z=np.array(rates['Z']);mass=np.array(rates['mass_MeV']);idx={n:i for i,n in enumerate(rates['species'])}
    for r in rates['rates']:
        ck('baryon '+r['name'],sum(A[i] for i in r['react'])==sum(A[i] for i in r['product']))
        dz=sum(Z[i] for i in r['product'])-sum(Z[i] for i in r['react'])
        ck('charge/leptons '+r['name'],abs(dz)<=1 if r['weak'] else dz==0)
    for f in P.glob('*.json'):ck('strict JSON '+f.name,isinstance(read(f),dict))
    plant=read(P/'plant-results.json')
    area=4*math.pi*(1.1e9)**2
    for c in plant['cases']:
        e=c['energy_W_m2'];m=c['mass_breakdown_kg_m2'];q=c['cold_heat_W_m2']
        ck(c['name']+' mass sum',abs(sum(m.values())/c['installed_kg_m2']-1)<1e-9)
        ck(c['name']+' global energy',abs(e['nuclear_nonnu']-e['net']-sum(q.values()))<1e-7)
        ck(c['name']+' electrical export',abs(e['thermal_electric']+e['direct_after_storage']-e['driver_electric']-q['ordinary']-q['direct_auxiliary']-q['repair']-q['support']-e['net'])<1e-7)
        ck(c['name']+' radiator',abs(c['radiator_fraction']*.9*5.670374419e-8*400**4-e['waste']-11.0254)<1e-7)
        ck(c['name']+' specific power',abs(c['net_W']/c['mass_kg']/c['net_W_kg']-1)<1e-12)
        ck(c['name']+' round trip bookkeeping',abs(e['direct_gross']*.9-e['direct_after_storage'])<1e-9)
    ck('overloaded thermal20 is flagged',not next(c for c in plant['cases'] if c['name']=='thermal20')['radiator_fits'])
    for r in plant['cold_equipment_requirements']:
        ck('cold equipment component sum '+str((r['helium_active_area_multiplier'],r['guarded_200K_intercepts'])),abs(sum(r['components_W_m2'].values())-r['named_electric_W_m2'])<1e-10)
    guarded=next(r for r in plant['cold_equipment_requirements'] if r['helium_active_area_multiplier']==16 and r['guarded_200K_intercepts'])
    actual=next(c for c in plant['cases'] if c['name']=='spherical_large_N15_16')
    ck('guarded cold requirements fit central auxiliary allocation',guarded['named_electric_W_m2']<=actual['cold_heat_W_m2']['ordinary'])
    ck('unguarded cold equipment fails radiator area',not next(c for c in plant['cases'] if c['name']=='spherical16_unguarded_cold_area16')['radiator_fits'])
    for fn in ['large-n15-16x-events.csv','spherical-large-n15-16x-events.csv']:
        events=list(csv.DictReader((P/fn).open()))
        ck(fn+' serialization',all(float(b['launch_s'])>=float(a['capture_end_s'])-1e-6 for a,b in zip(events,events[1:])))
        ck(fn+' lateness',all(float(e['delay_s'])>=-1e-6 for e in events))
    ck('archived schedule reproduced',any(s['rate_multiplier']==1 and s['n15_mass_multiplier']==1 and s['event_count']==4228 and abs(s['max_delay_s']-5010.857824806124)<1e-5 for s in plant['schedules']))
    ck('old pulse 16x overload detected',any(s['rate_multiplier']==16 and s['n15_mass_multiplier']==1 and s['centres']==1 and not s['overlapping_legs'] and not s['stable_average_load'] for s in plant['schedules']))
    probes={p.name:read(p) for p in P.glob('pulse-n15-*.json')}
    for fn,r in probes.items():
        ck(fn+' complete',r['success'] and r['termination_reason']=='final_time')
        ck(fn+' pulse energy residual',abs(r['energy_residual_fraction_initial'])<3e-5)
        ck(fn+' baryons',r['max_cell_baryon_residual']<1e-6)
        ck(fn+' raw trace hash',Path(r['trace_path']).exists() and sha(r['trace_path'])==r['trace_sha256'])
    base=read(P.parent/'mixture-design/selected-plant/point.json')
    control=probes['pulse-n15-m1-v1-z128-t1.json'];big=probes['pulse-n15-m100-v1-z128-t1.json']
    ck('N15 reference control',abs(control['desired_product_per_initial_catalyst']-base['cycle']['output_catalyst_matrix_M'][0][1])<1e-7)
    ck('large pulse mass',abs(big['mass_kg']/control['mass_kg']/100-1)<1e-12)
    ck('large pulse prepared energy',abs(big['initial_energy_J']/control['initial_energy_J']/100-1)<1e-12)
    late=probes['pulse-n15-m100-v1-z128-t2.json']
    ck('large pulse nuclear end-time check',abs(late['desired_product_per_initial_catalyst']-big['desired_product_per_initial_catalyst'])<1e-10)
    ck('slower pulse failure recorded',probes['pulse-n15-m100-v0p8-z128-t1.json']['desired_product_per_initial_catalyst']<.1)
    fine=probes['pulse-n15-m100-v1-z256-t1.json']
    ck('256 relative yield change below 3pct',abs(fine['desired_product_per_initial_catalyst']/big['desired_product_per_initial_catalyst']-1)<.03)
    family=read(P.parent/'component-design/final-family-dispatch/sigma350/input.json')
    bulk=read(P/'bulk-16.json')
    for r in bulk['cases']:
        ck(r['stage']+'/'+r['case']+' abundance baryons',r['baryon_residual']<1e-10)
        if r['case']=='control':
            ck(r['stage']+' prescribed-history Q replay within 1pct',abs(r['hot_total_mass_Q_MeV_per_initial_heavy']/family['loaded_stages'][r['stage']]['Q_hot_MeV']-1)<.01)
    t32=read(P/'network-32.json');t64=read(P/'network-64.json')
    for r in t32['tracer']:
        ck(r['stage']+'/'+r['isotope']+' tagged probability',abs(sum(r['cold'].values())-1)<1e-9)
        other=next(s for s in t64['tracer'] if s['stage']==r['stage'] and s['isotope']==r['isotope'])
        if r['isotope']=='o16':
            p=1-r['cold']['o16'];q=1-other['cold']['o16']
            ck(r['stage']+' O16 radial quadrature within 3pct',abs(p/q-1)<.03)
    oxy=read(P/'oxygen-results.json');dest=oxy['final_nuclei_per_introduced_O16'];led=oxy['total_reaction_ledger']
    da=sum(A[idx[k]]*v for k,v in dest.items())-16+led['net_protons']+4*led['net_alpha']+led['net_neutrons']
    dz=sum(Z[idx[k]]*v for k,v in dest.items())-8+led['net_protons']+2*led['net_alpha']-led['net_nuclear_charge']
    dq=mass[idx['o16']]-sum(mass[idx[k]]*v for k,v in dest.items())-led['net_protons']*mass[idx['p']]-led['net_alpha']*mass[idx['he4']]-led['net_neutrons']*mass[idx['n']]-led['total_Q_MeV']
    ck('closed oxygen baryons',abs(da)<1e-7);ck('closed oxygen charge/leptons',abs(dz)<1e-7);ck('closed oxygen mass/Q',abs(dq)<1e-4)
    ck('oxygen pool restores quartet',oxy['quartet_return_probability']>.999)
    for r in read(P/'quarantine-results.json')['results']:
        ck(r['isotope']+' onward competition',r['return_to_extended_pool_probability']<1e-4)
        ck(r['isotope']+' absorption matrix residual',r['conservation_matrix_residual']<1e-10)
    for f in P.glob('*.md'):
        for link in re.findall(r'\]\(([^)]+)\)',f.read_text()):
            if not re.match(r'^[a-z]+:',link) and not link.startswith('#'):
                ck(f.name+' link '+link,(f.parent/link.split('#')[0]).exists())
    summary=dict(checks=len(checks),passed=sum(c['passed'] for c in checks),all_passed=all(c['passed'] for c in checks),
        scope=__doc__,oxygen_residuals=dict(baryons=da,charge=dz,Q_MeV=dq),
        numerical_differences=dict(N15_128_to256_conversion_relative=fine['desired_product_per_initial_catalyst']/big['desired_product_per_initial_catalyst']-1),
        checks_detail=checks,code_sha256=sha(__file__))
    (P/'validation.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k!='checks_detail'},indent=2))
    if not summary['all_passed']:raise SystemExit(1)

if __name__=='__main__':main()
