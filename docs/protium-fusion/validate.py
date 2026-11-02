"""Physical invariants and numerical checks for the screening models."""
import json
from pathlib import Path
import numpy as np
import pynucastro as p
from nuclear_screen import KB, C, RATES as CN_RATES, screen
from pulse_model import AA, ZZ, STOICH, MASS, Q, RATES, NUCS, flows

checks={}
checks['max_baryon_residual_per_reaction']=float(np.max(np.abs(AA@STOICH)))
assert checks['max_baryon_residual_per_reaction']==0
strong=[j for j,r in enumerate(RATES) if not r.weak_type]
checks['max_charge_residual_strong_reactions']=float(np.max(np.abs(ZZ@STOICH[:,strong])))
assert checks['max_charge_residual_strong_reactions']==0

yc=np.ones(len(NUCS));yc/=AA@yc
errors=[]
for t in [2e8,5e8,1e9]:
    direct=np.array([r.eval(t)*r.prefactor*1e5**r.dens_exp*
                     np.prod([yc[NUCS.index(n)] for n in r.reactants]) for r in RATES])
    errors.append(float(np.max(np.abs(flows(t,1e5,yc)-direct)/np.maximum(direct,1e-300))))
checks['max_relative_rate_error']=max(errors)
assert max(errors)<1e-9

cycle_names=['n14(p,g)o15','o15(,)n15','n15(p,a)c12','c12(p,g)n13',
             'n13(p,g)o14','o14(,)n14']
from nuclear_screen import LIB
cycle_stoich=np.zeros(len(NUCS))
for name in cycle_names:
    r=LIB.get_rate_by_name(name)
    for n in r.reactants: cycle_stoich[NUCS.index(n)]-=1
    for n in r.products: cycle_stoich[NUCS.index(n)]+=1
expected=np.zeros(len(NUCS));expected[NUCS.index(p.Nucleus('p'))]=-4
expected[NUCS.index(p.Nucleus('he4'))]=1
assert np.all(cycle_stoich==expected)
checks['closed_hot_cycle_Q_MeV']=float(-MASS@cycle_stoich)
assert abs(checks['closed_hot_cycle_Q_MeV']-26.731)<.002

# Integrating 4*pi*alpha_ff*B_nu over frequency gives this coefficient.
checks['free_free_bolometric_coefficient_SI']=float(8*np.pi*3.692e-2*KB/C**2)
checks['free_free_emission_absorption_relative_mismatch']=float(
 abs(checks['free_free_bolometric_coefficient_SI']/1.426e-40-1))
assert checks['free_free_emission_absorption_relative_mismatch']<.002

beta_wait=1/CN_RATES['bo14'].eval(1e9)+1/CN_RATES['bo15'].eval(1e9)
high=screen(3e8,1e10)
checks['hot_beta_wait_s']=beta_wait
checks['high_density_cycle_over_beta_wait']=high['cycle_s']/beta_wait
assert abs(checks['high_density_cycle_over_beta_wait']-1)<.001

eos=json.loads(Path(__file__).with_name('eos-validation-extended.json').read_text())
checks['EOS_max_relative_errors']={k:max(x[k] for x in eos) for k in ['p','e','c_v','gamma1']}
assert max(checks['EOS_max_relative_errors'].values())<2e-4

checked=[]
for path in sorted(Path(__file__).parent.glob('pulse-*-fast.json')):
    record=json.loads(path.read_text())
    assert record['integrator_success'],path
    assert abs(record['mass_fraction_residual'])<1e-8,path
    assert abs(record['energy_residual_fraction_initial'])<2e-4,path
    checked.append(path.name)
checks['pulse_files_checked']=checked
radial=[]
for path in sorted(Path(__file__).parent.glob('radial-*.json')):
    record=json.loads(path.read_text())
    assert record['success'],path
    assert abs(record['max_cell_baryon_residual'])<1e-8,path
    assert abs(record['energy_residual_fraction_initial'])<2e-4,path
    radial.append(path.name)
checks['radial_files_checked']=radial
compression=[]
for path in sorted(Path(__file__).parent.glob('compression-*.json')):
    record=json.loads(path.read_text())
    if 'nuclear_enabled' not in record: continue
    assert record['success'],path
    assert abs(record['max_cell_baryon_residual'])<1e-8,path
    assert abs(record['energy_residual_fraction_initial'])<2e-4,path
    if 'initial_kinetic_energy_J' in record:
        assert abs((record['initial_internal_energy_J']+record['initial_kinetic_energy_J'])/
                   record['initial_energy_J']-1)<1e-12,path
    compression.append(path.name)
checks['compression_files_checked']=compression
analytic=[]
for nz in [16,32,64]:
    record=json.loads(Path(__file__).with_name(f'hydro-test-{nz}.json').read_text())
    analytic.append(record['analytic_error']['temperature_mass_weighted_relative_L1'])
assert analytic[1]<analytic[0]/3 and analytic[2]<analytic[1]/3
checks['smooth_radial_analytic_T_errors']=analytic
network_diffs=[]
for nz in [64,128,256]:
    full=json.loads(Path(__file__).with_name(f'radial-A-300m-{nz}-rad.json').read_text())
    reduced=json.loads(Path(__file__).with_name(f'radial-A-300m-{nz}-reduced.json').read_text())
    network_diffs.append(abs(full['desired_product_per_initial_catalyst']-
                             reduced['desired_product_per_initial_catalyst']))
assert max(network_diffs)<3e-6
checks['reduced_A_absolute_conversion_errors']=network_diffs
coarse=json.loads(Path(__file__).with_name('radial-A-300m-512-reduced.json').read_text())
fine=json.loads(Path(__file__).with_name('radial-A-300m-1024-reduced.json').read_text())
checks['A_512_to_1024_conversion_relative_change']=abs(
    coarse['desired_product_per_initial_catalyst']/fine['desired_product_per_initial_catalyst']-1)
assert checks['A_512_to_1024_conversion_relative_change']<.001
ledger_data=json.loads(Path(__file__).with_name('cycle-ledger.json').read_text())
checks['cycle_ledger_max_baryon_residual']=max(abs(d['baryon_residual_per_hydrogen']) for d in ledger_data['cases'])
checks['cycle_ledger_max_energy_residual']=max(abs(d['mass_defect_energy_relative_residual']) for d in ledger_data['cases'])
assert checks['cycle_ledger_max_baryon_residual']<1e-10
assert checks['cycle_ledger_max_energy_residual']<1e-8
for d in ledger_data['cases']:
    components=sum(d[k] for k in ['final_bulk_kinetic_MeV_per_helium',
      'final_hot_internal_MeV_per_helium','cold_decay_heat_MeV_per_helium'])
    assert abs(components-d['available_heat_MeV_per_helium']-d['fuel_preparation_MeV_per_helium'])<1e-7
collector=json.loads(Path(__file__).with_name('collector-screen.json').read_text())
checks['collector_max_ballistic_energy_quadrature_error']=max(abs(d['energy_quadrature_relative_error']) for d in collector['cases'])
assert checks['collector_max_ballistic_energy_quadrature_error']<.002
assert abs(collector['rectangular_pulse_relative_error'])<1e-12
plant=json.loads(Path(__file__).with_name('plant-budget.json').read_text())
for d in plant['cases']:
    if not d['closes']: continue
    assert abs(sum(d['mass_budget_kg'].values())/d['total_mass_kg']-1)<1e-10
    op=d['operations']
    net=op['gross_electric_generation_W']-op['driver_electric_power_W']-op['support_loss_W']-op['other_auxiliary_reserve_W']
    assert abs(net/op['net_useful_power_W']-1)<1e-10
checks['conditional_plant_closing_cases']=sum(d['closes'] for d in plant['cases'])
checks['scope']='Conservation and numerical verification of reduced models; not reactor validation.'
Path(__file__).with_name('validation.json').write_text(json.dumps(checks,indent=2)+'\n')
print(json.dumps(checks,indent=2))
