"""Independent unit, conservation, quadrature, and resolution checks."""
from pathlib import Path
import json
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
from delivery_screen import EPS0,SIGMA,charge_accelerator,survival

HERE=Path(__file__).resolve().parent
x=json.loads((HERE/'delivery-screen.json').read_text())
checks=[]
def check(name, measured, tolerance, detail=''):
    checks.append({'name':name,'error_or_value':float(measured),'maximum':float(tolerance),
                   'passed':bool(measured<=tolerance),'detail':detail})

# Charge from sphere capacitance and surface potential, not charge/mass formula.
a,rho,Es,E,L=.017,1700.,3e8,2e6,4e8
m=4*np.pi*a**3*rho/3
C=4*np.pi*EPS0*a
q=C*(Es*a)
v=np.sqrt(2*q*E*L/m)
s=charge_accelerator(v,rho,a,Es,E,L)
check('sphere capacitance -> accelerator work',abs(s['minimum_length_m']/L-1),2e-14)
check('charge-mass scaling',abs(s['maximum_charge_per_mass_C_kg']/(q/m)-1),2e-14)

# Integrate an isotropic blackbody radiance over a receiving hemisphere.
T=734.; radiance=SIGMA*T**4/np.pi
hemisphere=2*np.pi*quad(lambda theta:radiance*np.cos(theta)*np.sin(theta),0,np.pi/2)[0]
emit_absorb=.02;dt=4321.;u=1e6
P=4*np.pi*a*a*hemisphere*emit_absorb
speed=2e10/dt
ss=survival(speed,rho,a,T,emit_absorb,u)
check('isotropic cavity hemispheric radiance integral',abs(ss['radiative_absorbed_energy_J_kg']/(P*dt/m)-1),2e-14)

# Sweep intercepted gas volume, integrate kinetic energy, compare thermal budget.
column=ss['maximum_gas_column_kg_m2_at_full_energy_deposition']
gas_mass=column*np.pi*a*a
check('stationary gas interception work',abs(.5*gas_mass*speed**2/(m*u)-1),2e-14)

# Exact uniform homologous solution under ballistic backward mapping.
M=1e16;R=1000.;H=2e3;lookback=.004
r=np.linspace(0,R,65);shell_mass=M*np.diff((r/R)**3)
rr=r*(1+H*lookback)
rhomap=shell_mass/(4*np.pi/3*np.diff(rr**3))
analytic=(3*M/(4*np.pi*R**3))/(1+H*lookback)**3
check('homologous ballistic density',np.max(abs(rhomap/analytic-1)),2e-13)

# Independent uniform-sphere optical-depth integral: average chord times n pi a^2.
n=M/((4*np.pi/3*R**3)*m)
mean_tau=quad(lambda b:2*np.sqrt(R*R-b*b)*n*np.pi*a*a*2*b/R**2,0,R,epsabs=1e-12)[0]
projected_total=(M/m)*np.pi*a*a/(np.pi*R*R)
check('grain optical depth area / chord integral',abs(mean_tau/projected_total-1),2e-12)

for c in x['cases']:
    label=c['stage']+str(c['zones'])
    check(label+' node kinetic-energy reconstruction',abs(c['kinetic_ledger_relative_error']),5e-14)
    for p in c['contact_backmaps']:
        check(label+' first-contact root '+str(p['material_density_kg_m3']),
              abs(p['maximum_cell_density_kg_m3']/p['material_density_kg_m3']-1),1e-9)
    for p in c['cold_pellet_acceleration_survival_trade']:
        aa=p['maximum_accelerable_pellet_radius_m']
        ss=survival(c['energy_rms_speed_m_s'],p['pellet_density_kg_m3'],aa,
                    p['wall_temperature_K'],p['maximum_Planck_absorptivity_for_cold_flight'],
                    p['allowable_absorbed_heat_J_kg'])
        check(label+' coupled thermal acceleration equality '+str(p['wall_temperature_K'])+'/'+str(p['accelerator_length_m']),
              abs(ss['radiative_absorbed_energy_J_kg']/p['allowable_absorbed_heat_J_kg']-1),2e-14)
    for p in c['uncorrected_pointing_requirements']:
        t=p['distance_after_last_correction_m']/c['energy_rms_speed_m_s']
        # Cartesian propagation with a perturbation in transverse velocity.
        offset=p['uncorrelated_transverse_velocity_error_m_s']*t
        check(label+' perturbed Cartesian trajectory '+str(p['footprint_tolerance_m'])+'/'+str(p['distance_after_last_correction_m']),
              abs(offset/p['footprint_tolerance_m']-1),2e-14)
    for p in c['granular_preheat_requirements']:
        mass=c['mass_kg'];radius=p['pellet_radius_m'];density=p['pellet_material_density_kg_m3']
        grain_mass=4*np.pi/3*radius**3*density
        area=(mass/grain_mass)*np.pi*radius**2
        Rtau=p['mean_projected_grain_optical_depth_one_radius_m']
        check(label+' sum of individual grain cross-sections '+str(radius),abs(area/(np.pi*Rtau**2)-1),2e-14)

coarse,fine=x['cases'][:2]
for key,tol in [('energy_rms_speed_m_s',.001)]:
    check('A128/A512 '+key,abs(coarse[key]/fine[key]-1),tol,'Input-profile resolution comparison, not delivery-engineering accuracy.')
for quantile in ('0.05','0.5','0.95'):
    a=coarse['launch_schedule']['mass_weighted_flight_time_quantiles_s'][quantile]
    b=fine['launch_schedule']['mass_weighted_flight_time_quantiles_s'][quantile]
    check('A128/A512 flight quantile '+quantile,abs(a/b-1),.02)
pa=coarse['launch_schedule']['peak_kinetic_source_power_W'];pb=fine['launch_schedule']['peak_kinetic_source_power_W']
check('A128/A512 peak launch power',abs(pa/pb-1),.02,'Finite-volume launch-time bins; narrow time structure is unresolved.')

result={'checks':checks,'passed':sum(c['passed'] for c in checks),'total':len(checks),
        'scope':'Verifies formulas and saved-input accounting. Does not validate cold EOS, absorber properties, focusing, or accelerator performance.'}
(HERE/'delivery-validation.json').write_text(json.dumps(result,indent=2)+'\n')
print(result['passed'],'/',result['total'])
for c in checks:
    if not c['passed']:print('FAIL',c)
if result['passed']!=result['total']:raise SystemExit(1)
