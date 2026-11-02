"""Independent analytic and numerical checks of the shared-field allocation."""
from bundle_solenoid import *
import gzip
from math import exp
from scipy.special import erfi
from scipy.integrate import quad
from scipy.optimize import brentq
from ion_radiation_and_ends import orbit_integral
checks=[]
def check(name,value,tolerance):
 checks.append({'name':name,'error':float(abs(value)),'tolerance':tolerance,'passed':bool(abs(value)<tolerance)})

for B,a,J,s in [(40,10,1e8,2e9),(100,105,1e8,8e9),(180,1005,1e8,16e9)]:
 w=winding(B,a,J,s)
 scale=sqrt(2*MU0*s)
 def r(b):
  return exp(-b*b/(2*MU0*s))*(a*exp(B*B/(2*MU0*s))+sqrt(pi*MU0*s/2)/(MU0*J)*(erfi(B/scale)-erfi(b/scale)))
 ro=r(0)
 mc=quad(lambda b:2*pi*r(b)*8000/(MU0*J),0,B,epsabs=0,epsrel=1e-10)[0]
 mh=quad(lambda b:2*pi*r(b)**2*1600*b/(MU0*s),0,B,epsabs=0,epsrel=1e-10)[0]
 ue=B*B*pi*a*a/(2*MU0)+quad(lambda b:pi*r(b)*b*b/MU0*(1/(MU0*J)+b*r(b)/(MU0*s)),0,B,epsabs=0,epsrel=1e-10)[0]
 for k,v in [('outer_radius_m',ro),('conductor_mass_kg_per_m',mc),('hoop_mass_kg_per_m',mh),('field_energy_J_per_m',ue)]:
  check(f'analytic_radius_quadrature_B{B}_{k}',w[k]/v-1,2e-9)
 check(f'hoop_virial_B{B}',w['hoop_virial_relative_error'],1e-8)
 wi=winding(B,a,1e22,s)
 limit=pi*s*a*a/2*np.expm1(B*B/(MU0*s))
 check(f'infinite_J_energy_B{B}',wi['field_energy_J_per_m']/limit-1,1e-9)
 # A single-radius winding has an exact elementary on-axis field.
 wt={'_r':np.full(1001,a),'_B':np.linspace(B,0,1001),'B_T':B}
 L=20*a
 check(f'on_axis_thin_sheet_B{B}',finite_axis_factor(wt,L,0)/(L/sqrt(L*L+4*a*a))-1,1e-12)

rows=json.loads(gzip.decompress((ROOT/'bundle-scan.json.gz').read_bytes()))['cases']
check('all_cases_mass_ledger',max(abs(r['mass_closure_kg'])/max(1,r['total_mass_kg']) for r in rows),1e-12)
check('all_cases_energy_ledger',max(abs(r['energy_closure_W'])/max(1,abs(r['net_W']),r['waste_heat_W']) for r in rows),1e-12)
check('all_cases_Euler_target_margin',max(abs(r['global_buckling_margin']/1.2-1) for r in rows),1e-12)
check('all_cases_source_B_derating',max(abs(r['assumptions']['effective_source_B_T']/(r['assumptions']['nominal_B_T']*r['on_axis_field_fraction_at_active_end'])-1) for r in rows),1e-12)
for r in rows:
 assert r['total_wall_pressure_Pa']>=r['electric_wall_pressure_Pa']
 assert r['assumptions']['heat_efficiency']<=1-r['assumptions']['radiator_K']/r['assumptions']['hot_K']

b=json.loads((ROOT/'titanium-source-refined.json').read_text())['best'][0]
eq=equilibrium(22,49,b['h0'],b['c0'],b['xi'],rtol=1e-12,nout=12001)
rr=at_wall(eq,b['xwall'],'49Ti',T=1.7,heat=.6)
check('refined_source_overlap_tight_tolerance',rr['reaction_rate_per_m_s']/b['reaction_rate_per_m_s']-1,1e-6)
check('refined_source_net_tight_sampling',rr['net_MeV_per_pn']/b['net_MeV_per_pn']-1,1e-5)
check('classical_Coulomb_orbit_integral',quad(orbit_integral,0,np.inf,epsabs=1e-11)[0]/(2/3)-1,1e-9)

endrows=json.loads((ROOT/'end-transfer-sensitivity.json').read_text())['cases']
check('end_cases_mass_ledger',max(abs(r['mass_closure_kg'])/r['total_mass_kg'] for r in endrows),1e-12)
check('end_cases_energy_ledger',max(abs(r['energy_closure_W'])/max(1,abs(r['net_W']),r['waste_heat_W']) for r in endrows),1e-12)
check('end_cases_coherent_barrier_energy',max(abs(r['additional_end_circulation_MeV_per_pn']/920-1) for r in endrows),1e-12)
for r in endrows:
 a=r['assumptions'];eta=a['end_transfer_in'];eh=a['heat_efficiency']
 expected=920*((1-eta)/eta+(1-eta)*(1-eh))
 check(f"end_recovery_B{a['nominal_B_T']}_heat{eh}_eta{eta}",r['additional_end_transfer_loss_MeV_per_pn']/expected-1,1e-10)

r0=plant(annual_replacement_fraction=0,sigma_hoop=16e9,sigma_comp=8e9)
check('zero_maintenance_bus_equals_lifecycle',r0['plant_bus_net_W_kg']/r0['net_W_kg']-1,1e-12)
xopt=brentq(lambda x:x-2*(1-np.exp(-x)),.1,3)
out={'scope':'Algebra, quadrature, convergence and internal ledgers only; no validation of the assumed reactor hardware or nuclear side factory.',
 'checks':checks,'all_passed':all(r['passed'] for r in checks),'ideal_structural_field_optimum_x':xopt,
 'ideal_optimum_field_T_at_16GPa':sqrt(xopt*MU0*16e9)}
(ROOT/'bundle-validation.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
print(json.dumps(out,indent=2))
assert out['all_passed']
