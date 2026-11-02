"""Conditional charged-grain neutralizer and scalable condensed-source ledger.
No nonlinear sheath/secondary-emission simulation or grain lifetime claim.
"""
from pathlib import Path
import json,math,hashlib
E=1.602176634e-19;EP=8.8541878128e-12;ME=9.1093837139e-31;MP=1.67262192595e-27;C=299792458.;R=1.1e9;V=2.4e6;M=1.1309733552923255e13;DUR=700
checks=[];rows=[]
for rho in [1000.,132.]:
 for a in [.001,.01,.05]:
  Es=1e9;qm=3*EP*Es/(rho*a);mg=4*math.pi*rho*a**3/3;q=qm*mg;U=.5*qm*Es*a
  for ne in [1e12,1e13,1e14,1e15]:
   te=1.;ve=math.sqrt(8*te*E/(math.pi*ME));j=E*ne*ve/4;dt=EP*Es/j;length=V*dt
   debye=math.sqrt(EP*te/(ne*E));I=M*qm/DUR/1e5;npgr=M/DUR/1e5/(V*math.pi*100**2)/mg
   column=ne*MP*length;heat=.5*column*V*V/(4*rho*a/3)
   residual=5/(Es*a);Qtot=M*qm*residual;rc=2000.;coulomb=3*Qtot**2/(20*math.pi*EP*rc);K=.5*M*V*V
   rows.append(dict(rho_kg_m3=rho,grain_radius_m=a,surface_field_V_m=Es,grain_mass_kg=mg,grain_charge_C=q,q_m_C_kg=qm,grain_potential_V=Es*a,field_energy_J_kg=U,plasma_ne_m3=ne,plasma_Te_eV=te,geometric_thermal_electron_current_A_m2=j,neutralizer_time_s=dt,neutralizer_length_m=length,Debye_length_m=debye,grain_radius_over_Debye=a/debye,plasma_supply_volume_per_grain_m3=q/E/ne,plasma_supply_radius_per_grain_m=(3*q/E/ne/(4*math.pi))**(1/3),mean_grain_spacing_m=npgr**(-1/3),mean_current_per_1e5_outlet_A=I,advective_electron_current_through_100m_radius_A=ne*E*V*math.pi*100**2,stationary_H_plasma_sweep_mass_column_kg_m2=column,stationary_H_plasma_geometric_energy_upper_J_kg=heat,assumed_final_floating_potential_magnitude_V=5,residual_net_grain_charge_fraction_at_5V=residual,uniform_assembly_residual_charge_energy_over_K_at_R2000m=coulomb/K,field_energy_over_c_directed_momentum_upper_m_s=U/C,neutralization_electron_count=q/E))
   checks += [math.isclose(4*math.pi*a*a*j*dt,q,rel_tol=1e-12),math.isclose(U,.5*q*Es*a/mg,rel_tol=1e-12),dt>0]
# Recipes are volume-additive illustrations. N matrices densities are parameters.
source=[]
for cat,A,capacity,rm in [('C12',12,2,1000.),('C13',13,2,1000.),('N14',14,3,850.),('N15',15,3,850.)]:
 for h in [.5,1,2,3,4,8,16]:
  # Below capacity: corresponding matrix and catalyst-rich filler absorbed into chosen rho.
  excess=max(0,h-capacity);frac=excess/(A+h);rho=1/((1-frac)/rm+frac/70.8);vf=frac*rho/70.8
  alpha=.02182;sigma=5.670374419e-8;a=.01;dose=3*alpha*sigma*1000**4*R/(rho*a*V)
  surface_scale=math.sqrt(rho/1000);Es=1e9*surface_scale
  source.append(dict(catalyst=cat,H_atoms_per_catalyst=h,condensed_matrix_H_capacity=capacity,chosen_matrix_density_kg_m3=rm,excess_H2_mass_fraction=frac,excess_H2_volume_fraction=vf,volume_additive_density_kg_m3=rho,fuel_volume_m3_per_kg=1/rho,H2_component_density_kg_m3=70.8,reference_Al_IR_dose_1cm_J_kg=dose,grain_radius_for_same_IR_dose_as_1cm_at1000kg_m3_m=.01*1000/rho,surface_field_for_constant13p3kJ_perkg_V_m=Es,acceleration_length_ratio_at_that_field_and_fixedradius=surface_scale,skin_mass_ratio_if_1GV_and1GPa_skin_alone=3*2000/rho*(.5*EP*1e18)/(2*1e9)))
  checks += [0<=frac<1,0<=vf<1,rho<=max(rm,70.8)*(1+1e-12)]
out={'scope':'Parametric geometric-flux neutralization; initial nonlinear sheath, MeV electron stopping/secondary emission, actual final charge and collector transport uncomputed. Cold mixture volume additive at chosen matrix densities; no H2 containment/phase closure.','code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'neutralizer_rows':rows,'source_recipe_rows':source,'source_urls':['https://doi.org/10.1063/1.5128521','https://doi.org/10.1103/PhysRevLett.113.035002','https://nonneutral.pppl.gov/pdfpapers2010/Kaganovich.PoP.2010.pdf','https://www.energy.gov/sites/default/files/2021-10/lh2-storage-handling-demonstrations.pdf'],'checks':{'count':len(checks),'passed':sum(checks),'all_passed':all(checks)}}
Path(__file__).with_name('launcher-neutralizer.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out['checks']))
