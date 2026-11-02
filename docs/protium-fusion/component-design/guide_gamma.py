"""Guide-interception, local hot shielding and finite65K buffer acceptance screen.
No radiation Monte Carlo or quench simulation. All shield transfer factors inherit
localized charged-secondary/<=20MeV qualifications from photon-containment.
"""
from pathlib import Path
import json,math,hashlib
import numpy as np
P=Path(__file__).parent;trade=Path('/home/tkwa/src/protium-fusion/output/current-rate128-receiver-tradeoff.json');gp=P/'guide-component.json';sp=P/'receiver-tradeoff-support.json';cp=P/'photon-containment/containment.json'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
g=json.loads(gp.read_text());s=json.loads(sp.read_text());data=json.loads(trade.read_text());table=json.loads(cp.read_text())['nist_table'];kappa=min(r[2]for r in table)*.1
line=g['parameters']['total_bore_length_per_receiver_area_m_m2'];coil=g['coil_current_and_mass']['rows'][1]['finite_cos2theta_winding_outer_radius_m'];oldrad=g['cryogenic_radiation'][1]['reference_outer_cryostat_radius_m'];cop=g['cryogenic_radiation'][1]['cryocooler_COP'];old_ir=g['cryogenic_radiation'][1]['electric_W_m2_receiver']
coef=[-1.91844,-.15973,8.61013,-18.996,21.9661,-12.7328,3.54322,-.3797]
T=np.linspace(65,70,1001);cpCu=10**sum(c*np.log10(T)**i for i,c in enumerate(coef));buffer_J_kg=float(np.trapezoid(cpCu,T))*.95 #5% fit margin; not composite heat capacity.
rows=[]
for r in data['rows']:
 if r['receiver_sigma_kg_m2']not in[500,750,1000] or r['additional_missed_hot_fraction_HplusU']!=.005 or r['missed_hot_routing']!='cold_equipment':continue
 p=r['plant'];F=r['heat_in_cold_equipment_W']/p['area_m2'];sig=r['receiver_sigma_kg_m2']
 # Loose pulse reserve: all future non-neutrino decay heat shifted to<=20MeV
 # photons with1.1 Doppler energy reserve; prompt also<=20MeV. Extra missed
 # .005 heat is carried as a separate photon-like acceptance scenario.
 pulses={name:(1.1*e['delayed_J_m2']+e['prompt_J_m2'])*math.exp(-kappa*sig)+.005*e['total_smooth_J_m2'] for name,e in s['pulse_fluences'].items()}
 Emax=max(pulses.values());configs=[]
 for shield in[0.,1000.,2000.]:
  for ang in[1.,math.pi/2,3.]:
   def geometry(buffer):
    cold_outer=math.sqrt(coil*coil+buffer/(math.pi*8960.*line))+.05
    outer=cold_outer+shield/2000
    fraction=min(1.,2*outer*line*ang)
    return cold_outer,outer,fraction
   B=0.
   for i in range(100):
    cold,outer,coverage=geometry(B);target=Emax*coverage*math.exp(-kappa*shield)/buffer_J_kg
    if abs(target-B)<1e-10*max(1,target):B=target;break
    B=target
   cold,outer,coverage=geometry(B);assert abs(B-Emax*coverage*math.exp(-kappa*shield)/buffer_J_kg)<1e-8*max(1,B);trans=math.exp(-kappa*shield);Q65=F*coverage*trans
   shield_mass=math.pi*2000*(outer**2-cold**2)*line
   ir_increment=old_ir*((.25+cold)/(.25+oldrad)-1)
   configs.append({'shield_column_kg_m2_radial':shield,'angular_factor':ang,'cold_cryostat_outer_radius_m_with_buffer':cold,'shield_outer_radius_m':outer,'projected_outer_intercept_fraction':coverage,'conditional_transmission_factor':trans,'local_hot_shield_kg_m2_receiver':shield_mass,'65K_copper_buffer_kg_m2_receiver':B,'shield_plus_buffer_kg_m2_receiver':shield_mass+B,'mean_photon_load65K_W_m2':Q65,'additional_refrigerator_electric_W_m2':Q65/cop,'increment_existing200to65K_IR_refrigeration_W_m2_due_buffer_area':ir_increment,'buffer_energy_J_m2':B*buffer_J_kg,'hot_shield_to200K_leak_W_m2_at_effective_emissivity1em4':1e-4*5.670374419e-8*(1000**4-200**4)*2*math.pi*cold*line if shield else 0.,'scope':'Annular copper buffer enlarges cryostat and interception; geometry solved self-consistently. Copper stabilizer already present is not credited as free buffer. Continuous closed-ring local graphite shield; holes/joins/secondary ranges and energy transfer to buffer require qualification.'})
 bare=[]
 for ang in[1.,math.pi/2,3.,10.]:
  f=min(1.,2*oldrad*line*ang);bare.append({'angular_factor':ang,'intercept_fraction':f,'Q65_W_m2':F*f,'refrigerator_electric_W_m2':F*f/cop,'maximum_effective_intercept_fraction_for1W_m2_refrigerator':cop/F})
 rows.append({'receiver_kg_m2':sig,'mean_missed_hot_to_cold_W_m2':F,'bare_existing_geometry_average_only':bare,'pulse_photon_upper_fluence_J_m2_after_main_wall':pulses,'local_shield_buffer_scenarios':configs})
out={'scope':__doc__,'code_sha256':sha(__file__),'input_sha256':{str(p):sha(p)for p in[trade,gp,sp,cp]},'guide_line_density_m_m2':line,'original_cryostat_radius_m':oldrad,'radial_projected_coverage':2*oldrad*line,'outward_isotropic_intensity_angular_factor':math.pi/2,'angular_derivation':'For cylinder axis tangent to a plane, intercepted isotropic outward intensity per length is pi²*r*I while plane flux is pi*I; equivalent width pi*r versus2*r for radial rays. Grazing/concentrated fields need their actual distribution.','COP65to200':cop,'carbon_minimum_mass_energy_absorption_m2_kg_through20MeV':kappa,'copper_cp65_J_kg_K':float(cpCu[0]),'copper65to70_sensible_J_kg_after5pct_margin':buffer_J_kg,'copper_source':'https://trc.nist.gov/cryogenics/materials/OFHC%20Copper/OFHC_Copper_rev1.htm','carbon_source':'https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z06.html','rows':rows,'limitations':['All intercepted energy deposited at65K is conservative for mean heat; the existing200K shields may absorb part.','The continuous local shield factor is an energy-transport scenario under local charged deposition, no upscatter, sources outside the shield and minimum radial column; it is not a completed secondary-particle or joint-leakage calculation.','Pulse reserve is deliberately broader than the parent cycle-average photon bound. Time-resolved refrigeration removal and existing coil heat capacity are not credited, but ideal buffer thermalization is assumed.','Five-kelvin rise is a chosen acceptance margin, not measured magnet quench headroom; conductor critical-current margin and local diffusion need qualification.','Local hot-shield mass and cold buffer are finite service allocations. Keeping parent missed heat charged at200K while returning some to hot side is conservative; do not double credit its heat.','The additional electrical demand must be itemized within ordinary auxiliaries or debited once; at radiator saturation additional unallocated electrical dissipation changes cadence.']}
(P/'guide-gamma.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
print(json.dumps({'geometry':out['radial_projected_coverage'],'cp65':out['copper_cp65_J_kg_K'],'buffer_Jkg':buffer_J_kg,'rows':[{'sigma':x['receiver_kg_m2'],'bare':x['bare_existing_geometry_average_only'][:2],'shield2000_diffuse':next(y for y in x['local_shield_buffer_scenarios']if y['shield_column_kg_m2_radial']==2000 and y['angular_factor']==math.pi/2)}for x in rows]},indent=2))
