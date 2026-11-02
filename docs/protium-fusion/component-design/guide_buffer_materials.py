"""Abundant-material calorimetric branch of the frozen guide-gamma screen.
Uses measured 65--70 K polyethylene Cp, not room-temperature properties.
This is a heat-capacity/geometry acceptance calculation, not quench validation.
"""
from pathlib import Path
import json,math,hashlib
P=Path('/home/tkwa/.codex/worktrees/protium-fusion/website/docs/protium-fusion/component-design')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
gamma=json.loads((P/'guide-gamma.json').read_text())
g=json.loads((P/'guide-component.json').read_text())
line=gamma['guide_line_density_m_m2']; coil=g['coil_current_and_mass']['rows'][1]['finite_cos2theta_winding_outer_radius_m']
oldrad=gamma['original_cryostat_radius_m'];cop=gamma['COP65to200'];old_ir=g['cryogenic_radiation'][1]['electric_W_m2_receiver']
kappa=gamma['carbon_minimum_mass_energy_absorption_m2_kg_through20MeV'];area=4*math.pi*(1.1e9)**2
# Chang & Bestul (1973), Table 1, SRM1475 as received, series V/II.
# Cp per base mole CH2 (14.027 g). These points bracket the whole interval.
data=[(61.44,6.008),(67.41,6.662),(73.12,7.252)];molar=.014027

def cp(T):
 for (t0,c0),(t1,c1) in zip(data[:-1],data[1:]):
  if t0<=T<=t1:return (c0+(T-t0)/(t1-t0)*(c1-c0))/molar
 raise ValueError('No extrapolation')
knots=[65.,67.41,70.]
pe_h=sum((b-a)*(cp(a)+cp(b))/2 for a,b in zip(knots[:-1],knots[1:]));margin=.95
materials=[
 {'name':'OFHC_copper_benchmark','density_kg_m3':8960.,'usable_delta_h_J_kg':gamma['copper65to70_sensible_J_kg_after5pct_margin'],'capacity_fraction':1.,'basis':'Frozen guide-gamma NIST fit with its existing 5% margin'},
 {'name':'linear_polyethylene_effective_density900','density_kg_m3':900.,'usable_delta_h_J_kg':pe_h*margin,'capacity_fraction':1.,'basis':'Measured SRM1475 Cp interpolation; 5% engineering reserve; effective 900 kg/m3 includes volume allowance, not a measured 65K density'},
 {'name':'linear_polyethylene_effective_density600','density_kg_m3':600.,'usable_delta_h_J_kg':pe_h*margin,'capacity_fraction':1.,'basis':'Same measured Cp; explicitly lower packing-density sensitivity'},
 {'name':'linear_polyethylene_half_accessible_capacity','density_kg_m3':900.,'usable_delta_h_J_kg':pe_h*margin*.5,'capacity_fraction':.5,'basis':'Only half calorimetric capacity accessible during pulse; sensitivity, not a measured diffusion/damage loss'},
]

def solve(row,mat,ang,shield=2000.):
 Emax=max(row['pulse_photon_upper_fluence_J_m2_after_main_wall'].values());F=row['mean_missed_hot_to_cold_W_m2'];trans=math.exp(-kappa*shield)
 def geom(B):
  cold=math.sqrt(coil**2+B/(math.pi*mat['density_kg_m3']*line))+.05
  outer=cold+shield/2000.; f=min(1.,2*outer*line*ang)
  return cold,outer,f
 B=0.
 for it in range(1000):
  cold,outer,f=geom(B);new=Emax*f*trans/mat['usable_delta_h_J_kg']
  if abs(new-B)<1e-12*max(1,new):B=new;break
  B=new
 else:raise AssertionError('no convergence')
 cold,outer,f=geom(B);E=Emax*f*trans;assert abs(B*mat['usable_delta_h_J_kg']-E)<1e-9*E
 Mshield=math.pi*2000*(outer**2-cold**2)*line;Q=F*f*trans
 ir=old_ir*((.25+cold)/(.25+oldrad)-1)
 # Assuming all absorbed heat resides in the buffer overstates buffer dose if
 # coils/spreaders take a share; calorimetry alone does not predict radiolysis.
 dose=Q/B*365.25*86400
 return {'receiver_kg_m2':row['receiver_kg_m2'],'material':mat['name'],'angular_factor':ang,'shield_column_kg_m2_radial':shield,'cold_radius_m':cold,'outer_shield_radius_m':outer,'outer_intercept_fraction':f,'buffer_kg_m2':B,'buffer_total_kg':B*area,'local_graphite_shield_kg_m2':Mshield,'shield_plus_buffer_kg_m2':B+Mshield,'whole_pulse_absorbed_J_m2':E,'mean_load65K_W_m2':Q,'photon_refrigeration_electric_W_m2':Q/cop,'extra_IR_refrigeration_electric_W_m2':ir,'total_extra_electric_W_m2':Q/cop+ir,'local_hot_to200K_leak_W_m2_at_effective_emissivity1em4':1e-4*5.670374419e-8*(1000**4-200**4)*2*math.pi*cold*line,'buffer_mean_dose_if_all_heat_is_ionizing_Gy_year':dose,'buffer_recycle_kg_m2_year_at_chosen_dose_limit':{str(lim):B*dose/lim for lim in[1e5,1e6,1e7]},'buffer_storage_identity_relative_error':abs(B*mat['usable_delta_h_J_kg']-E)/E}
rows=[solve(r,m,a)for r in gamma['rows'] for m in materials for a in[math.pi/2,3.]]
checks=[]
for r in gamma['rows']:
 old=next(x for x in r['local_shield_buffer_scenarios']if x['shield_column_kg_m2_radial']==2000 and x['angular_factor']==math.pi/2)
 new=next(x for x in rows if x['receiver_kg_m2']==r['receiver_kg_m2']and x['material']=='OFHC_copper_benchmark'and x['angular_factor']==math.pi/2)
 for oldkey,newkey in [('65K_copper_buffer_kg_m2_receiver','buffer_kg_m2'),('shield_plus_buffer_kg_m2_receiver','shield_plus_buffer_kg_m2')]:
  err=abs(old[oldkey]-new[newkey])/old[oldkey];assert err<1e-9;checks.append(err)
# Diffusion-length acceptance envelope: l <= sqrt(k*t/(rho*cp)).
# k is a declared required effective composite conductivity, not polymer data.
diff=[{'effective_k_W_m_K':k,'time_s':t,'max_diffusion_length_m':math.sqrt(k*t/(900*(pe_h/5)))}for k in[.1,1.,10.]for t in[1.,10.,100.]]
out={'scope':__doc__,'code_sha256':sha(__file__),'inputs':{str(P/n):sha(P/n)for n in['guide-gamma.json','guide-component.json']},'polyethylene_primary_source':'https://nvlpubs.nist.gov/nistpubs/jres/77A/jresv77An4p395_A1b.pdf','polyethylene_source_details':'S.S.Chang and A.B.Bestul (1973), Table 1, SRM1475 as received, measured Cp; base mole CH2=14.027 g. Reported sample room-temperature density954 kg/m3.','tabulated_T_K_Cp_J_mol_K':data,'Cp65_J_kg_K':cp(65),'Cp70_J_kg_K':cp(70),'integrated_Cp65to70_J_kg':pe_h,'engineering_capacity_fraction':margin,'usable_delta_h65to70_J_kg':pe_h*margin,'materials':materials,'rows':rows,'diffusion_acceptance_scenarios':diff,'checks':{'copper_benchmark_reproduction_checks':len(checks),'maximum_relative_error':max(checks),'all_storage_identity_max_relative_error':max(r['buffer_storage_identity_relative_error']for r in rows)},'limits':['This demonstrates an abundant solid C/H calorimetric branch. It does not establish pulse heat access, conductor quench margin, thermal contact, cryogenic mechanical integrity or radiation lifetime.','All original guide-gamma photon transport assumptions, source geometry, chosen 5K margin and cooling COP remain unchanged.','Effective packing densities900/600 and half-accessible capacity are explicit design sensitivities, not measurements or credible intervals.','Heat is transferred to an abundant buffer through finite local spreaders/contacts, whose mass, conductivity and electrical insulation must be allocated within services. Copper in conductor stabilizer is retained; only the extra calorimetric stock is replaced.','Cumulative absorbed dose is reported as an acceptance curve. The chosen0.1/1/10MGy limits are not asserted polymer lifetimes; radiation chemistry/outgassing and heat capacity must be qualified. All retired material may be recycled; fresh feed is not equated with gross replacement throughput.','No methane numerical heat capacity is imported from an inaccessible historical table; solid polyethylene supplies directly verified primary calorimetry.','No manufacturing energy is inferred from atomization; any restoration process retains the user\'s >=50% useful free-energy efficiency premise and needs its own throughput/cost boundary.']}
(P/'guide-buffer-materials.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
print(json.dumps({'Cp':[cp(65),cp(70)],'integral_Jkg':pe_h,'usable_Jkg':pe_h*margin,'checks':out['checks'],'main_rows':[r for r in rows if r['angular_factor']==math.pi/2 and r['material']in['OFHC_copper_benchmark','linear_polyethylene_effective_density900']]},indent=2))
