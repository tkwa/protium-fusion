"""Finite-placement, closed-flux superconducting cleanup-array acceptance screen.
Exact periodic alternating-line-current vacuum field; engineering cable, buffer,
shield and finite module qualifications are explicitly separate.
"""
from pathlib import Path
import json,math,hashlib
import numpy as np
P=Path(__file__).resolve().parent;C=P.parents[1]
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
source=C/'final-family-dispatch/sigma350/input.json';d=json.loads(source.read_text());plant=d['plant']
place=json.loads((P/'placement.json').read_text());B=place['selected_actual_duty']['B_T']
buf=json.loads((C/'guide-buffer-materials.json').read_text());gamma=json.loads((C/'guide-gamma.json').read_text())
mu=4*math.pi*1e-7;SB=5.670374419e-8;kappa=gamma['carbon_minimum_mass_energy_absorption_m2_kg_through20MeV'];h=buf['usable_delta_h65to70_J_kg'];rhoPE=900.;rhoC=2000.;rhoCable=8000.;cop=.15*65/(200-65)
area=plant['area_m2'];mean=plant['cold_heat_components_W']['captured_source_ineligible_for_hot_engine']/area
# Conservative full-stage delayed-heat photon reserve, including beta kinetic
# energy and a10% Doppler reserve; prompt energy kept separately. This is not a
# time-resolved secondary transport calculation. Main-wall transmission uses the
# same local-deposition <=20MeV energy-attenuation scenario as the guide packet.
pulses={}
for key,s in d['loaded_stages'].items():
 nonnu_delayed=s['Q_cold_MeV']-s['nu_cold_MeV'];prompt=s['radiated_MeV'];E=(1.1*nonnu_delayed+prompt)*s['attempts_per_pulse']*1.602176634e-13/area
 pulses[s['label']]=E*math.exp(-kappa*350.)
Emax=max(pulses.values());checks=[]

def field(p,depth,Je):
 k=2*math.pi/p
 if k*depth>15:return None
 I=B*p*math.cosh(k*depth)/mu;radius=math.sqrt(I/(math.pi*Je));line=2/p
 if radius>=p/4:return None
 U=mu*I*I/(2*math.pi*p)*(math.log(p/(math.pi*radius))+.25)
 return {'current_A_per_rail':I,'cable_radius_m_uniform_bundle':radius,'cable_mass_kg_m2':I/Je*line*rhoCable,'return_field_energy_J_m2':U,'energy_equivalent_field_volume_m':2*mu*U/B**2,'front_Bmin_T':B,'front_Bmax_T':B/math.tanh(k*depth),'cable_surface_max_field_T_uniform_bundle':mu*I/(p*math.sin(k*radius)),'line_density_m_m2':line}

def solve(p,Je,shield,capacity=1.):
 # Back of retained175mm wall,50mm warm clearance, shield outer radius,50mm
 # cold-to-hot insulation envelope. All annuli must fit at their actual depth.
 column_clearance=.175+.05;depth=column_clearance+.05+shield/rhoC;buffer=0.
 status='no_fixed_point_within_bounds'
 for iteration in range(10000):
  f=field(p,depth,Je)
  if f is None:break
  rc=math.sqrt(f['cable_radius_m_uniform_bundle']**2+buffer/(math.pi*rhoPE*f['line_density_m_m2']))+.05
  rout=rc+shield/rhoC
  # Radial rays and an outward-isotropic crossing reference are distinct.
  # Overlap is rejected, not repaired by silently capping an impossible geometry.
  if 2*rout>=p/2:status='neighbor_shield_envelopes_overlap';break
  intercept=min(1.,math.pi*rout*f['line_density_m_m2'])
  trans=math.exp(-kappa*shield);newB=Emax*intercept*trans/(h*capacity);newd=column_clearance+rout
  if max(abs(newB-buffer)/max(1,newB),abs(newd-depth)/max(1,newd))<1e-11:
   buffer=newB;depth=newd;status='finite_fixed_point';break
  buffer=.5*(buffer+newB);depth=.5*(depth+newd)
 base={'pitch_m':p,'Je_A_mm2':Je/1e6,'shield_radial_column_kg_m2':shield,'buffer_accessible_capacity_fraction':capacity,'status':status,'iterations':iteration+1}
 if status!='finite_fixed_point':return base
 f=field(p,depth,Je);rc=math.sqrt(f['cable_radius_m_uniform_bundle']**2+buffer/(math.pi*rhoPE*f['line_density_m_m2']))+.05;rout=rc+shield/rhoC;line=f['line_density_m_m2'];intercept=min(1.,math.pi*rout*line);trans=math.exp(-kappa*shield)
 shieldmass=math.pi*rhoC*(rout*rout-rc*rc)*line;Q=mean*intercept*trans;ir=.01*SB*(200**4-65**4)*2*math.pi*rc*line
 # Finite racetrack: straight length100p, semicircular end returns, radiusp/4.
 endfrac=math.pi/400;wire=f['cable_mass_kg_m2']*(1+endfrac)
 total=wire+shieldmass+buffer
 checks.extend([abs(buffer*h*capacity/(Emax*intercept*trans)-1)<1e-9,abs(depth-(column_clearance+rout))<1e-8,2*rout<p/2,rc>f['cable_radius_m_uniform_bundle']])
 return dict(base,**f,coil_center_depth_below_hot_face_m=depth,cold_envelope_radius_m=rc,hot_shield_outer_radius_m=rout,outward_isotropic_outer_intercept_fraction=intercept,radial_projected_outer_fraction=2*rout*line,conditional_photon_transmission=trans,PE_buffer_kg_m2=buffer,graphite_halo_kg_m2=shieldmass,cable_with_semicircular_returns_kg_m2=wire,subtotal_cable_shield_buffer_kg_m2=total,subtotal_fraction_existing_plant_mass=total/(plant['mass_kg']/area),mean_photon_heat65K_W_m2=Q,photon_refrigeration_electric_W_m2=Q/cop,IR200to65K_W_m2_at_emissivity0p01=ir,IR_refrigeration_electric_W_m2=ir/cop,hot_halo_to200K_W_m2_at_emissivity1em4=1e-4*SB*(1000**4-200**4)*2*math.pi*rc*line,maximum_pulse_absorbed65K_J_m2=Emax*intercept*trans,buffer_ionizing_dose_Gy_year_if_all_heat_deposits_in_buffer=Q/buffer*365.25*86400.,stored_field_energy_whole_array_J=f['return_field_energy_J_m2']*area,finite_module_length_m=100*p,end_field_region_fraction_reference=.02)

rows=[solve(p,Je,shield)for p in[1.,2.,5.,10.,20.,50.,100.]for Je in[30e6,100e6]for shield in[0.,500.,1000.,2000.,4000.]]
# Independent Biot-Savart truncated-array check against the csc field formula.
for p,depth in[(2.,.5),(10.,1.),(50.,2.)]:
 I=B*p*math.cosh(2*math.pi*depth/p)/mu
 for x in[0.,p/8,p/4]:
  n=np.arange(-20000,20001);dx=x-n*p/2;sgn=(-1.)**n;Bx=float(np.sum(sgn*depth/(dx*dx+depth*depth)))*mu*I/(2*math.pi);Bz=float(np.sum(sgn*dx/(dx*dx+depth*depth)))*mu*I/(2*math.pi)
  exact=mu*I/p/math.sqrt(math.sin(2*math.pi*x/p)**2+math.sinh(2*math.pi*depth/p)**2)
  checks.append(abs(math.hypot(Bx,Bz)/exact-1)<2e-8)
out={'scope':__doc__,'code_sha256':sha(__file__),'input_sha256':{str(p):sha(p)for p in[source,P/'placement.json',C/'guide-buffer-materials.json',C/'guide-gamma.json']},'parameters':{'Bmin_front_T':B,'complete_cable_density_kg_m3':rhoCable,'PE_effective_density_kg_m3':rhoPE,'PE_usable65to70_J_kg':h,'COP65to200':cop,'minimum_carbon_mu_en_m2_kg_up_to20MeV':kappa,'mean_current_family_missed_heat_W_m2':mean,'pulse_upper_fluences_after350kg_m2_J_m2':pulses,'maximum_pulse_upper_fluence_J_m2':Emax,'existing_plant_mass_kg_m2':plant['mass_kg']/area,'existing_export_W_m2':plant['net_W']/area},'equations':{'field':'B(x,z)=mu0 I/p / sqrt(sin(2pi x/p)^2+sinh(2pi z/p)^2); alternate +/-I rails spacedp/2. MinB atx=p/4; maxB atx=0. Flux returns in alternating cells, decay lengthp/(2pi).','I':'I=Bmin*p*cosh(2pi*d/p)/mu0','field_energy':'U/A=mu0*I^2/(2pi*p)*(ln(p/(pi*r))+1/4). Uniform-current circular cable, r<p/4; includes both half-spaces and internal field.','cable':'r=sqrt(I/(pi Je)); line density2/p. Semicircular finite-module returns addpi/400 to straight-cable mass for length100p.','placement':'d=.175+.05+router; router=sqrt(rcable^2+Mbuffer/(pi*rhoPE*line))+.05+Sigmahalo/rhoC.','pulse':'Mbuffer*delta_h=Emax*min(1,pi*router*line)*exp(-kappa*Sigmahalo). Coupled to coil current and placement.'},'rows':rows,'checks':{'number':len(checks),'passed':sum(checks)},'limits':['This is a periodic vacuum-field and calorimetric geometry screen. A complete magnet needs finite end fields, winding stress, cryostats, thermal straps, joints, protection, feeds, local radiation transport and quench validation. These are not included in the subtotal.','Je30/100A/mm2 are required complete-cable operating ratings at65K and the computed local field/angle, not measured ratings assigned universally to REBCO. Senatore etal2016 measures strong temperature, field, angle and manufacturer dependence.','Local electron cells must follow the nonuniform field vector or have explicitly qualified cross-field confinement. The field direction is not everywhere normal to the wall; near horizontal field regions need oriented cells and actual access geometry.','Cooled rails are behind retained wall+clearance and the entire local shield/buffer radius. No full cold plane views the hot cavity, and no1000K superconductivity is assumed.','Photon mean is the actual350kg/m2 family cold-equipment allocation. The pulse upper reserve treats all non-neutrino delayed heat as<=20MeV photons with10%Doppler allowance; not measured pulse leakage. Applying energy-absorption exponential requires local charged-secondary deposition and joint/streaming control.','Outer-halo interception gives a deliberately conservative photon load; actual cold projected area is smaller. The local halo may return heat to the hot side, but this calculation takes no hot-engine recovery credit.','Only a single largest-pulse stored-energy upper bound is used. Real time-resolved cooling between and during pulses, adjacent coils and common buffers may alter the stock. A5K excursion is an acceptance condition, not qualified quench margin.','The full calorimetric capacity is assumed pulse-accessible. Cryogenic polymer heat diffusion, local dose/radiolysis and finite thermal spreaders need qualification.','Dense overlapping shield envelopes are rejected. Other layouts, shared planar shielding and distributed remote buffers can improve geometry but require a new field/heat account.','Any extra refrigerator work and hot-to-cold heat bypass must be charged once in the radiator-limited station ledger. Existing auxiliary/services budgets are allocations, not automatically spare capacity.']}
assert all(checks);(P/'coil-array.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
finite=[x for x in rows if x['status']=='finite_fixed_point'];print(json.dumps({'parameters':out['parameters'],'checks':out['checks'],'finite_count':len(finite),'lightest':sorted(finite,key=lambda x:x['subtotal_cable_shield_buffer_kg_m2'])[:8]},indent=2))
