"""Conditional diamagnetic carrier / superconducting quadrupole guide screen.
Uses measured-susceptibility scale; no prediction of unmeasured carrier drag.
"""
from pathlib import Path
import json, math, hashlib
P=Path(__file__).parent; src=P/'support-topology.json'; top=json.loads(src.read_text())
mu=4*math.pi*1e-7; SB=5.670374419e-8; rho=2200.;chi=85e-6;grad=8.; bore=.25;pitch=1000.;v=1e4;R=1.1e9
acc=v*v/R; om=math.sqrt(chi*grad*grad/(mu*rho));offset=acc/om**2
lam=max(top['geometry']['axis_rail_linear_kg_m_at_1km_pitch'].values())/2
fill=.5; body_radius=math.sqrt(lam/(rho*math.pi*fill)); margin=bore-body_radius-offset
assert margin>0
rows=[]
for length in [100.,1000.,3000.,10000.]:
 t=length/v; dx=.5*acc*t*t; dv=acc*t; amplitude=math.sqrt(dx*dx+(dv/om)**2)
 rows.append({'failed_contiguous_guide_length_m':length,'transit_s':t,'tangent_departure_m':dx,'transverse_velocity_m_s':dv,'linear_recapture_oscillation_amplitude_m':amplitude,'within_reference_margin':amplitude<margin,'scope':'Complete local guide field loss; immediate ideal return to same quadrupole at next segment. End fields, real geometry and scattering excluded.'})
# Ideal cylindrical current sheet Kz=K0 cos2theta: G=mu K0/(2a).
K0=2*grad*bore/mu; Iabs=4*bore*K0
line_density=6/pitch # two counterstream apertures for each of three families
coil=[]
for Je in [1e7,3e7,1e8]:
 mass_m=Iabs/Je*8000
 outer=bore*math.exp(2*grad/(mu*Je))
 finite_area=2*(outer*outer-bore*bore)
 peak=max(mu*Je*(rr*math.log(outer/rr)/2+(rr-bore**4/rr**3)/8) for rr in [bore+(outer-bore)*i/1000 for i in range(1001)])
 coil.append({'engineering_operating_current_A_m2':Je,'conductor_kg_per_bore_m':mass_m,'conductor_kg_m2_receiver':mass_m*line_density,'finite_cos2theta_winding_outer_radius_m':outer,'finite_winding_conductor_kg_m2_receiver':finite_area*8000*line_density,'sampled_peak_winding_field_T':peak,'scope':'Complete conductor cross-section density8000. End turns, insulation, cryostat, redundancy and field-peak margins additional. Je is a requirement at65K and actual coil field/angle, not its77K self-field rating.'})
Uinside=math.pi*grad*grad*bore**4/(4*mu); Utotal=2*Uinside
cold=65.;warm=200.;cryostat_outer=coil[1]['finite_cos2theta_winding_outer_radius_m']+.05;area_per_length=2*math.pi*(.25+cryostat_outer); COP=.15*cold/(warm-cold)
cryo=[]
for ee in [.003,.01,.03]:
 Q=ee*SB*(warm**4-cold**4)*area_per_length*line_density
 cryo.append({'effective_radial_emissivity':ee,'radiative_heat_at65K_W_m2_receiver':Q,'cryocooler_COP':COP,'electric_W_m2_receiver':Q/COP,'reference_outer_cryostat_radius_m':cryostat_outer,'scope':'Inner+outer65K coil surfaces behind200K shields. Supports, leads, joints and AC losses additional; no zero-cost refrigeration.'})
# Small conducting sphere in a spatially uniform harmonic B field: induced
# eddy loss P/V=sigma_e omega² B_amp² r² /20, valid below skin-effect regime.
# This is a material-ripple diagnostic, not full guide loss.
eddy=[]
for bead in [1e-5,1e-4,1e-3,.01]:
 freq=v/100.;w=2*math.pi*freq; Bamp=.001;se=1e5
 ppkg=se*w*w*Bamp*Bamp*bead*bead/(20*rho)
 skin=math.sqrt(2/(mu*se*w))
 eddy.append({'conducting_domain_radius_m':bead,'ripple_period_m':100.,'B_ripple_amplitude_T':Bamp,'conductivity_S_m':se,'skin_depth_m':skin,'eddy_W_kg':ppkg,'below_skin_depth':bead<skin/3})
result={'scope':'Candidate support guide, independently from source launch/recovery. Scalar susceptibility lower scale, ideal quadrupole, static transverse confinement and material-ripple bounds; not a qualified10km/s bearing.',
'input_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
'parameters':{'graphite_density_kg_m3':rho,'susceptibility_magnitude_assumed':chi,'gradient_T_m':grad,'bore_radius_m':bore,'longitudinal_carrier_fill_fraction':fill,'equatorial_speed_m_s':v,'total_bore_length_per_receiver_area_m_m2':line_density},
'confinement':{'centripetal_acceleration_m_s2':acc,'transverse_omega_rad_s':om,'equilibrium_offset_m':offset,'heaviest_family_each_counterstream_kg_m':lam,'carrier_radius_m_at_fill_fraction':body_radius,'available_transverse_margin_m':margin,'bore_edge_T':grad*bore},
'fault_recapture':rows,'coil_current_and_mass':{'ideal_sheet_amplitude_A_m':K0,'sum_absolute_axial_current_A':Iabs,'rows':coil,'field_energy_inside_J_per_bore_m':Uinside,'field_energy_inside_and_outside_J_per_bore_m':Utotal,'field_energy_J_m2_receiver':Utotal*line_density},
'cryogenic_radiation':cryo,'eddy_ripple_diagnostic':eddy,
'sources':{'susceptibility_and_levitation':'https://pmc.ncbi.nlm.nih.gov/articles/PMC10304365/','insulating_diamagnetic_composite_experiment':'https://arxiv.org/abs/2312.01660','eddy_damping_patterning_experiment':'https://arxiv.org/abs/2211.08764','REBCO_manufacturer':'https://superpower-inc.com/specification.aspx','REBCO_field_temperature_measurements':'https://arxiv.org/abs/1512.01930'}}
assert abs(mu*K0/(2*bore)-grad)<1e-12
Path(__file__).with_name('guide-component.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
print(json.dumps({k:result[k] for k in ['confinement','fault_recapture','coil_current_and_mass','cryogenic_radiation','eddy_ripple_diagnostic']},indent=2))
