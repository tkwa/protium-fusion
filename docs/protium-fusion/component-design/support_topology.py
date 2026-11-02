"""Analytic three-axis latitude-ring support and common breathing screens.
No claim of complete guide or nonradial stability design.
"""
from pathlib import Path
import json, math, hashlib
P=Path(__file__).resolve().parents[1]
base=P/'component-design/sector-account.json'; j=json.loads(base.read_text());r=j['reference'];s=j['support']
R=r['radius_m'];A=r['area_m2'];sigma=s['total_areal_mass_kg_m2'];v=1e4;omega=v/R
aself=s['self_acceleration_m_s2']/R;n2=s['maximum_tidal_acceleration_m_s2']/(3*R)
# Static Hill-frame external field diag(3n²,0,-n²). Desired support is
# diag(aself-3n²,aself,aself+n²)*r. Stream families supply b_i(I-aa')*r.
b=[(aself+4*n2)/2,(aself-2*n2)/2,(aself-4*n2)/2]
assert min(b)>0
sig=[sigma*x/omega**2 for x in b];D=[aself-3*n2,aself,aself+n2]
checks=[]
for lat in [-1.,-.6,0.,.4,1.]:
 for angle in [0.,.4,1.5,2.7]:
  u=[math.sqrt(1-lat*lat)*math.cos(angle),math.sqrt(1-lat*lat)*math.sin(angle),lat]
  F=[sum(sig[i]*omega**2*R*u[k] for i in range(3) if i!=k) for k in range(3)]
  target=[sigma*D[k]*R*u[k] for k in range(3)]
  checks.append(max(abs(F[k]-target[k]) for k in range(3)))
Ms=sum(sig)*A;K=Ms*v*v/3;period=2*math.pi/omega
incoming={'A':2.5945081e19,'B1':2.45497515e17,'B2':2.80190470e19,'C13':2.53252699e19}
rows=[]
for k,x in s['outgoing_capture'].items():
 T=x['interval_s'];imp=x['momentum_fluence_Pa_s']+incoming[k]/A;dv=imp/sigma
 rows.append({'stage':k,'incoming_reaction_impulse_Pa_s':incoming[k]/A,'outgoing_capture_impulse_Pa_s':x['momentum_fluence_Pa_s'],'total_radial_impulse_Pa_s':imp,'delta_v_m_s':dv,'period_s':T,'constant_mean_support_free_mass_displacement_range_m':dv*T/8,'mean_added_radial_pressure_Pa':imp/T})
# Uniform radial homologous motion with individually conserved ring angular
# momenta has centrifugal a~R^-3 and self gravity~R^-2: w0²=aself.
# For periodic delta-v impulses on linear oscillator about its mean, each
# Fourier acceleration harmonic is 2*dv/T. Range <=2 sum |a_n|/|w_n²-w0²|.
# Tail estimate bounds 1/(n²-r²) by 1/[n²(1-r²/(N+1)²)].
w0=math.sqrt(aself);N=10000
for x in rows:
 wf=2*math.pi/x['period_s'];rr=(w0/wf)**2
 assert rr<1
 total=sum(1/(n*n-rr) for n in range(1,N+1))
 tail=1/N/(1-rr/(N+1)**2)
 x['linear_radial_oscillator_range_bound_m']=4*x['delta_v_m_s']/x['period_s']/wf**2*(total+tail)
result={'scope':'Three cold guide layers, latitude rings around orthogonal Hill-frame axes, equal opposite streams within each ring. Exact continuum support force in linear tidal field; global shear modes, guides, fault transfer and discrete coverage remain open.',
'input_sha256':hashlib.sha256(base.read_bytes()).hexdigest(),'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
'geometry':{'equatorial_speed_m_s':v,'common_angular_speed_rad_s':omega,'counterstreams_per_ring':2,'radially_separate_axis_families':3,'axis_surface_stream_kg_m2':dict(zip(['solar_radial_x','orbital_tangent_y','orbital_normal_z'],sig)),'total_stream_mass_kg':Ms,'kinetic_energy_J':K,'loss_W_at_0_001_per_revolution':.001*K/period,'stream_mass_delta_from_old_kg':Ms-s['stream_mass_kg'],'stream_mass_delta_kg_m2':Ms/A-s['stream_areal_mass_kg_m2'],'axis_rail_linear_kg_m_at_1km_pitch':dict(zip(['x','y','z'],[z*1000 for z in sig])),'axis_max_normal_force_N_m_at_1km_pitch':dict(zip(['x','y','z'],[z*1000*omega**2*R for z in sig])),'max_fractional_Coriolis_force_away_from_poles':2*math.sqrt(n2)/omega,'ring_count_at_1km_meridional_pitch':3*math.pi*R/1000,'all_family_mass_positive':min(sig)>0},
'breathing':{'incoming_source':'opaque_driver message 2026-10-31, nodal quadrature of archived A/B1/B2/C13 initial profiles; values rounded to8digits','incoming_total_scope':'Launch reactions and capture separated in time; triangle bounds do not assume coincident timing. No launched radiation momentum added.','rows':rows,'free_mass_periodic_range_sum_m':sum(x['constant_mean_support_free_mass_displacement_range_m'] for x in rows),'linear_selfgravity_radial_mode_period_s':2*math.pi/w0,'selfgravity_to_solar_n2_ratio':aself/n2,'directional_tidal_stiffness_screen_s2':dict(zip(['x','y','z'],[aself-12*n2,aself,aself+4*n2])),'angular_average_common_mode_stiffness_s2':aself-8*n2/3,'angular_average_common_mode_period_s':2*math.pi/math.sqrt(aself-8*n2/3),'linear_oscillator_range_sum_bound_m':sum(x['linear_radial_oscillator_range_bound_m'] for x in rows),'scope':'Tide-free common spherical mode; impulses projected on uniform shell mass. Conserved stream angular momentum and linear self gravity; no damping is credited. Solar/shape modes require the next network model. This range is known ephemeris motion, not random differential aiming error.'},
'checks':{'max_sample_vector_force_residual_Pa':max(checks),'energy_equals_trace_work_relative':K/(A*R*R*sigma*sum(D)/6)-1}}
assert max(checks)<1e-14
assert abs(result['checks']['energy_equals_trace_work_relative'])<1e-12
Path(__file__).with_name('support-topology.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
print(json.dumps(result,indent=2))
