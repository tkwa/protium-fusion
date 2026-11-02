"""Bounded local warm-ion-pump requirements, not demonstrated pump performance."""
from pathlib import Path
import math,json,hashlib
from scipy.integrate import quad
from scipy.optimize import brentq
P=Path(__file__).resolve().parent
E=1.602176634e-19;KB=1.380649e-23;ME=9.1093837139e-31;EPS=8.8541878188e-12;MU0=1.25663706212e-6;MU=1.66053906660e-27
MHE=4.00260325413*MU;MW=183.84*MU;YEAR=365.25*86400
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
src=Path('/home/tkwa/src/protium-fusion/output/catcher-followup/gas/inputs.json');plant=json.loads(src.read_text())['plant']
R=plant['radius_m'];AREA=plant['area_m2'];V=4*math.pi*R**3/3;T=1000.;f=.001;g=.01
vbar=math.sqrt(8*KB*T/(math.pi*MHE));sv=math.sqrt(KB*T/MHE);tau0=4*R/(3*vbar);Gamma=g/tau0;Q=plant['helium_per_s']*f;density=Q/(Gamma*V);flux=density*vbar/4;net=Q/AREA
baseline=dict(radius_m=R,area_m2=AREA,volume_m3=V,gas_K=T,He_source_nuclei_s=plant['helium_per_s'],residual_fraction=f,residual_He_nuclei_s=Q,residual_He_kg_s=Q*MHE,geometric_capture_g=g,He_mean_speed_m_s=vbar,unremoved_mean_wall_time_s=tau0,net_capture_rate_s=Gamma,mean_residence_s=1/Gamma,steady_n_m3=density,steady_partial_pressure_Pa=density*KB*T,radial_mass_column_kg_m2=density*MHE*R,neutral_wall_flux_m2_s=flux,net_capture_flux_m2_s=net,required_pumping_speed_m3_s=Gamma*V,pumping_speed_per_shell_area_m_s=Gamma*V/AREA,net_capture_current_A_m2=net*E,plant_net_export_W_m2=plant['net_export_W_m2'],plant_total_mass_per_shell_area_kg_m2=plant['mass_kg']/AREA)

def single_survival(beta):
 if beta==0:return 1.
 return quad(lambda x:x*math.exp(-x*x/2-beta/x),0,14,epsabs=1e-11)[0]
def probability(beta):return 1-single_survival(beta)**2

def case(a,L,eta,duty=1.,sigma=3.5e-21,admission=1.):
 target=g/(a*admission*eta*duty)
 if target>=1:return dict(aperture_fraction=a,depth_m=L,net_collection_given_ionization=eta,live_fraction=duty,feasible_ionization_probability=False,required_Pion=target,entrance_admission=admission)
 beta=brentq(lambda z:probability(z)-target,1e-14,100.);Ee=100.;ve=math.sqrt(2*Ee*E/ME);k=sigma*ve;ne=beta*sv/(L*k);column=a*L*ne
 Bmin=math.sqrt(2*ME*ne/EPS);B=max(.2,Bmin/math.sqrt(.5));frac=(Bmin/B)**2;r=.001
 rg=math.sqrt(2*ME*Ee*E)/(E*B);phisc=ne*E*r*r/(4*EPS);uB=a*L*B*B/(2*MU0)
 magmin=2*uB*7500/1e5 # ideal load-line bound at assumed(BH)max100kJ/m3
 ionrate=net/eta;u_electron=column*Ee*E
 renew=[{'effective_replenishment_eV':er,'electron_loss_time_s':te,'net_loss_energy_recovery_fraction':recover,'power_W_m2_time_average':duty*column*er*E/te*(1-recover),'electron_makeup_current_A_m2_time_average':duty*column*E/te} for er in [100.,1000.,10000.] for te in [.001,.01,.1,1.] for recover in [0.,.9]]
 resistive=2*6e-8*B*B*L*a/(MU0*MU0*.5*.001*r)
 return dict(aperture_fraction=a,entrance_admission=admission,depth_m=L,net_collection_given_ionization=eta,live_fraction=duty,feasible_ionization_probability=True,required_Pion=target,dimensionless_ionization_optical_parameter=beta,electron_energy_eV=Ee,sigma_ionization_m2=sigma,electron_speed_m_s=ve,ionization_coefficient_m3_s=k,electron_density_m3=ne,electron_inventory_per_shell_area_m2=column,mean_two_crossing_neutral_dwell_s=4*L/vbar,Brillouin_minimum_B_T=Bmin,chosen_B_T=B,chosen_fraction_Brillouin=frac,cell_radius_m=r,electron_gyroradius_m=rg,cell_radius_to_gyroradius=r/rg,uniform_column_self_potential_V=phisc,electron_kinetic_inventory_J_m2=u_electron,electron_self_field_inventory_J_m2=.25*phisc*E*column,magnetic_gap_energy_J_m2=uB,ideal_permanent_magnet_mass_kg_m2_at_BH100kJ_m3=magmin,assumed_magnet_BHmax_J_m3=1e5,assumed_magnet_density_kg_m3=7500.,isolated_micro_solenoid_resistive_W_m2=resistive,microcoil_scope='Each1mm-radius cell separately wound, winding thickness1mm and fill.5, conductor resistivity6e-8ohmm. Deliberate unfavorable architecture diagnostic, not a bound on shared magnets.',electron_ionizing_collision_time_s=1/(density*k),electron_renewal_sensitivities=renew,ion_formation_flux_m2_s=ionrate,ion_target_current_A_m2=ionrate*E,implantation_power_W_m2={str(ev):ionrate*ev*E for ev in [1000.,10000.,100000.]},ionization_threshold_power_W_m2=ionrate*24.587387*E,check_g_reproduced=a*admission*eta*duty*probability(beta))
rows=[case(a,L,eta,duty) for a in [.03,.1,.3] for L in [.01,.1,1.] for eta in [.35,.5,.8] for duty in [1.,.5,.25]]
selected=case(.1,.1,.5)
traps=[]
for nu in [1e12,1e13]:
 for hold in [1.,100.,1000.,86400.]:
  eb=KB*T/E*math.log(nu*hold/(-math.log(.95)))
  traps.append({'attempt_frequency_s':nu,'hold_s':hold,'barrier_eV_for95pct_survival_at1000K':eb})
regen=[]
for barrier in [3.1,4.6,4.7]:
 for temp in [1000.,1500.,1800.,2000.]:
  rate=1e13*math.exp(-barrier*E/(KB*temp));regen.append({'barrier_eV':barrier,'temperature_K':temp,'attempt_frequency_s':1e13,'mean_first_order_escape_time_s':1/rate,'time95pct_release_s':-math.log(.05)/rate})
coatings=[]
for thick in [1e-8,1e-7,1e-6]:
 for capacity in [1e-6,1e-4,.01]:
  for interval in [1000.,86400.]:
   mass=.1*thick*19300;atoms=mass/MW
   coatings.append({'coating_thickness_m':thick,'coated_area_per_shell_area':.1,'coating_mass_kg_m2':mass,'assumed_effective_single_He_trap_capacity_per_W_atom':capacity,'regeneration_interval_s':interval,'captured_He_per_available_trap':net*interval/(atoms*capacity),'whole_coating_processing_W_m2_at2MJ_kg':mass*2e6/interval,'whole_coating_processing_W_m2_at100MJ_kg':mass*1e8/interval,'sensible_heat_W_m2_at150J_kgK_800K_rise':mass*150*800/interval})
# Charge-neutralization circuit and ballistic extraction at prescribed potential.
ions=[]
for volts in [1000.,10000.,100000.]:
 gap=.01;vel=math.sqrt(2*E*volts/MHE);time=math.sqrt(2*MHE*gap*gap/(E*volts));ions.append({'potential_V':volts,'gap_m':gap,'field_V_m':volts/gap,'He_exit_speed_m_s':vel,'uniform_field_transit_s':time,'He_gyro_angle_during_transit_at_selected_B_rad':E*selected['chosen_B_T']/MHE*time,'max_W_recoil_eV':4*MHE*MW/(MHE+MW)**2*volts})
checks=[]
def ck(name,ok):checks.append(dict(check=name,passed=bool(ok)));assert ok,name
ck('source area',abs(AREA/(4*math.pi*R*R)-1)<1e-12);ck('steady capture closure',abs(flux*g/net-1)<1e-12)
for j,row in enumerate(rows):
 if row['feasible_ionization_probability']:
  ck('capture target '+str(j),abs(row['check_g_reproduced']/g-1)<1e-9);ck('density positive '+str(j),row['electron_density_m3']>0);ck('Brillouin necessary condition '+str(j),row['chosen_fraction_Brillouin']<=.50000001)
ck('quadrature collisionless limit',abs(single_survival(0)-1)<1e-12)
for tr in traps:
 S=math.exp(-tr['attempt_frequency_s']*tr['hold_s']*math.exp(-tr['barrier_eV_for95pct_survival_at1000K']*E/(KB*T)))
 ck('trap survival '+str(tr),abs(S-.95)<1e-12)
result={'scope':__doc__,'input_path':str(src),'input_sha256':sha(src),'code_sha256':sha(__file__),'baseline':baseline,'local_model':'Planar two-crossing electron layer; first incoming and then diffuse1000K reflected normal speed independently follow flux-weighted Maxwellian. Aperture fraction is physical active cross section; entrance_admission separately represents neutral conductance loss through guards/baffles. eta folds ion routing and eventual captured transfer. This is a bounded column requirement, not a solved Penning-cell field or kinetic stability model.','ionization_cross_section_basis':'3.5e-21m2 near100eV from primary NASA1977Table2-1 and NIST helium ionization literature; monochromatic energy is a preparation assumption.','selected_conditional_case':selected,'cases':rows,'baffled_selected_geometry_sensitivities':[case(.1,.1,.5,admission=c) for c in [1.,.5,.25]],'trap_requirements':traps,'single_barrier_release_scenarios':regen,'coating_inventory_sensitivities':coatings,'ion_extraction_scenarios':ions,'historical_replication':{'NASA1977_mass_kg':1.5,'measured_He_speed_m3_s':.0006,'measured_pressure_Pa':1e-6*133.322368,'mass_kg_m2_if_naively_replicated':baseline['pumping_speed_per_shell_area_m_s']*1.5/.0006,'scope':'Dimensional extrapolation of measured ambient apparatus, not a pump prediction at1000K or10^-10Pa.'},'global_uncompensated_electron_cloud':{'electron_density_m3_for_same_capture_at_eta05':Gamma/(.5*selected['ionization_coefficient_m3_s']),'center_to_surface_potential_V_uniform_sphere':Gamma/(.5*selected['ionization_coefficient_m3_s'])*E*R*R/(6*EPS),'scope':'Uniform nonneutral diagnostic; a quasineutral volume plasma changes the physics and can discharge charged source grains. The proposed candidate instead uses local shielded cells.'},'uncertainty':'All eta,duty,magnet energy product,retention,regeneration and recovery values are conditional requirement/sensitivity inputs. No measured overall1000Kpump efficiency or probabilistic credible interval is assigned.','checks':checks}
(P/'screen.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps({'baseline':baseline,'selected':{k:v for k,v in selected.items() if k!='electron_renewal_sensitivities'},'checks':len(checks)},indent=2))
