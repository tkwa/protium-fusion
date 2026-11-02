"""Neutral-access and nonneutral-rotation discriminator for warm local ion cells."""
from pathlib import Path
import json,math,hashlib
from scipy.optimize import brentq
P=Path(__file__).resolve().parent
s=json.loads((P/'screen.json').read_text());b=s['baseline'];e=1.602176634e-19;me=9.1093837139e-31;eps=8.8541878188e-12;mu0=1.25663706212e-6;v=b['He_mean_speed_m_s'];k=s['selected_conditional_case']['ionization_coefficient_m3_s'];g=b['geometric_capture_g'];rows=[]
for label,a,L,C,duty in [('sparse_deep_open',.1,.1,1.,1.),('sparse_deep_baffled',.1,.1,.5,1.),('broad_shallow_baffled',1.,.01,.5,1.),('broad_shallow_quarter_duty',1.,.01,.5,.25)]:
 radius=.001;DK=2*radius*v/3;c=C*v/4;Sreq=g*v/(4*a*.5*duty)
 def internal(nu):return math.sqrt(DK*nu)*math.tanh(L*math.sqrt(nu/DK))
 def effective(nu):
  si=internal(nu);return si*c/(si+c)
 nu=brentq(lambda n:effective(n)-Sreq,1e-12,1e10);ne=nu/k;N=ne*a*L
 # Separate cold equilibrium rotation from the~100eV electron energy assumption.
 omega_edge=math.sqrt(2*40*e/me)/radius
 B_rotation=me/e*(omega_edge+ne*e*e/(2*eps*me*omega_edge))
 B_brill=math.sqrt(2*me*ne/eps);B=max(.2,B_brill/math.sqrt(.5),B_rotation)
 wc=e*B/me;nB=eps*B*B/(2*me);w=.5*wc*(1-math.sqrt(1-ne/nB));Erot=.5*me*(w*radius)**2/e
 UB=a*L*B*B/(2*mu0);J=N*1000*e/.1*duty
 rows.append({'case':label,'physical_active_fraction':a,'depth_m':L,'cell_radius_m':radius,'entrance_admission_C':C,'live_duty':duty,'net_transfer_eta':.5,'Knudsen_diffusivity_m2_s':DK,'uniform_electron_ionization_frequency_s':nu,'reaction_diffusion_penetration_m':math.sqrt(DK/nu),'electron_density_m3':ne,'electron_inventory_per_wall_area_m2':N,'entrance_He_density_over_cavity':1-Sreq/c,'volume_mean_He_density_over_cavity':effective(nu)/(nu*L),'net_g_reproduced':4*a*.5*duty*effective(nu)/v,'chosen_B_T':B,'electron_fraction_Brillouin':ne/nB,'edge_rotation_energy_eV':Erot,'uniform_cloud_axis_to_edge_potential_V':ne*e*radius*radius/(4*eps),'gyroradius_at100eV_m':math.sqrt(2*me*100*e)/(e*B),'magnetic_gap_energy_J_m2':UB,'ideal_magnet_mass_kg_m2_at_BH100kJ_m3':2*UB*7500/1e5,'electron_replenishment_W_m2_at1keV_0p1s_no_recovery':J,'electron_loss_time_s_for1W_m2_at1keV':J*.1,'field_scope':'Cold rigid rotation diagnostic; B also chosen to keep edge rotation≤40eV. Does not establish warm cloud stability,100eV distribution, transport or hardware realization.'})
checks=[]
for row in rows:
 checks.append({'case':row['case'],'capture_closure':abs(row['net_g_reproduced']/g-1)<1e-10,'positive_densities':0<row['volume_mean_He_density_over_cavity']<=row['entrance_He_density_over_cavity']<=1,'rotation_cap':row['edge_rotation_energy_eV']<=40.000001,'below_Brillouin':row['electron_fraction_Brillouin']<=.50000001})
for row in checks:assert all(v for k0,v in row.items() if k0!='case'),row
out={'scope':__doc__,'screen_sha256':hashlib.sha256((P/'screen.json').read_bytes()).hexdigest(),'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'model':'Diffuse reflecting cylindrical sidewalls,closed rear,uniform volumetric ionizationnu. D_K=2r*vbar/3; S_internal=sqrt(D_K*nu)tanh[Lsqrt(nu/D_K)]; aperture conductanceC*vbar/4 in series. Netg=4a*eta*duty*S_effective/vbar. Knudsen diffusion is a homogenized long-tube diagnostic, approximate forL/r10. Ion transport and hot capture remain in assumedeta.','reference':'Cylindrical molecular-flow conductance and thermal speed conventions: CERN JUAS Vacuum Technology lecture. No measured complete pump geometry claimed.','rows':rows,'checks':checks,'limitations':'The separate electron cloud radius/length and neutral channel need a shared field/structure design. Return flux extends beyond active electron volume; ideal gap magnet mass is a lower hardware diagnostic. Sidewall material retention and electron-induced gas loads are not solved.'}
(P/'capillary.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(rows,indent=2))
