"""First launcher/contact checkpoint. Conserved shock states from primary FPEOS tables.
No hydro, no yield prediction. Run from any directory; NumPy/SciPy required.
"""
from pathlib import Path
import json, re, math, hashlib
import numpy as np
from scipy.interpolate import LinearNDInterpolator
from scipy.optimize import least_squares
P=Path(__file__).resolve().parent;ROOT=P.parent
EPS=8.8541878128e-12;MU=1.66053906660e-27;HA=4.3597447222071e-18;QE=1.602176634e-19;R=1.1e9

def digest(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
class Table:
 def __init__(self,name,A):
  self.path=P/'launcher-data'/name; txt=self.path.read_text();self.A=A;self.m=A*MU
  self.E0=float(re.search(r'E0\[Ha\]=\s*([\-\d.]+)',txt)[1])*HA/self.m
  self.V0=float(re.search(r'V0\[A\^3\]=\s*([\-\d.]+)',txt)[1])*1e-30
  rows=[]
  for line in txt.splitlines():
   if not line.startswith('f='):continue
   a=line.split();rows.append([float(a[7])*1e-30,float(a[9]),float(a[11])*1e9,float(a[14])*HA/self.m])
  self.data=np.array(rows);x=np.log(self.data[:,:2]);self.pressure=LinearNDInterpolator(x,np.log(self.data[:,2]));self.energy=LinearNDInterpolator(x,self.data[:,3]);self.bounds=(x.min(0),x.max(0))
 def state(self,logVT):return float(np.exp(self.pressure(logVT)).item()),float(self.energy(logVT).item())
 def shock(self,q):
  # Two equal streams with relative normal speed2u have deposited energyq=u²/2.
  # Rankine–Hugoniot: e-e0=q=.5*P*(V0-V)/m, neglectinitialP.
  def fun(x):
   p,e=self.state(x);V=np.exp(x[0])
   if not np.isfinite(p+e):return np.array([1e6,1e6])
   return np.array([(e-self.E0-q)/max(q,1e7),(.5*p*(self.V0-V)/self.m-q)/max(q,1e7)])
  candidates=[]
  for T in [1e4,3e4,1e5,3e5,1e6,3e6,1e7]:
   x0=np.log([self.V0/4,T]);x0=np.maximum(self.bounds[0]+1e-8,np.minimum(self.bounds[1]-1e-8,x0))
   fit=least_squares(fun,x0,bounds=self.bounds,xtol=1e-11,ftol=1e-11,gtol=1e-11,max_nfev=200)
   err=float(np.max(np.abs(fun(fit.x))));V,T=np.exp(fit.x);p,e=self.state(fit.x)
   if err<1e-5 and V<self.V0:candidates.append((err,V,T,p,e))
  if not candidates:return {'success':False,'reason':'No physical root in tabulated convex hull; no extrapolation','q_J_kg':q}
  err,V,T,p,e=min(candidates)
  up=math.sqrt(2*q);us=up/(1-V/self.V0)
  return {'success':True,'q_J_kg':q,'rho_kg_m3':self.m/V,'rho_initial_kg_m3':self.m/self.V0,'compression':self.V0/V,'temperature_K':T,'pressure_Pa':p,'particle_speed_m_s':up,'shock_speed_m_s':us,'residual_relative':err,'scope':'Primary-table P/e interpolation and Rankine-Hugoniot root. Pure CH/CH2; chosen carbon isotope changes inertial mass only. No radiation, interface mixing or subsequent isentrope.'}

def main():
 stages={};checks=[]
 for stage,path in [('B2',ROOT/'b2-completion/B2-64-compression-cap.initial.json'),('C13',ROOT/'reassessment-2049/small-stages/scaled-hot-C13-64-R30-initial.json')]:
  d=json.loads(path.read_text());m=np.array(d['mass_fraction_of_total']);v=np.array(d['velocity_m_s']);r=np.array(d['outer_radius_m']);nm=np.r_[.5*(m[:-1]+m[1:]),.5*m[-1]];k=float(.5*nm@(v*v));u=math.sqrt(2*k)
  stages[stage]={'profile_sha256':digest(path),'specific_kinetic_J_kg':k,'rms_speed_m_s':u,'angular_mergers':[]}
  tab=Table('C14H28_EOS_09-18-20.txt' if stage=='B2' else 'C18H18_EOS_09-18-20.txt',14)
  for n in [100,1000,10000,100000,1000000]:
   f=2/n-4/(3*n*n);q=k*f;shock=tab.shock(q)
   if shock['success']:checks.append(shock['residual_relative']<1e-5)
   stages[stage]['angular_mergers'].append({'independent_direction_count':n,'cap_half_angle_deg':math.degrees(math.acos(1-2/n)),'typical_axis_spacing_deg':math.degrees(math.sqrt(4*math.pi/n)),'transverse_energy_fraction':f,'q_transverse_J_kg':q,'shock':shock})
 grains=[]
 u=stages['C13']['rms_speed_m_s'];k=.5*u*u;M=1.130973355e13;duration=700.;rho=1000.;Es=1e9
 for a in [1e-4,1e-3,1e-2]:
  gm=4*math.pi/3*rho*a**3;q=4*math.pi*EPS*a*a*Es;qom=q/gm
  for Ea in [1e6,1e7]:
   acc=qom*Ea;L=k/acc;t=u/acc;coat=3*2000*1e-7/(rho*a)
   grains.append({'radius_m':a,'core_density_kg_m3':rho,'grain_mass_kg':gm,'surface_field_V_m':Es,'surface_pressure_Pa':.5*EPS*Es*Es,'grain_potential_V':Es*a,'charge_mass_C_kg':qom,'acceleration_field_V_m':Ea,'acceleration_length_m':L,'acceleration_time_s':t,'integrated_acceleration_voltage_V':k/qom,'grain_charge_C':q,'grains_per_large_shot':M/gm,'whole_pulse_charge_C':M*qom,'mean_whole_source_current_A':M*qom/duration,'co_moving_electron_neutralization_energy_upper_J_kg':.5*qom*Es*a,'100nm_carbon_coat_mass_per_core_mass':coat,'scope':'Constant surface charge and effective accelerating gradient. Stagedtime-varying cells prevent assigning the integrated voltage to one terminal. Beamloading/focusing/materialsurvival unqualified; coating uses planar thinfilm limit.'})
   checks.extend([abs(qom-3*EPS*Es/(rho*a))<1e-12,abs(.5*u*u-qom*Ea*L)/k<1e-12])
 # Thermal expansion is internal to rigidgrains untilrelease. Molecular/ion beams requirelowtransverseT orlatefocusing.
 beam=[]
 kb=1.380649e-23
 for T in [1,10,100,1000]:
  for A in [1,14,1000,1e6]:
   sigma=math.sqrt(kb*T/(A*MU));beam.append({'particle_mass_u':A,'transverse_temperature_K':T,'rms_transverse_speed_m_s':sigma,'freeflight_rms_spread_m':sigma*R/u})
 # Exactglobalneutralization fraction fromCoulombenergy foruniformsphere; geometrydiagnostic.
 qa=next(g['whole_pulse_charge_C'] for g in grains if g['radius_m']==.001)
 neut=[]
 for rc in [100,1000,10000]:
  uc=3*qa*qa/(20*math.pi*EPS*rc);fraction=math.sqrt(.01*M*k/uc)
  neut.append({'assembly_radius_m':rc,'residual_net_charge_fraction_for_1pct_energy':fraction,'scope':'Uniformsphere charge diagnostic;1pct chosen allowance, not a derived yield tolerance.'})
 columns=[]
 g=next(g for g in grains if g['radius_m']==.001 and g['acceleration_field_V_m']==1e6)
 for n in [1000,10000,100000,1000000]:
  current=g['mean_whole_source_current_A']/n
  for bore in [100.,1000.]:
   density=M/duration/n/(u*math.pi*bore*bore)
   eself=current/(2*math.pi*EPS*u*bore)
   columns.append({'grain_radius_m':.001,'effective_acceleration_field_V_m':1e6,'organizational_complexes':1000,'spatially_distinct_outlets':n,'outlets_per_complex':n/1000,'typical_outlet_spacing_on_receiver_m':R*math.sqrt(4*math.pi/n),'rms_assembly_patch_radius_at_1755m_m':1755*2/math.sqrt(n),'column_beam_radius_m':bore,'mean_charged_current_per_outlet_A':current,'exit_space_charge_edge_field_V_m':eself,'exit_beam_bulk_density_kg_m3':density,'illustrative_focusing_gradient_V_m2':3*eself/bore,'tube_length_m':g['acceleration_length_m'],'tube_side_area_m2':2*math.pi*bore*g['acceleration_length_m'],'all_tube_skin_mass_at_10kg_m2':n*2*math.pi*bore*g['acceleration_length_m']*10,'surface_pressure_skin_fraction_at_1GPa_2000kgm3':3*2000/1000*(.5*EPS*1e18)/(2*1e9),'scope':'Finite spatially distributed columns, not renamed rays from one point. Uniform unneutralized exit beam for charge-field diagnostic. Injector speed and current, alternating-gradient stability, emitter duty, beam halo, cooling and external structures still require design. Tube skin is an area/mass sensitivity, not total accelerator mass.'})
 #1cm/10MVm has sameaccelerationlength but10xlesscharge/current than1mm/1MVm.
 for row in list(columns):
  row=row.copy();row['grain_radius_m']=.01;row['effective_acceleration_field_V_m']=1e7
  for key in ['mean_charged_current_per_outlet_A','exit_space_charge_edge_field_V_m','illustrative_focusing_gradient_V_m2']:row[key]*=.1
  columns.append(row)
 out={'scope':'Launcher architecture screening and nonideal CH/CH2 shock states; no physical driver-to-burn oryield validation.','source_urls':{'FPEOS':'https://militzer.berkeley.edu/FPEOS/','CH':'https://doi.org/10.1103/PhysRevE.96.013204','N_shock':'https://doi.org/10.1103/PhysRevLett.129.015701','dust_accelerator':'https://impact.colorado.edu/dustaccelerator'},'code_sha256':digest(__file__),'data_sha256':{p.name:digest(p) for p in (P/'launcher-data').glob('*EOS*.txt')},'stages':stages,'grain_acceleration':grains,'distributed_column_candidates':columns,'molecular_beam_freeflight':beam,'net_charge_diagnostic':neut,'checks':{'count':len(checks),'passed':sum(checks),'all_passed':all(checks)}}
 (P/'launcher-screen.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps(out,indent=2))
if __name__=='__main__':main()
