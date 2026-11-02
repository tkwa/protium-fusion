"""Conservation, spherical-kernel and discretization checks of the gas packet."""
from pathlib import Path
import json,math,hashlib
import numpy as np
from scipy.integrate import quad
import chamber as c
P=Path(__file__).resolve().parent;d=json.loads((P/'inputs.json').read_text());fine=json.loads((P/'chamber-dt150.0-q128.json').read_text());coarse=json.loads((P/'chamber-dt300.0-q64.json').read_text());pump=json.loads((P/'pumping-interfaces.json').read_text());checks=[]
def ck(name,ok,diagnostic=None):
 checks.append({'name':name,'passed':bool(ok),'diagnostic':diagnostic})
 if not ok:raise AssertionError((name,diagnostic))
A=4*math.pi*d['R_m']**2;V=A*d['R_m']/3;Q=[sum(s[key]/d['plant']['pulse_interval_s'][cat] for cat,s in d['stages'].items()) for key in ['cold_H_atomic_kg','cold_He_atomic_kg']]
ck('actual stage H closes plant H',math.isclose(Q[0],d['plant']['recovered_H_kg_s'],rel_tol=1e-12));ck('actual stage He closes net He',math.isclose(Q[1],d['plant']['helium_per_s']*d['atomic_mass_kg']['he4'],rel_tol=1e-9))
for st in fine['stationary']:
 for i in range(2):
  q=st['Q_atomic_kg_s'][i];m=st['gas_particle_mass_kg'][i];T=st['T_K'];v=st['vbar_m_s'][i];tau=st['unit_capture_residence_s'][i];M=st['unit_capture_inventory_kg'][i];sig=st['unit_capture_columns_kg_m2'][i];p=st['unit_capture_pressures_Pa'][i]
  ck(f'{T}/{m} wall flux',abs(M/V*v/4*A/q-1)<1e-12);ck(f'{T}/{m} column',abs(3*M/A/sig-1)<1e-12);ck(f'{T}/{m} pressure',abs(M/m*c.KB*T/V/p-1)<1e-12)
# Independent quadratures over emitted speed and direction, avoiding closed forms.
for z in [.01,.1,.5,1,2,5]:
 S=quad(lambda w:2*w**3*math.exp(-w*w)*max(0,1-(z*w/2)**2),0,2/z)[0]
 E=quad(lambda w:4/3*w**5*math.exp(-w*w)*max(0,1-(z*w/2)**2),0,2/z)[0]
 ck(f'S kernel z{z}',abs(float(c.firstflight_S(z))-S)<1e-10);ck(f'pressure kernel z{z}',abs(float(c.firstflight_pressure_kernel(z))-E)<1e-10)
for x in [.01,.2,.8,1.2,1.9]:
 I=quad(lambda mu:2*mu/(1+x*x-2*x*mu),x/2,1)[0];ck(f'angular column kernel x{x}',abs(c.firstflight_I(x)-I)<1e-9)
IS=quad(lambda z:float(c.firstflight_S(z)),0,1)[0]+quad(lambda z:float(c.firstflight_S(z)),1,100)[0]
IP=quad(lambda z:float(c.firstflight_pressure_kernel(z)),0,1)[0]+quad(lambda z:float(c.firstflight_pressure_kernel(z)),1,100)[0]
ck('survival steady normalization',abs(IS/(2*math.sqrt(math.pi)/3)-1)<1e-6);ck('pressure steady normalization',abs(IP/(2*math.sqrt(math.pi)/3)-1)<1e-6)
changes=[]
for a,b in zip(coarse['cases'],fine['cases']):
 vals={k:b[k]/a[k]-1 for k in ['mean_radial_column_kg_m2','peak_radial_column_kg_m2','peak_pressure_Pa']};vals['maximum_encounter_relative_change']=max(abs(y['maximum_column_kg_m2']/x['maximum_column_kg_m2']-1) for cat in a['source_profile_pellet_encounters'] for x,y in zip(a['source_profile_pellet_encounters'][cat],b['source_profile_pellet_encounters'][cat]));changes.append(vals);ck('refinement '+str((b['model'],b['T_K'],b['effective_wall_capture_g'],b['material_release_time_s'])),max(abs(x) for x in vals.values())<.001,vals)
for j,sp in enumerate(['p','he4']):
 t,src,_=c.emission(d,sp,150,128);lam=1/fine['stationary'][0]['unit_capture_residence_s'][j];M=c.uniform_response(t,src,Q[j],lam,0);loss=(1-math.exp(-lam*150))*sum(M[:-1]);res=sum(src[1:])-loss-(M[-1]-M[0]);ck(sp+' exact discrete mass balance',abs(res)/sum(src)<1e-11,res/sum(src));ck(sp+' nonnegative inventory',min(M)>0)
 ck(sp+' all arrivals after ignition',all(min(c.quantile_time(s,(np.arange(128)+.5)/128,sp,d['R_m']))>0 for s in d['stages'].values()))
for r in pump['continuous_extraction_scenarios']:
 S=np.array(r['channel_plus_pump_speed']);phi=np.array(r['front_equivalent_open_fraction_for_f_0p001_species']);v=np.sqrt(8*c.KB*1000/(math.pi*np.array([2*d['atomic_mass_kg']['p'],d['atomic_mass_kg']['he4']])));back=S;front=phi*v/4;ck('front competition '+str(r['downstream_pump_speed_m3_s_per_wall_m2_each_species']),np.max(abs(front/(front+back)-.001))<1e-12)
out={'checks_passed':len(checks),'checks':checks,'coarse_to_fine':changes,'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'Numerical identities and refinement only; no pumping capture, surface release population, gas tightness or global cleanup hardware is experimentally qualified.'};(P/'validation.json').write_text(json.dumps(out,indent=2)+'\n');print('passed',len(checks))
