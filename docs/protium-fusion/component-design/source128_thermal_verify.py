"""Independent identities for the new source128 finite-backing screen."""
from pathlib import Path
import ast,hashlib,json,math
import numpy as np
from scipy.integrate import quad
P=Path(__file__).parent;s=P/'source128_thin_thermal.py';f=P/'source128-thin-thermal.json';d=json.loads(f.read_text());sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();checks=[]
assert d['code_sha256']==sha(s)
for p,h in d['input_sha256'].items():assert sha(p)==h
q=sum(p['fluence']/p['period']for p in d['pulse_fluences'].values());assert abs(q/d['reference']['all_source_mean_W_m2']-1)<1e-13
# Neumann front-impulse temperature at the back, multiplied by area heat
# capacity, is the CDF of a positive delay kernel. Tau=alpha*t/L^2.
def cdf(tau):
 if tau<=0:return 0.
 if tau<.2:return 2/math.sqrt(math.pi*tau)*sum(math.exp(-(2*j+1)**2/(4*tau))for j in range(20))
 return 1+2*sum((-1)**j*math.exp(-(j*math.pi)**2*tau)for j in range(1,100))
t=np.geomspace(1e-5,10,10000);v=np.array([cdf(x)for x in t]);assert v.min()>=-1e-14 and v.max()<=1+1e-14 and min(np.diff(v))>=-1e-13
lag=quad(lambda x:1-cdf(x),0,20,epsabs=1e-12,epsrel=1e-12)[0];assert abs(lag/(1/6)-1)<1e-11
# Independently integrate the analytic odd-mode moment of the cosine basis.
x=np.linspace(0,1,100001)
for n in[1,2,3,4,5]:
 m=float(np.trapezoid((x-.5)*np.cos(n*math.pi*x),x));exact=((-1)**n-1)/(n*math.pi)**2
 assert abs(m-exact)<2e-11;checks.append({'mode':n,'moment_error':m-exact})
# Re-run only the pure curvature step, no pulse reconstruction or output.
tree=ast.parse(s.read_text());ns={'np':np,'math':math,'RHO':2000.};exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef)and n.name=='curvature_step'],type_ignores=[]),str(s),'exec'),ns)
curv=float(ns['curvature_step'](np.array([1e9]),.25,50.,1800.,5e-6)[0]);assert abs(curv/(5e-6/(2*50))-1)<1e-8
rows=[r for r in d['thermal_rows']if r['k_W_m_K']==50.]
full=[r['isolated_peaks']['C13']['isolated_slow_full_front_rise_K_including_its_bulk_mode']for r in rows];assert max(full)/min(full)-1<1e-10
out={'code_sha256':sha(__file__),'input_sha256':{str(f):sha(f),str(s):sha(s)},'all_source_hashes_match':True,'back_delay_CDF_monotone_and_normalized':True,'back_delay_mean_dimensionless':lag,'expected_back_delay_mean':1/6,'back_temperature_offset_identity':'T_back = positive-delay-average(T_bulk) - qmean*L/(2k), after constant-extraction startup transients','curvature_long_time_per_unit_flux':curv,'expected_curvature_long_time':5e-6/(2*50),'cosine_moment_checks':checks,'C13_isolated_front_rise_across_thickness_K':full,'passed':True,'scope':'Independent analytic identities and provenance; does not validate omitted material physics or an uncomputed dispatch schedule.'}
(P/'source128-thin-thermal-validation.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
