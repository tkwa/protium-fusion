"""Small analytic and provenance checks; no hydrodynamic calculation."""
from pathlib import Path
import hashlib,json,math
import numpy as np
import collector as c
HERE=Path(__file__).resolve().parent
checks=[]
def check(name,error,tolerance):
    checks.append({'name':name,'relative_or_absolute_error':float(error),'tolerance':tolerance,'passed':bool(abs(error)<=tolerance)})
q=np.ones(1000)*2e5;dt=.01;p=c.PROPS
expected=2*q[0]*math.sqrt(10/(math.pi*p['k_W_m_K']*p['rho_kg_m3']*p['cp_J_kg_K']))
t=c.old.temperature_rise(q,dt,p)[-1]
check('Rectangular-flux semi-infinite analytic temperature',t/expected-1,1e-12)
t2=c.old.temperature_rise(q/8,dt*2,p)[-1]
check('Ballistic radius factor2 gives temperature factor2^-2.5',t2/t*2**2.5-1,1e-12)
profile={'velocity_m_s':[10.],'mass_fraction_of_total':[1.]}
check('Uniform dm/dv mass quantile',c.velocity_quantile(profile,1.,.001,'mass')/.01-1,1e-12)
check('Uniform dm/dv kinetic-energy quantile',c.velocity_quantile(profile,1.,.001,'energy')/1.-1,1e-12)
out=json.loads((HERE/'collector.json').read_text());inp=json.loads((HERE/'collector-inputs.json').read_text());proc=json.loads((HERE/'process.json').read_text())
for file,h in out['input_sha256'].items():check('Current input hash '+file,int(hashlib.sha256((HERE/file).read_bytes()).hexdigest()!=h),0)
for file,h in out['code_sha256'].items():check('Current code hash '+file,int(hashlib.sha256((c.ROOT/file).read_bytes()).hexdigest()!=h),0)
for mode,got in out['occupied_s_per_He_at_R1e9'].items():
    expected=sum(out['stages'][k]['timing_at_R1e9'][mode]['occupied_s']*p['pulses_per_He'] for k,p in proc['selected']['stage_pulse_data'].items())
    check('Reconstruct combined cycle occupancy '+mode,got/expected-1,1e-12)
for k,s in out['stages'].items():
    check('Initial EOS reconstruction '+k,s['source_rating']['initial_internal_reconstruction_relative_error'],1e-10)
    check('Ballistic energy integration '+k,s['ballistic']['energy_quadrature_relative_error'],2e-3)
    check('Reduced record matches canonical process '+k,int(s['source_sha256']['record']!=proc['selected']['stage_pulse_data'][k]['sha256']),0)
result={'scope':'Analytic numerical checks and provenance, not independent physical validation of collector assumptions.','all_passed':all(x['passed'] for x in checks),'count':len(checks),'checks':checks,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(HERE/'collector-check.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'all_passed':result['all_passed'],'count':len(checks)}))
assert result['all_passed']
