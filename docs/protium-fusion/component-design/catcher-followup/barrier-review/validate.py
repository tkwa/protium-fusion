import json,math
from pathlib import Path
import calculate as c
out=[]
def check(n,v):
 assert v,n
 out.append({'name':n,'passed':True})
for T in [1000,1500,1700]:
 for k,m in c.M.items():
  for P in [.001,1,1000,1e5]:
   check(f'kinetic flux identity {T} {k} {P}',abs(c.wall_flux(P,T,m)/(P/(c.KB*T)*c.speed(T,m)/4)-1)<1e-14)
  for t in [.01,1,100,7698.198541942018]:
   phi=c.hold_phi_max(.001,50e-6,t,T,m)
   remain=c.pressure_decay(1.,50e-6,phi,0,T,m,t)
   check(f'finite hold {T} {k} {t}',abs(remain-.999)<1e-14)
for phi in [1e-17,1e-9,1e-4]:
 beta=3.35e-6;h=50e-6;T=1000;m=c.M['H2'];N0=1e19;P0=N0*c.KB*T/h
 loss=c.fraction_leaked(phi,beta)
 check(f'partition {phi}',abs(loss+beta/(phi+beta)-1)<1e-14)
 # Exact finite pressure-time integral closes all outlets plus remaining inventory.
 k=c.aperture_speed(phi+beta,T,m)/h;t=10/k
 P=c.pressure_decay(P0,h,phi,beta,T,m,t)
 intP=P0*(-math.expm1(-k*t))/k
 emitted=c.aperture_speed(phi+beta,T,m)*intP/(c.KB*T)
 check(f'dynamic inventory {phi}',abs((emitted+P*h/(c.KB*T))/N0-1)<1e-14)
check('serial equal barriers',c.serial_effective(2.,2.)==1.)
check('pumped interspace reduces leak',c.serial_effective(2.,2.,1000)<.004)
for z in [.001,1,10,100,1000,10000]:
 # Independent integration over sqrt(P/P0): integrand2u/(1+z*u).
 n=100000;vals=[2*(i/n)/(1+z*i/n) for i in range(n+1)]
 num=(vals[0]+vals[-1]+4*sum(vals[1:-1:2])+2*sum(vals[2:-1:2]))/(3*n)
 check(f'sqrt law integral {z}',abs(num/c.sqrt_pressure_leak_fraction(z)-1)<1e-8)
snap=json.loads((c.ROOT/'input-snapshot.json').read_text());d=c.evaluate(snap)
check('continuous leakage',d['checks']['continuous_fraction_identity'])
check('mass cap repair area',abs(d['refurbishment_sensitivities']['mouth_caps_kg_m2_per_event']/d['refurbishment_sensitivities']['whole_100nm_face_kg_m2_per_event']-.005)<1e-14)
check('all stage peaks converged',all(x['charged_peak32_to64_relative']<1e-4 for x in snap['stages'].values()))
(c.ROOT/'validation.json').write_text(json.dumps({'scope':'analytical and bookkeeping verification; no material leakage or healing performance validated','all_passed':True,'check_count':len(out),'checks':out},indent=2)+'\n')
print(len(out),'checks passed')
