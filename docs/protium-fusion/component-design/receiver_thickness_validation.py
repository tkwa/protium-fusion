"""Independent finite-volume check of the two slab kernels and grid refinement.
No hydrodynamics; does not modify the receiver-thickness packet.
"""
from pathlib import Path
import json,hashlib,math
import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.integrate import quad
from scipy.special import erfc
P=Path(__file__).parent;source=P/'receiver_thickness.py'
# Load definitions and the four archived pulse inputs, ending before the family
# evaluation/output. This reuses source geometry, while kernel validation below
# uses an independently assembled finite-volume thermal operator.
s=source.read_text();prefix=s[:s.index('\nbig=[n for n in stages')]
ns={'__file__':str(source)};exec(compile(prefix,str(source),'exec'),ns)
front=ns['front_step'];back=ns['back_pdf'];pb=ns['pulse_bins'];stages=ns['stages'];RHO=ns['RHO']
# Set L=1,k=rho*cp so thermal diffusivity is1; normalized temperature k*T.
times=np.geomspace(.001,10.,100);fvrows=[]
for N in [128,512,1024]:
 dx=1/N
 for bc in ['D','N']:
  diag=np.full(N,-2/dx**2);diag[0]=-1/dx**2;diag[-1]=(-3 if bc=='D' else -1)/dx**2
  lam,V=eigh_tridiagonal(diag,np.full(N-1,1/dx**2));forcing=V[0]/dx
  zero=abs(lam)<1e-7;lam[zero]=0.
  gain=np.empty((N,len(times)));gain[zero]=times;gain[~zero]=-np.expm1(lam[~zero,None]*times)/(-lam[~zero,None])
  # The last minus expression is positive: integral exp(lambda*t)dt.
  temps=V@(forcing[:,None]*gain)
  face=temps[0]+dx/2
  if bc=='N':face-=times
  analytic=front(times,1.,RHO,1.,bc)*RHO
  err=float(max(abs(face/analytic-1)))
  fvrows.append({'cells':N,'boundary':bc,'max_front_step_relative_difference':err})
  if N==1024:assert err<.0003,(bc,err)
# Back-face impulse response must transmit exactly the injected energy.
integral=quad(lambda t:float(back(np.array([t]),1.,RHO,1.)[0]),0.,100.,epsabs=1e-11,epsrel=1e-11,limit=300)[0]
assert abs(integral-1)<1e-9
# Actual slow-wave pulse grid refinement at the most responsive thin receiver.
grid=[];L=.5;k=50.;cp=1800.;fft=ns['fftconvolve']
for name,pulse in stages.items():
 peaks=[]
 for dt in [.5,.25,.125]:
  edges=np.arange(round(20000/dt)+1)*dt;qparts=pb(pulse,edges,True);q=sum(qparts[:3]);kernel=np.diff(front(edges,L,k,cp,'D'))
  rise=fft(q,kernel)[:len(q)];peaks.append(float(max(rise)))
 grid.append({'stage':name,'dt_s':[.5,.25,.125],'slow_front_peak_K':peaks,'finest_over_reference_relative':peaks[-1]/peaks[1]-1})
 assert abs(peaks[-1]/peaks[1]-1)<.001
out={'scope':__doc__,'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'finite_volume_checks':fvrows,'back_kernel_energy_integral':integral,'pulse_grid_refinement':grid,'passed':True}
Path(__file__).with_name('receiver-thickness-validation.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
print(json.dumps(out,indent=2))
