"""Physical invariant checks for the corrected scalar viscosity prescription."""
import hashlib,json
from pathlib import Path
import numpy as np

def work(r,v,rho,cs,mode):
 dr=np.diff(r);dvdt=4*np.pi*np.diff(r*r*v);vol=4*np.pi/3*np.diff(r**3)
 c=np.maximum(-np.diff(v),0)
 if mode=='cap':c=np.minimum(c,np.maximum(-dvdt/vol*dr,0))
 elif mode=='divergence':c=np.maximum(-dvdt/vol*dr/3,0)
 q=rho*(2*c*c+.5*cs*c)
 force=4*np.pi*r[1:]**2*(q-np.r_[q[1:],0])
 return q,-q*dvdt,float(force@v[1:])
checks=[]
def ck(name,yes,details=None):
 checks.append({'check':name,'passed':bool(yes),'details':details});assert yes,name
r=np.array([1.,2.]);v=np.array([4.,3.])
ql,hl,wl=work(r,v,np.ones(1),np.ones(1),'legacy');qn,hn,wn=work(r,v,np.ones(1),np.ones(1),'cap')
ck('Spherical expanding-cell counterexample legacy cools',hl[0]<0,float(hl[0]))
ck('Corrected expanding-cell viscosity vanishes',qn[0]==0 and hn[0]==0)
r=np.linspace(0,4,17);rho=np.ones(16);cs=np.ones(16)
for slope in [-2.,2.]:
 v=slope*r;q,h,w=work(r,v,rho,cs,'cap');ql,hl,wl=work(r,v,rho,cs,'legacy')
 ck('Homologous '+str(slope)+' preserves legacy q',np.allclose(q,ql,atol=1e-14,rtol=1e-14))
 ck('Homologous '+str(slope)+' nonnegative dissipation',np.all(h>=0))
rng=np.random.default_rng(32901)
min_heat=0.;residual=0.
for k in range(200):
 r=np.r_[0.,np.cumsum(rng.uniform(.01,1,32))];v=np.r_[0.,rng.normal(0,2,32)];rho=np.exp(rng.normal(0,1,32));cs=np.exp(rng.normal(0,1,32))
 for mode in ['cap','divergence']:
  q,h,w=work(r,v,rho,cs,mode);min_heat=min(min_heat,float(h.min()));scale=max(abs(w),float(np.abs(h).sum()),1);residual=max(residual,abs(w+h.sum())/scale)
ck('400 random spherical states nonnegative dissipation',min_heat>=0,min_heat)
ck('Force/PdV mechanical exchange cancels',residual<1e-14,residual)
source=Path(__file__).with_name('radial_completion.py')
Path(__file__).with_name('viscosity-checks.json').write_text(json.dumps({'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'checks':checks,'passed':len(checks),'failed':0},indent=2)+'\n')
print('Passed',len(checks),'physical viscosity checks')
