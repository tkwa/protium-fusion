import sys,json
from pathlib import Path
import numpy as np
ROOT=Path('/home/tkwa/.codex/worktrees/protium-fusion/website/docs/protium-fusion');sys.path.insert(0,str(ROOT))
import pulse_model as pm
from fast_eos import FastEOS,KB_CGS,NA
p=Path('/tmp/protium-small-stages/scaled-hot-B2-64-R30-tight');d=json.loads(p.with_suffix('.json').read_text());trace=np.load(p.with_suffix('.npz'));z=trace['state'].reshape(64,32,-1);mass=d['mass_kg']*1000*np.diff(np.linspace(0,1,65)**3);eos=FastEOS();y0=np.zeros(27);y0[1]=2/14;y0[3]=1/14;center=eos.state(1e5,3e8,float(y0@pm.ZZ),float(y0.sum()));cs0=np.sqrt(center.gamma1*center.p/1e5)
rows=[]
for k in [-10,-7,-5,-3,-2,-1]:
 yy=z[:,:27,k];r=np.r_[0,z[:,28,k]*3000];v=np.r_[0,z[:,29,k]*cs0];rho=mass/(4*np.pi/3*np.diff(r**3));T=z[:,27,k]*3e8;ycl=np.maximum(yy,0);norm=ycl@pm.AA;ye=ycl@pm.ZZ/norm;yi=ycl.sum(axis=1)/norm;s=eos.state(rho,T,ye,yi);comp=np.maximum(-np.diff(v),0);cs=np.sqrt(s.gamma1*s.p/rho);qv=rho*(2*comp**2+.5*cs*comp);rhodot=-rho**2*4*np.pi*np.diff(r*r*v)/mass
 i=-1;row={'time_s':float(trace['time_s'][k]),'T':T[i],'rho':rho[i],'rhodot':rhodot[i],'expansion_Tdot':float((s.p[i]/rho[i]**2-s.de_drho[i])*rhodot[i]/s.c_v[i]),'viscosity_Tdot':float(qv[i]/rho[i]**2*rhodot[i]/s.c_v[i]),'cv':float(s.c_v[i])}
 for mode in ['raw','floor']:
  tt=max(T[i],1e7) if mode=='floor' else T[i];t=tt/1e9;tf=np.array([1.,1/t,t**(-1/3),t**(1/3),t,t**(5/3),np.log(t)]);rates=np.bincount(pm.RATE_INDEX,weights=np.exp(pm.COEFF@tf),minlength=pm.NR);ypad=np.r_[yy[i],1];flow=rates*pm.PREF*rho[i]**pm.DEXP*np.prod(ypad[pm.REACT],axis=1);dy=pm.STOICH@flow;qtot=pm.Q@flow*pm.MEV/pm.MU*1e4;qnu=pm.QNU@flow*pm.MEV/pm.MU*1e4;duyi=1.5*KB_CGS*NA*T[i];duye=(s.e[i]+rho[i]*s.de_drho[i]-duyi*yi[i])/ye[i];top=np.argsort(abs(pm.Q*flow))[-5:][::-1]
  row[mode]={'nuclear_Tdot':float(qtot/s.c_v[i]),'weak_Tdot':float(-qnu/s.c_v[i]),'composition_Tdot':float((-duyi*dy.sum()-duye*(dy@pm.ZZ))/s.c_v[i]),'top':[(str(pm.RATES[j]),float(flow[j]),float(pm.Q[j])) for j in top]}
 rows.append(row)
Path('/tmp/protium-b2-completion/rhs-audit.json').write_text(json.dumps(rows,indent=2)+'\n')
for row in rows:print(json.dumps(row),flush=True)
