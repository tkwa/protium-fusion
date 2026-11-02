"""Single screened-Coulomb hard backscatter, exact recoil and two-leg CSDA escape.
Normal incidence, amorphous carbon, one hard collision. Not full reflection physics.
"""
from pathlib import Path
import sys,json,math,hashlib,importlib.util
import numpy as np
from numpy.polynomial.legendre import leggauss
P=Path(__file__).resolve().parent;S=P.parent;sys.path.insert(0,str(S));import spectra as sp;import calculate as stopping
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
e=1.602176634e-19;mu=1.66053906660e-27;hb=1.054571817e-34;c=299792458.;alpha=7.2973525643e-3;a0=5.29177210544e-11;ke=1.43996454784255e-9*e
class Model:
 def __init__(self,tab,Z,A,na=64,nx=48,screened=True):
  self.tab=tab;self.Z=Z;self.A=A;self.screened=screened;r=A/12.011;g,w=leggauss(na);self.cm=(g+1)/2*(1-r)-1;self.wcm=w/2*(1-r);self.u=-(self.cm+r)/np.sqrt(1+r*r+2*r*self.cm);self.k=(1+r*r+2*r*self.cm)/(1+r)**2;g,w=leggauss(nx);self.y=(g+1)/2;self.wy=w/2;self.ascr=.88534*a0/(6**(1/3) if Z==1 else Z**.23+6**.23)
 def R(self,E):return np.interp(E,self.tab['E'],self.tab['R'])
 def inv(self,R):return np.where(R>0,np.interp(np.maximum(R,0),self.tab['R'],self.tab['E']),0.)
 def xmax(self,E0,cut=0):
  R0=self.R(E0);lo=np.zeros_like(self.u);hi=np.full_like(lo,R0);rm=self.R(cut) if cut else 0;possible=self.R(self.k*E0)>rm
  for _ in range(44):
   x=(lo+hi)/2;f=self.R(self.k*self.inv(R0-x))-x/self.u-rm;lo=np.where(f>0,x,lo);hi=np.where(f<=0,x,hi)
  return np.where(possible,(lo+hi)/2,0)
 def integ(self,E0,lo,hi):
  x=lo[:,None]+(hi-lo)[:,None]*self.y;Ei=self.inv(self.R(E0)-x);Ec=Ei*12.011/(self.A+12.011)*1e6*e;mr=mu*self.A*12.011/(self.A+12.011);pr=np.sqrt(2*mr*Ec);beta=pr/mr/c;corr=(10*math.pi*self.Z*6*alpha)**.12 if self.Z>=2 else 1
  As=(hb/(2*pr*self.ascr))**2*(1.13+3.76*corr*(alpha*self.Z*6/beta)**2) if self.screened else 0.;ds=(ke*self.Z*6/(4*Ec))**2/(As+(1-self.cm[:,None])/2)**2
  w=(hi-lo)[:,None]/100/(12.011*mu)*self.wy*self.wcm[:,None]*2*math.pi*ds;Eout=self.inv(self.R(Ei*self.k[:,None])-x/self.u[:,None]);return w,Eout,np.broadcast_to(self.u[:,None],Eout.shape),Ei
 def row(self,E0,routing=None):
  xx=self.xmax(E0);w,Eo,u,Ei=self.integ(E0,np.zeros_like(xx),xx);R=float(w.sum());out={'E_MeV':E0,'R':R,'RE':float(np.sum(w*Eo)/E0),'mean_exit_E_MeV':float(np.sum(w*Eo)/R),'mean_exit_cosine':float(np.sum(w*u)/R),'weighted_first_flight_s':float(np.sum(w*2*1.1e9*u/np.sqrt(2*np.maximum(Eo,1e-20)*1e6*e/(self.A*mu)))),'collision_below50keV_u_R':float(w[Ei/self.A<.05].sum()),'collision_above1MeV_R':float(w[Ei>1].sum())}
  for cut in[1e-6,1e-4,.001,.01,.1]:out['Rout_below_'+str(cut)+'MeV']=float(self.integ(E0,self.xmax(E0,cut),xx)[0].sum())
  if routing:
   depth=np.interp(Eo,self.tab['E'],self.tab['projected'])/100/2000*1e6*u;out['one_rehit_absorbed_mobile_front_R']=float(np.sum(w*(1-np.interp(depth,*routing))))
  return out

def main():
 rs,pm,old,decay,weak,hashes=sp.setup();tables=stopping.Tables(old);rf=sp.R/'component-design/catcher-followup/channel-routing/calculate.py';spec=importlib.util.spec_from_file_location('reflection_route',rf);route=importlib.util.module_from_spec(spec);spec.loader.exec_module(route);z,col,t=route.solve(20.,.1,150.,1e5,320,128);routing=(np.r_[0,z,150],np.r_[0,col,col[-1]])
 Eg=np.unique(np.r_[np.geomspace(.001,20.,150),.002,.005,.01,.02,.05,.1,.2,.3,.5,1.,2.,5.,10.]);kernels={};controls={};checks=[]
 for name,Z,A in[('p',1,1),('he4',2,4)]:
  tab=tables.get(Z,A);m=Model(tab,Z,A);kernels[name]=[m.row(float(E),routing)for E in Eg];fine=Model(tab,Z,A,128,96);bare=Model(tab,Z,A,64,48,False);controls[name]=[]
  for E in[.001,.005,.02,.1,.3,1.,3.]:
   lo=m.row(E);hi=fine.row(E);err=abs(lo['R']/hi['R']-1);controls[name].append({'E_MeV':E,'base':lo,'fine':hi,'bare_Rutherford':bare.row(E),'relative_quadrature_error':err});checks.append(err<.002)
 s=json.loads((S/'surface-spectra.json').read_text());rows=[]
 for stage,x in s['stages'].items():
  for name,A in[('p',1),('he4',4)]:
   E=np.array(x['packet_energy_MeV_u'])*A;n=np.array(x['species'][name]['packet_nuclei_m2']);N=n.sum();ok=E>=.001;kernel=kernels[name]
   def avg(key):return float(n[ok]@np.interp(np.log(E[ok]),np.log(Eg),[v[key]for v in kernel])/N)
   row={'stage':stage,'species':name,'nuclei_m2':float(N),'mean_incident_flux_m2_s':float(N/x['period_s']),'unmodeled_incident_below1keV_fraction':float(n[~ok].sum()/N),'incident_KE_J_m2':float(n@E)*1e6*e}
   for key in ['R','weighted_first_flight_s','collision_below50keV_u_R','collision_above1MeV_R','one_rehit_absorbed_mobile_front_R']+[k for k in kernel[0]if k.startswith('Rout')]:row[key]=avg(key)
   re=np.array([v['RE']*v['E_MeV']for v in kernel]);row['reflected_KE_J_m2']=float(n[ok]@np.interp(np.log(E[ok]),np.log(Eg),re))*1e6*e;row['mean_first_reflected_E_MeV']=row['reflected_KE_J_m2']/(N*row['R']*1e6*e);row['mean_first_rehit_s']=row['weighted_first_flight_s']/row['R'];rows.append(row)
 means={}
 for name in['p','he4']:
  sel=[v for v in rows if v['species']==name];tot=sum(v['mean_incident_flux_m2_s']for v in sel);keys=['R','unmodeled_incident_below1keV_fraction','weighted_first_flight_s','collision_below50keV_u_R','collision_above1MeV_R','one_rehit_absorbed_mobile_front_R']+[k for k in sel[0]if k.startswith('Rout')];means[name]={key:sum(v[key]*v['mean_incident_flux_m2_s']for v in sel)/tot for key in keys};means[name]['mean_incident_flux_m2_s']=tot
 out={'scope':__doc__,'code_sha256':sha(__file__),'input_sha256':{str(p):sha(p)for p in[S/'calculate.py',S/'surface-spectra.json',S/'sources/nist-graphite-tables.json',rf]},'kernel':kernels,'controls':controls,'stages':rows,'cadence_weighted':means,'checks':{'number':len(checks),'passed':sum(checks)},'limits':['Normalfirstincidence, onehardcollision, straightmean stopping onbothlegs. No multipleangles,crystal channeling,roughness,pores,surfacebarrier or fields.','Wentzel-Moliere followsGeant4 primarymanual with quoted50keV/u lowerwindow; lower-energy collision contributions are explicitlyreported and arenot validated.','BareRutherford is a larger Coulombsinglecollision comparison, not a universalbound on nuclear elastic interference/resonances or multiple collisions.','Incidentenergybelow1keV excluded and carriedas independent[0,1]response. Nearthermalexit numbers depend on unmodeledatomicsurfacephysics.','One-rehitabsorption assumes nextimpactstops at cosine timesmeanprojectedrange thenidealalready-mobile sidecollection. It is a conditionalrouting sensitivity, not fullsecondreflection.','Closed sphere: firstreflection remainsinventory, followschord2Rcosine and hitsnextwall atsame incidencecosine. Firstreflection probability isnot ultimatewarmgasrelease.']}
 (P/'single-collision.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps({'means':means,'checks':out['checks'],'benchmark':{s:[{'E':v['E_MeV'],'R':v['base']['R'],'Eout':v['base']['mean_exit_E_MeV']}for v in controls[s]]for s in controls}},indent=2));assert all(checks)
if __name__=='__main__':main()
