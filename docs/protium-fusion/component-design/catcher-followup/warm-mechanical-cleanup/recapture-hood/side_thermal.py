"""Additional uniform axial deposition from a small primary-beam tilt.
This is a full-retained-side-energy sensitivity on one actualC13waveform;
no scattering-loss fraction or material failure temperature is assigned.
"""
from pathlib import Path
import json,hashlib,math
P=Path(__file__).resolve().parent;orig=P/'thermal.py';src=orig.read_text().split('rows=[];profileout=[]',1)[0]
src=src.replace('def solve(q,B,k,N):','def solve(q,B,k,N,side_factor=0.):')
src=src.replace('vq=U[0,:]/(rho*cp*dx);','vq=U[0,:]/(rho*cp*dx)+np.sum(U,axis=0)*side_factor/(rho*cp*H);')
src=src.replace('init=B[0]+q[0]*(H-x)/k;','init=B[0]+q[0]*(H-x)/k+q[0]*side_factor/H*(H*H-x*x)/(2*k);')
src=src.replace('Ein=sum(q[:-1])*dt;','Ein=sum(q[:-1])*dt*(1+side_factor);')
ns={'__file__':str(orig),'__name__':'hood_side_readonly'};exec(compile(src,str(orig),'exec'),ns)
np=ns['np'];h=ns['h'];dt=ns['dt'];H=ns['H'];pr=next(z for z in h['profiles'] if z['stage']=='C13' and z['selection']=='event_1');t=np.arange(-20000.,4000.+dt,dt);absolute=pr['prompt_arrival_s']+t;b,ff,m,de=ns['ns']['evaluate'](absolute);B=b+ns['floor']+ff+ns['offset'];q=m+de
# The centralprompt is retained in both definedtopandtiltedsideinputs.
for tm,name in zip(ns['ns']['et'],ns['ns']['enames']):
 ix=int(np.searchsorted(absolute,tm))
 if 0<=ix<len(q):q[ix]+=ns['ns']['stages'][name]['prompt_fluence']/dt
rows=[];pitch=.01
for thick in [1e-5,1e-4]:
 f=1-(pitch/(pitch+thick))**2
 for theta in [1e-5,1e-4,1e-3]:
  for az in [0.,math.pi/4]:
   dx=H*math.tan(theta)*math.cos(az);dy=H*math.tan(theta)*math.sin(az);intercept=1-(pitch-dx)*(pitch-dy)/(pitch+thick)**2;extra=intercept-f
   # Axisalignedbeamloadsonlyonewallfamily;diagonalsharesequally.
   loaded_fraction=f/2 if az==0 else f;ratio=extra/loaded_fraction
   tip=ns['solve'](q,B,50.,100,ratio);mask=t>=.25;im=np.argmax(np.where(mask,tip,-np.inf));rows.append({'wall_thickness_m':thick,'angle_rad':theta,'azimuth_rad':az,'wholeface_interception':intercept,'extra_side_energy_per_loaded_tip_area_relative_to_tip_energy':ratio,'max_tip_K_conditional_full_side_retention':float(tip[im]),'relative_time_at_peak_s':float(t[im]),'q_at_peak_W_m2':float((m+de)[im])})
out={'scope':__doc__,'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'dependency_sha256':hashlib.sha256(orig.read_bytes()).hexdigest(),'selected_profile':pr,'rows':rows,'checks':ns['checks'],'limits':['Uniformaxialsideheatingcorrespondstoauniformincomingbeamwitha smallconstanttilt; narrowedges,nonuniformbeamprofiles andlateraltransportcanalterlocalhotspots.','Foraxialazimuthonlyonewallfamilytakes thesideload, soitscrosssectionishalfthetotal.','Thefullretentioninput is a deliberatelyboundingenergyallocation; no lateralescape orsecondarycooling is assigned.','Wallshape,thermalbow,anisotropy and radiationexchange betweenfinfacesareseparate.','Do not hold these peaktemperatures for anentirepulse.']};(P/'side-thermal.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(rows,indent=2))
