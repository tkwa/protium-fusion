from pathlib import Path
import sys,json,math
import numpy as np
from scipy.optimize import least_squares
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P));from launcher_screen import Table,MU
n=Table('N_EOS_09-18-20.txt',14);h=Table('H_EOS_09-18-20.txt',1)
z=json.loads((P/'launcher-mixture-entropy.json').read_text())['cases']['AH3'];sh=z['mean_shock'];T=sh['T_K'];rho=sh['rho_kg_m3'];V=17*MU/rho
# Equal totalpressure and additivevolume at fixedT.
def f(x):
 vn,vh=np.exp(x);pn,en=n.state([math.log(vn),math.log(T)]);ph,eh=h.state([math.log(vh),math.log(T)])
 if not np.isfinite(pn+ph+en+eh):return [1e6,1e6]
 return [math.log(pn/ph),(vn+3*vh)/V-1]
res=least_squares(f,np.log([.6*V,.4*V/3]),bounds=([n.bounds[0][0],h.bounds[0][0]],[n.bounds[1][0],h.bounds[1][0]]),xtol=1e-12,ftol=1e-12,gtol=1e-12)
vn,vh=np.exp(res.x);pn,en=n.state([math.log(vn),math.log(T)]);ph,eh=h.state([math.log(vh),math.log(T)]);eex=(14*(en-n.E0)+3*(eh-h.E0))/17
out={'scope':'WarmNH3additivevolumeapproximationfrompureNandHfirstprinciplesEOS, commonT/P. Noexplicitchemicalmixtureinteraction orNH3coldEOS; coldmolecularbindingdifference omitted (tensMJ/kg vsGJ/kg).','T_K':T,'rho_kg_m3':rho,'partial_volume_N_m3_per_atom':vn,'partial_volume_H_m3_per_atom':vh,'pressure_Pa':.5*(pn+ph),'internal_above_elemental_cold_reference_J_kg':eex,'FEOS_pressure_Pa':sh['p'],'FEOS_internal_J_kg':sh['e'],'FEOS_over_mixture_P':sh['p']/(.5*(pn+ph)),'FEOS_over_mixture_E':sh['e']/eex,'root_residual':float(max(abs(np.array(f(res.x))))),'success':bool(res.success)}
print(json.dumps(out,indent=2))
(P/'launcher-ammonia-check.json').write_text(json.dumps(out,indent=2)+'\n')
