"""Independent warm B1 N/H additive-volume FPEOS comparison; no cold hydrazine model."""
from pathlib import Path
import sys,json,math,hashlib
import numpy as np
from scipy.optimize import least_squares
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P.parent));from launcher_screen import Table,MU
n=Table('N_EOS_09-18-20.txt',15);h=Table('H_EOS_09-18-20.txt',1)
src=Path('/home/tkwa/src/protium-fusion/output/b1-source/B1-H2-constant-contact-T0p9.inputs.json');data=json.loads(src.read_text());sh=data['contact'];T=sh['T_K'];rho=sh['rho_kg_m3'];V=17*MU/rho

def f(x):
 vn,vh=np.exp(x);pn,en=n.state([math.log(vn),math.log(T)]);ph,eh=h.state([math.log(vh),math.log(T)])
 if not np.isfinite(pn+ph+en+eh):return [1e6,1e6]
 return [math.log(pn/ph),(vn+2*vh)/V-1]
r=least_squares(f,np.log([.65*V,.35*V/2]),bounds=([n.bounds[0][0],h.bounds[0][0]],[n.bounds[1][0],h.bounds[1][0]]),xtol=1e-12,ftol=1e-12,gtol=1e-12);vn,vh=np.exp(r.x);pn,en=n.state([math.log(vn),math.log(T)]);ph,eh=h.state([math.log(vh),math.log(T)]);ex=(15*(en-n.E0)+2*(eh-h.E0))/17
out=dict(scope='Warm N15:H2 additive-volume approximation from elemental first-principlesEOS at equalpressure/temperature; isotope15inertialmass. No chemicalmixtureinteraction, cold hydrazineEOS or realcontacttrajectory. Cold binding reference difference omitted.',T_K=T,rho_kg_m3=rho,pressure_Pa=.5*(pn+ph),internal_above_elemental_cold_reference_J_kg=ex,FEOS_pressure_Pa=sh['p'],FEOS_internal_J_kg=sh['e'],FEOS_over_primary_mixture_P=sh['p']/(.5*(pn+ph)),FEOS_over_primary_mixture_E=sh['e']/ex,root_residual=float(max(abs(np.array(f(r.x))))),success=bool(r.success),source_inputs_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),primary_table_hashes={name:hashlib.sha256((P.parent/'launcher-data'/name).read_bytes()).hexdigest() for name in ['N_EOS_09-18-20.txt','H_EOS_09-18-20.txt']})
assert out['root_residual']<1e-6
(P/'contact-primary-check.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
