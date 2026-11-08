"""Check the compiled EOS adapter and exported reaction coefficients independently."""
import argparse,json,subprocess,sys
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True)
p.add_argument('--inputs',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
sys.path.insert(0,str(a.root));from fast_eos import FastEOS
import pulse_model as pm
code=r'''
#include "pulse_eos.hpp"
#include <iostream>
#include <iomanip>
namespace pulse {int nr,nt,mode;double gam;std::vector<double>lr,lt,dat;}
int main(int argc,char**argv){
 using namespace pulse;std::ifstream f(argv[1],std::ios::binary);
 f.read((char*)&nr,4);f.read((char*)&nt,4);lr.resize(nr);lt.resize(nt);dat.resize(nr*nt*3);
 f.read((char*)lr.data(),8*nr);f.read((char*)lt.data(),8*nt);f.read((char*)dat.data(),8*dat.size());
 double r,t;std::cout<<std::setprecision(17);
 while(std::cin>>r>>t){double p=lookup(r,t,0),e=lookup(r,t,1),g=lookup(r,t,2);
 std::cout<<p<<" "<<e<<" "<<g<<" "<<invert(r,e,1)<<" "<<invert(r,p,0)<<"\n";}
}
'''
src=a.out/'eos_probe.cpp';src.write_text(code);exe=a.out/'eos_probe'
subprocess.run(['g++','-O2','-std=c++11','-I',str(Path(__file__).parent),str(src),'-o',str(exe)],check=True)
rng=np.random.default_rng(20261030)
rho=10**rng.uniform(-1,8.5,10000);T=10**rng.uniform(np.log10(2e6),np.log10(2e9),10000)
text=''.join(f'{r/1e8:.17g} {t:.17g}\n' for r,t in zip(rho,T))
result=subprocess.run([str(exe),str(a.inputs/'eos.bin')],input=text,text=True,capture_output=True,check=True)
v=np.array([list(map(float,line.split())) for line in result.stdout.splitlines()])
manifest=json.loads((a.inputs/'manifest.json').read_text())
state=FastEOS().state(rho/1000,T,8/15,manifest['yi'])
errors={
 'pressure_max_relative':float(np.max(abs(v[:,0]/(state.p/1e21)-1))),
 'energy_max_relative':float(np.max(abs(v[:,1]/(state.e*rho*1e-24)-1))),
 'gamma1_max_relative':float(np.max(abs(v[:,2]/state.gamma1-1))),
 'T_from_energy_roundtrip_max_relative':float(np.max(abs(v[:,3]/T-1))),
 'T_from_pressure_roundtrip_max_relative':float(np.max(abs(v[:,4]/T-1)))}
rates=(a.inputs/'rates.txt').read_text().split();it=iter(rates);q=float(next(it));coeffs=[]
for _ in range(2):coeffs.append(np.array([[float(next(it)) for _ in range(7)] for _ in range(int(next(it)))]))
tr=np.geomspace(1e7,3e9,200);t9=tr/1e9
tf=np.array([np.ones(200),1/t9,t9**(-1/3),t9**(1/3),t9,t9**(5/3),np.log(t9)])
for name,coeff in zip(['n14(p,g)o15','o15(,p)n14'],coeffs):
 calc=np.exp(coeff@tf).sum(axis=0);source_rate=pm.LIB.get_rate_by_name(name)
 true=np.array([source_rate.eval(t) for t in tr]);m=true>1e-280
 errors[name+'_rate_max_relative']=float(np.max(abs(calc[m]/true[m]-1)))
report=dict(sample_count=len(rho),range_rho_kg_m3=[float(rho.min()),float(rho.max())],
 range_T_K=[float(T.min()),float(T.max())],errors=errors,
 pass_adapter=all(v<1e-3 for v in errors.values()))
(a.out/'adapter-validation.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
assert report['pass_adapter']
