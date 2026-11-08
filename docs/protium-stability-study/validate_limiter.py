"""Test EOS admissibility and identity of the actual compiled slope limiter."""
import json,subprocess
from pathlib import Path
root=Path(__file__).resolve().parent
base=Path('/var/tmp/tkwa-chatgpt/protium-stability')
source=base/'validation/limiter_probe.cpp'
source.write_text(r'''
#include "pulse_eos.hpp"
#include <iostream>
#include <random>
namespace pulse {int nr,nt,mode;double gam;std::vector<double>lr,lt,dat;}
int main(int argc,char**argv) {
 using namespace pulse;std::ifstream f(argv[1],std::ios::binary);
 f.read((char*)&nr,4);f.read((char*)&nt,4);lr.resize(nr);lt.resize(nt);dat.resize(nr*nt*3);
 f.read((char*)lr.data(),8*nr);f.read((char*)lt.data(),8*nt);f.read((char*)dat.data(),8*dat.size());
 std::mt19937_64 rng(20261031);std::uniform_real_distribution<double> u(0.,1.);
 int limited=0;double min_T=1e99;
 for(int i=0;i<20000;++i) {
   double dc=std::pow(10.,-9+9.5*u(rng)),T=std::exp(lt.front()+(std::log(2e9)-lt.front())*u(rng));
   double pc=lookup(dc,T,0),d=dc*std::exp(4*u(rng)-2),p=pc*std::exp(4*u(rng)-2);
   bool original_valid=p>=minimum_pressure(d);
   double fraction=admissible_face_fraction(d,p,dc,pc);
   if(original_valid && fraction!=1)return 1;
   if(!(fraction>=0 && fraction<=1))return 2;
   if(fraction<1)++limited;
   double rd=dc+fraction*(d-dc),rp=pc+fraction*(p-pc),rt=invert(rd,rp,0);
   if(!(rt>=std::exp(lt.front())*(1-1e-9)))return 3;
   min_T=std::min(min_T,rt);
 }
 std::cout<<"{\"sample_count\":20000,\"limited_count\":"<<limited<<",\"minimum_recovered_T_K\":"<<min_T<<",\"passed\":true}\n";
}
''')
exe=source.with_suffix('')
subprocess.run(['g++','-O2','-std=c++11','-I',str(root),str(source),'-o',str(exe)],check=True)
r=subprocess.run([str(exe),str(base/'inputs-yi2/eos.bin')],capture_output=True,text=True,check=True)
result=json.loads(r.stdout)
(base/'validation/limiter-validation.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
