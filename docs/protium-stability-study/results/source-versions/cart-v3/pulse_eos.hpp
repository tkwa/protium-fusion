// Shared fixed-composition EOS adapter for Athena++. Code units documented in prepare.py.
#ifndef PULSE_EOS_HPP
#define PULSE_EOS_HPP
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
namespace pulse {
extern int nr, nt, mode;
extern double gam;
extern std::atomic<unsigned long long> face_limited, dense_face_limited;
extern std::vector<double> lr, lt, dat;
void load(const std::string &filename);
inline double field(int i,int j,int k) {return dat[(i*nt+j)*3+k];}
inline void density(double rho,int &i,double &f) {
  double x=std::log(rho);
  // Linear extrapolation in log density is only used for reconstruction/ambient
  // roundoff; actual fuel state range is independently audited from outputs.
  double u=(x-lr[0])/(lr[1]-lr[0]);
  i=std::max(0,std::min(nr-2,int(std::floor(u))));f=u-i;
}
inline double at_r(int i,double f,int j,int k) {
  return (1-f)*field(i,j,k)+f*field(i+1,j,k);
}
inline double minimum_pressure(double rho) {
  int i;double f;density(rho,i,f);return std::exp(at_r(i,f,0,0));
}
inline double minimum_energy(double rho) {
  int i;double f;density(rho,i,f);return std::exp(at_r(i,f,0,1));
}
inline double admissible_face_fraction(double d,double p,double dc,double pc) {
  if(d>0 && p>=minimum_pressure(d))return 1.;
  double lo=0.,hi=1.;
  if(dc>0 && pc>=minimum_pressure(dc)) {
    for(int n=0;n<24;++n) {
      double f=.5*(lo+hi),dd=dc+f*(d-dc),pp=pc+f*(p-pc);
      if(dd>0 && pp>=minimum_pressure(dd))lo=f;else hi=f;
    }
  }
  return lo;
}
inline double lookup(double rho,double T,int k) {
  if(!(rho>0 && T>0 && std::isfinite(T))) {
    std::fprintf(stderr,"EOS invalid lookup: rho=%.17g T=%.17g field=%d\n",rho,T,k);
    throw std::runtime_error("Invalid EOS lookup temperature");
  }
  int i;double f;density(rho,i,f);
  double u=(std::log(T)-lt[0])/(lt[1]-lt[0]);
  int j=std::max(0,std::min(nt-2,int(std::floor(u))));double g=u-j;
  double v=(1-g)*at_r(i,f,j,k)+g*at_r(i,f,j+1,k);
  return k==2?v:std::exp(v);
}
inline double invert(double rho,double val,int k) {
  if (!(rho>0 && val>0)) {
    std::fprintf(stderr,"EOS invalid input: rho=%.17g value=%.17g field=%d\n",rho,val,k);
    throw std::runtime_error("Nonpositive EOS input");
  }
  int i;double f;density(rho,i,f);
  double v=std::log(val);int lo=0,hi=nt-1;
  while(hi-lo>1) {int m=(lo+hi)/2;if(at_r(i,f,m,k)>v) hi=m;else lo=m;}
  double g=(v-at_r(i,f,lo,k))/(at_r(i,f,hi,k)-at_r(i,f,lo,k));
  double result=std::exp(lt[lo]+g*(lt[hi]-lt[lo]));
  if(!(result>0 && std::isfinite(result))) {
    std::fprintf(stderr,"EOS inverse outside numerical range: rho=%.17g value=%.17g field=%d Tmin_value=%.17g logT=%.17g\n",rho,val,k,std::exp(at_r(i,f,0,k)),lt[lo]+g*(lt[hi]-lt[lo]));
    throw std::runtime_error("EOS inverse outside numerical range");
  }
  return result;
}
inline double temperature(double rho,double eg) {return invert(rho,eg,1);}
}
#endif
