// Shared fixed-composition EOS adapter for Athena++. Code units documented in prepare.py.
#ifndef PULSE_EOS_HPP
#define PULSE_EOS_HPP
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
namespace pulse {
extern int nr, nt, mode;
extern double gam;
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
inline double lookup(double rho,double T,int k) {
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
  return std::exp(lt[lo]+g*(lt[hi]-lt[lo]));
}
inline double temperature(double rho,double eg) {return invert(rho,eg,1);}
}
#endif
