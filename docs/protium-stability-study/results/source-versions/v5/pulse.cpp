// Axisymmetric prescribed hot inflow, plus spherical verification problems.
// Athena++ code units: L=30 m, rho=1e8 kg/m3, v=1e6 m/s.
// Scalars: s0 = O15 product mass fraction; s1 = original fuel mass fraction.
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../eos/general/pulse_eos.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../scalars/scalars.hpp"

namespace {
struct Shell {double r0,r1,rho,T,v0,v1;};
std::vector<Shell> profile;
std::vector<std::vector<double>> forward,reverse;
double qheat=0, amplitude=0, norm=1, mean=0, variance=1, ambient=1e-9;
double source_accuracy=.02, transverse_norm=1; int ell=2, burn=1, test=0, perturbation=0;
double legendre(int l,double x) {
  if(l==0)return 1.;if(l==1)return x;
  double p=1.,v=x;for(int n=2;n<=l;++n){double w=((2*n-1)*x*v-(n-1)*p)/n;p=v;v=w;}return v;
}
double mode_average(double theta0,double theta1) {
  double x0=std::cos(theta1),x1=std::cos(theta0);
  if(ell==0)return 1.;
  return ((legendre(ell+1,x1)-legendre(ell-1,x1))-
          (legendre(ell+1,x0)-legendre(ell-1,x0)))/((2*ell+1)*(x1-x0));
}
double transverse_mode(double theta) {
  double mu=std::cos(theta);
  return -std::sin(theta)*ell*(mu*legendre(ell,mu)-legendre(ell-1,mu))/(mu*mu-1.);
}
double rate(const std::vector<std::vector<double>> &sets,double T) {
  double t=T/1e9,c=std::cbrt(t),tf[7]={1.,1/t,1/c,c,t,t*c*c,std::log(t)},r=0;
  for(auto &s:sets){double ex=0;for(int n=0;n<7;++n)ex+=s[n]*tf[n];r+=std::exp(ex);}return r;
}
double flow(double rho,double eg,double x,double f) {
  if(f<1e-10)return 0.;
  double T=pulse::temperature(rho,eg);
  // Nuclear rates are used only over the source library's nominal interval.
  if(T<1e7) return 0.;
  if(T>3e9) throw std::runtime_error("Reacting temperature exceeds source EOS/rate range");
  double left=std::max(0.,f-x);
  return 3e-5*(rate(forward,T)*(rho*1e5)*left*left/15.-rate(reverse,T)*x);
}
void nuclear(MeshBlock *pmb,const Real time,const Real dt,
 const AthenaArray<Real>& w,const AthenaArray<Real>& r,const AthenaArray<Real>& bcc,
 AthenaArray<Real>& u,AthenaArray<Real>& s) {
  if(!burn || test)return;
  for(int k=pmb->ks;k<=pmb->ke;++k)for(int j=pmb->js;j<=pmb->je;++j)
  for(int i=pmb->is;i<=pmb->ie;++i) {
    double rho=w(IDN,k,j,i),x=r(0,k,j,i),f=r(1,k,j,i);
    double eg=pmb->peos->EgasFromRhoP(rho,w(IPR,k,j,i));
    double dx=dt*flow(rho,eg,x,f);
    s(0,k,j,i)+=rho*dx;u(IEN,k,j,i)+=rho*qheat*dx;
  }
}
Real source_dt(MeshBlock *pmb) {
  if(!burn || test)return 1e30;
  double result=1e30;
  auto &w=pmb->phydro->w;auto &r=pmb->pscalars->r;
  for(int k=pmb->ks;k<=pmb->ke;++k)for(int j=pmb->js;j<=pmb->je;++j)
  for(int i=pmb->is;i<=pmb->ie;++i) {
    double rho=w(IDN,k,j,i),x=r(0,k,j,i),f=r(1,k,j,i);
    if(f<1e-8)continue;
    double eg=pmb->peos->EgasFromRhoP(rho,w(IPR,k,j,i));
    double T=pulse::temperature(rho,eg);if(T<1e7)continue;
    double rf=rate(forward,T)*rho*1e5/15.*3e-5,rr=rate(reverse,T)*3e-5;
    double d=std::abs(rf*(f-x)*(f-x)-rr*x);
    result=std::min(result,source_accuracy/(2*rf*f+rr+1e-100));
    result=std::min(result,source_accuracy*eg/(rho*qheat*d+1e-100));
  }
  return result;
}
Real history(MeshBlock *pmb,int n) {
  double sum=0;AthenaArray<Real> vol(pmb->ncells1);
  for(int k=pmb->ks;k<=pmb->ke;++k)for(int j=pmb->js;j<=pmb->je;++j) {
    pmb->pcoord->CellVolume(k,j,pmb->is,pmb->ie,vol);
    for(int i=pmb->is;i<=pmb->ie;++i) {
      auto &u=pmb->phydro->u;double d=u(IDN,k,j,i),v=vol(i),fuel=pmb->pscalars->s(1,k,j,i);
      if(n==0)sum+=v*(u(IEN,k,j,i)-qheat*pmb->pscalars->s(0,k,j,i));
      if(n==1)sum+=v*fuel;
      if(n==2)sum+=v*pmb->pscalars->s(0,k,j,i);
      if(n==3)sum+=v*fuel*d;
      if(n==4)sum+=v*fuel*SQR(pmb->pcoord->x1v(i));
      if(n==5)sum+=v*.5*SQR(u(IM2,k,j,i))/d;
      if(n==6)sum+=v*fuel*legendre(ell,std::cos(pmb->pcoord->x2v(j)));
      if(n>=7 && !test && fuel/d>1e-4) {
        double kinetic=.5*(SQR(u(IM1,k,j,i))+SQR(u(IM2,k,j,i))+SQR(u(IM3,k,j,i)))/d;
        double T=pulse::temperature(d,u(IEN,k,j,i)-kinetic);
        if(n==7 && (std::log(T)<pulse::lt.front() || std::log(T)>pulse::lt.back() ||
                     std::log(d)<pulse::lr.front() || std::log(d)>pulse::lr.back()))sum+=v*fuel;
        if(n==8)sum=std::max(sum,T);
      }
    }
  }return sum;
}
}
void Mesh::InitUserMeshData(ParameterInput *pin) {
  burn=pin->GetOrAddInteger("problem","burn",1);
  test=pin->GetOrAddInteger("problem","test",0);
  amplitude=pin->GetOrAddReal("problem","amplitude",0.);
  ell=pin->GetOrAddInteger("problem","ell",2);
  perturbation=pin->GetOrAddInteger("problem","perturbation",0);
  ambient=pin->GetOrAddReal("problem","ambient",1e-9);
  source_accuracy=pin->GetOrAddReal("problem","source_accuracy",.02);
  if(!test) {
    std::ifstream f(pin->GetString("problem","profile"));int n;f>>n;
    if(!f)throw std::runtime_error("Missing profile");
    profile.resize(n);for(auto &s:profile)f>>s.r0>>s.r1>>s.rho>>s.T>>s.v0>>s.v1;
    if(!f)throw std::runtime_error("Truncated profile");
    std::ifstream rates(pin->GetString("problem","rates"));rates>>qheat;
    for(auto *sets:{&forward,&reverse}) {rates>>n;sets->resize(n,std::vector<double>(7));
      for(auto &s:*sets)for(auto &c:s)rates>>c;}
    if(!rates)throw std::runtime_error("Missing nuclear rates");
  }
  int nang=mesh_size.nx2;
  double avg2=0,trans2=0;
  for(int j=0;j<nang;++j){double a=PI*j/nang,b=PI*(j+1)/nang,w=(std::cos(a)-std::cos(b))/2;
    double p=mode_average(a,b);mean+=w*p;avg2+=w*p*p;
    double theta=std::acos(.5*(std::cos(a)+std::cos(b)));
    trans2+=w*SQR(transverse_mode(theta));}
  variance=avg2-mean*mean;
  norm=std::sqrt(1+amplitude*amplitude);
  if(nang==1)variance=1;
  transverse_norm=std::sqrt(trans2);
  EnrollUserExplicitSourceFunction(nuclear);
  EnrollUserTimeStepFunction(source_dt);
  AllocateUserHistoryOutput(9);
  const char* names[9]={"E_minus_nuclear","fuel_mass","product_mass","fuel_rho","fuel_r2","theta_ke","fuel_mode","EOS_outside_fuel_mass","fuel_Tmax"};
  for(int i=0;i<9;++i)EnrollUserHistoryOutput(i,history,names[i],i==8?UserHistoryOperation::max:UserHistoryOperation::sum);
}
void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  for(int k=ks;k<=ke;++k)for(int j=js;j<=je;++j)for(int i=is;i<=ie;++i){
    double d=0,m=0,mt=0,E=0,fuel=0,r0=pcoord->x1f(i),r1=pcoord->x1f(i+1);
    double volume=(r1*r1*r1-r0*r0*r0)/3.;
    double angular=block_size.nx2>1?(mode_average(pcoord->x2f(j),pcoord->x2f(j+1))-mean)/std::sqrt(variance):0.;
    double vf=perturbation==0?(1+amplitude*angular)/norm:(perturbation==2?1/norm:1.);
    double density_factor=perturbation==1?1+amplitude*angular:1.;
    double theta=std::acos(.5*(std::cos(pcoord->x2f(j))+std::cos(pcoord->x2f(j+1))));
    double tf=perturbation==2?amplitude*transverse_mode(theta)/transverse_norm/norm:0.;
    if(density_factor<=0)throw std::runtime_error("Density perturbation is nonpositive");
    if(test==1){ // Spherical Noh, gamma=5/3. Initial near-cold radial inflow.
      d=1;m=-1;E=1e-6/(pulse::gam-1)+.5;fuel=1;
    }else if(test==2){ // Exact uniform-density, parabolic-pressure homologous expansion.
      double r2=.6*(std::pow(r1,5)-std::pow(r0,5))/(std::pow(r1,3)-std::pow(r0,3));
      d=1;m=0;E=std::max(1e-12,1-r2)/(pulse::gam-1);fuel=1;
    }else {
      double phase=perturbation==3?amplitude*angular:0.;
      for(auto &s:profile){
        double slope=(s.v1-s.v0)/(s.r1-s.r0),v0=s.v0-slope*s.r0;
        // Ballistic time translation r_new = r_old + v(r_old)*phase.
        // Integrate in ORIGINAL mass coordinates: exact mass and energy preservation.
        double jac=1+slope*phase;
        if(jac<=0 || s.r1+s.v1*phase<=0)throw std::runtime_error("Timing shift crosses shell trajectories");
        double a=std::max((r0-v0*phase)/jac,s.r0),b=std::min((r1-v0*phase)/jac,s.r1);if(b<=a)continue;
        double I2=(std::pow(b,3)-std::pow(a,3))/3,I3=(std::pow(b,4)-std::pow(a,4))/4,I4=(std::pow(b,5)-std::pow(a,5))/5;
        d+=s.rho*density_factor*I2;fuel+=s.rho*density_factor*I2;
        m+=s.rho*density_factor*vf*(slope*I3+v0*I2);
        mt+=s.rho*density_factor*tf*(slope*I3+v0*I2);
        E+=pulse::lookup(s.rho,s.T,1)*density_factor*I2+
          .5*s.rho*density_factor*(vf*vf+tf*tf)*(slope*slope*I4+2*slope*v0*I3+v0*v0*I2);
      }
      double a=std::max(r0,profile.back().r1+profile.back().v1*phase);
      if(r1>a){double v=(std::pow(r1,3)-std::pow(a,3))/3.;
        d+=ambient*v;E+=pulse::lookup(ambient,2e6,1)*v;}
      d/=volume;m/=volume;mt/=volume;E/=volume;fuel/=volume;
    }
    phydro->u(IDN,k,j,i)=d;phydro->u(IM1,k,j,i)=m;
    phydro->u(IM2,k,j,i)=mt;phydro->u(IM3,k,j,i)=0;phydro->u(IEN,k,j,i)=E;
    pscalars->s(0,k,j,i)=0;pscalars->s(1,k,j,i)=fuel;
  }
}
void MeshBlock::UserWorkInLoop() {
  for(int k=ks;k<=ke;++k)for(int j=js;j<=je;++j)for(int i=is;i<=ie;++i){
    double d=phydro->u(IDN,k,j,i),s=pscalars->s(0,k,j,i),f=pscalars->s(1,k,j,i);
    if(!std::isfinite(d) || !std::isfinite(phydro->u(IEN,k,j,i)) || d<=0 || s < -1e-8*d || s>f+1e-8*d)
      throw std::runtime_error("Invalid hydro state or species bound");
  }
}
