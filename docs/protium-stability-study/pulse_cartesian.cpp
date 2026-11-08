// Cartesian 3D cross-check: shaped Voronoi packets with one bulk direction.
// Athena++ code units: L=30 m, rho=1e8 kg/m3, v=1e6 m/s.
// Scalars: s0 = O15 product mass fraction; s1 = original fuel mass fraction.
#include <algorithm>
#include <array>
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
std::vector<std::array<double,3>> sites;
std::vector<double> packet_times;
bool phase_packet_directions=false;
std::vector<std::vector<double>> forward,reverse;
std::vector<double> phase_times,phase_radii,phase_states;
int phase_nt=0,phase_nr=0;double reference_time=11./30.;
double qheat=0, amplitude=0, norm=1, mean=0, variance=1, ambient=1e-9;
double source_accuracy=.02, transverse_norm=1; int ell=2, burn=1, test=0, perturbation=0, divisions=8, angular_family=0;
// Direction moments for shaped, contiguous meridional sectors. Each sector
// translates along its central direction. Polar caps point along the axis.
void packet_direction(double a,double b,double &vr,double &vt) {
  double step=PI/divisions,cc=0,ss=0;
  for(int n=0;n<=divisions;++n) {
    double c=n*step,lo=std::max(a,std::max(0.,c-step/2)),hi=std::min(b,std::min(PI,c+step/2));
    if(hi<=lo)continue;
    double i1=.5*(SQR(std::sin(hi))-SQR(std::sin(lo)));
    double i2=.5*(hi-lo)-.25*(std::sin(2*hi)-std::sin(2*lo));
    cc+=std::cos(c)*i1+std::sin(c)*i2;
    ss+=std::sin(c)*i1-std::cos(c)*i2;
  }
  double weight=std::cos(a)-std::cos(b);vr=cc/weight;vt=ss/weight;
}
double legendre(int l,double x) {
  if(l==0)return 1.;if(l==1)return x;
  double p=1.,v=x;for(int n=2;n<=l;++n){double w=((2*n-1)*x*v-(n-1)*p)/n;p=v;v=w;}return v;
}
double mode_average(double theta0,double theta1) {
  double x0=std::cos(theta1),x1=std::cos(theta0);
  if(ell==0)return 1.;
  if(angular_family==1 && ell!=1) {
    double hi=-.5*(std::cos((1+ell)*theta1)/(1+ell)+std::cos((1-ell)*theta1)/(1-ell));
    double lo=-.5*(std::cos((1+ell)*theta0)/(1+ell)+std::cos((1-ell)*theta0)/(1-ell));
    return (hi-lo)/(x1-x0);
  }
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
  if(n==9)return static_cast<double>(pulse::face_limited.load(std::memory_order_relaxed));
  if(n==10)return static_cast<double>(pulse::dense_face_limited.load(std::memory_order_relaxed));
  double sum=0;AthenaArray<Real> vol(pmb->ncells1);
  for(int k=pmb->ks;k<=pmb->ke;++k)for(int j=pmb->js;j<=pmb->je;++j) {
    pmb->pcoord->CellVolume(k,j,pmb->is,pmb->ie,vol);
    for(int i=pmb->is;i<=pmb->ie;++i) {
      auto &u=pmb->phydro->u;double d=u(IDN,k,j,i),v=vol(i),fuel=pmb->pscalars->s(1,k,j,i);
      if(n==0)sum+=v*(u(IEN,k,j,i)-qheat*pmb->pscalars->s(0,k,j,i));
      if(n==1)sum+=v*fuel;
      if(n==2)sum+=v*pmb->pscalars->s(0,k,j,i);
      if(n==3)sum+=v*fuel*d;
      if(n==4)sum+=v*fuel*(SQR(pmb->pcoord->x1v(i))+SQR(pmb->pcoord->x2v(j))+SQR(pmb->pcoord->x3v(k)));
      if(n==5) {
        double x=pmb->pcoord->x1v(i),y=pmb->pcoord->x2v(j),z=pmb->pcoord->x3v(k);
        double radial=(x*u(IM1,k,j,i)+y*u(IM2,k,j,i)+z*u(IM3,k,j,i))/std::sqrt(x*x+y*y+z*z);
        sum+=v*.5*(SQR(u(IM1,k,j,i))+SQR(u(IM2,k,j,i))+SQR(u(IM3,k,j,i))-radial*radial)/d;
      }
      if(n==6)sum+=v*fuel*pmb->pcoord->x3v(k);
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
  phase_packet_directions=pin->GetOrAddBoolean("problem","phase_packet_directions",false);
  divisions=pin->GetOrAddInteger("problem","divisions",8);
  angular_family=pin->GetOrAddInteger("problem","angular_family",0);
  reference_time=pin->GetOrAddReal("problem","reference_time",11./30.);
  if(perturbation==5) {
    std::ifstream f(pin->GetString("problem","phase_table"),std::ios::binary);
    f.read(reinterpret_cast<char*>(&phase_nt),4);f.read(reinterpret_cast<char*>(&phase_nr),4);
    if(!f || phase_nt<2 || phase_nr<2)throw std::runtime_error("Invalid phase table");
    phase_times.resize(phase_nt);phase_radii.resize(phase_nr+1);phase_states.resize(phase_nt*phase_nr*5);
    f.read(reinterpret_cast<char*>(phase_times.data()),8*phase_nt);
    f.read(reinterpret_cast<char*>(phase_radii.data()),8*(phase_nr+1));
    f.read(reinterpret_cast<char*>(phase_states.data()),8*phase_states.size());
    if(!f)throw std::runtime_error("Truncated phase table");
  }
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
  if(perturbation==4 || perturbation==5) {
    int count=pin->GetInteger("problem","packet_count");
    if(count<4 || count%2)throw std::runtime_error("Packet count must be even and >=4");
    const double golden=PI*(3-std::sqrt(5.));
    for(int n=0;n<count/2;++n) {
      double z=(n+.5)/(count/2),a=n*golden,r=std::sqrt(1-z*z);
      std::array<double,3> v={r*std::cos(a),r*std::sin(a),z};sites.push_back(v);
      for(auto &x:v)x=-x;sites.push_back(v);
    }
  }
  if(perturbation==5) {
    std::ifstream f(pin->GetString("problem","packet_times"));
    packet_times.resize(sites.size());for(auto &x:packet_times)f>>x;
    if(!f)throw std::runtime_error("Missing packet timing coefficients");
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
  AllocateUserHistoryOutput(11);
  const char* names[11]={"E_minus_nuclear","fuel_mass","product_mass","fuel_rho","fuel_r2","theta_ke","fuel_mode","EOS_outside_fuel_mass","fuel_Tmax","EOS_face_limits","EOS_dense_face_limits"};
  for(int i=0;i<11;++i)EnrollUserHistoryOutput(i,history,names[i],i>=8?UserHistoryOperation::max:UserHistoryOperation::sum);
}
void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  const double nodes[2]={.5-.5/std::sqrt(3.),.5+.5/std::sqrt(3.)};
  for(int k=ks;k<=ke;++k)for(int j=js;j<=je;++j)for(int i=is;i<=ie;++i) {
    double d=0,E=0,fuel=0;std::array<double,3> momentum={0,0,0};
    pscalars->s(0,k,j,i)=0;
    for(double a:nodes)for(double b:nodes)for(double c:nodes) {
      std::array<double,3> x={pcoord->x1f(i)+a*pcoord->dx1f(i),pcoord->x2f(j)+b*pcoord->dx2f(j),pcoord->x3f(k)+c*pcoord->dx3f(k)};
      double r=std::sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
      double rho=ambient,eg=pulse::lookup(ambient,2e6,1),v=0,f=0;
      auto it=std::lower_bound(profile.begin(),profile.end(),r,[](const Shell &s,double z){return s.r1<z;});
      if(it!=profile.end()) {
        rho=it->rho;eg=pulse::lookup(rho,it->T,1);f=rho;
        v=it->v0+(it->v1-it->v0)*(r-it->r0)/(it->r1-it->r0);
      }
      std::array<double,3> direction={x[0]/r,x[1]/r,x[2]/r};
      if((perturbation==4 && f>0) || perturbation==5) {
        double best=-1e99;int index=0;
        for(size_t n=0;n<sites.size();++n) {
          double dot=0;for(int q=0;q<3;++q)dot+=direction[q]*sites[n][q];
          if(dot>best){best=dot;index=n;}
        }
        if(perturbation==4 || phase_packet_directions)direction=sites[index];
        if(perturbation==5) {
          double t=reference_time+amplitude*packet_times[index];
          if(t<phase_times.front() || t>phase_times.back())throw std::runtime_error("3D phase seed outside reference trajectory");
          int it=std::max(0,std::min(phase_nt-2,int(std::upper_bound(phase_times.begin(),phase_times.end(),t)-phase_times.begin())-1));
          int ir=std::max(0,std::min(phase_nr-1,int(std::upper_bound(phase_radii.begin(),phase_radii.end(),r)-phase_radii.begin())-1));
          double fraction=(t-phase_times[it])/(phase_times[it+1]-phase_times[it]);
          double state[5];for(int q=0;q<5;++q)state[q]=(1-fraction)*phase_states[(it*phase_nr+ir)*5+q]+fraction*phase_states[((it+1)*phase_nr+ir)*5+q];
          rho=state[0];v=state[1]/rho;eg=state[2]-.5*rho*v*v;f=state[4];
          pscalars->s(0,k,j,i)+=state[3]/8;
        }
      }
      d+=rho/8;fuel+=f/8;E+=(eg+.5*rho*v*v)/8;
      for(int q=0;q<3;++q)momentum[q]+=rho*v*direction[q]/8;
    }
    phydro->u(IDN,k,j,i)=d;phydro->u(IEN,k,j,i)=E;
    for(int q=0;q<3;++q)phydro->u(IM1+q,k,j,i)=momentum[q];
    pscalars->s(1,k,j,i)=fuel;
  }
}
void MeshBlock::UserWorkInLoop() {
  for(int k=ks;k<=ke;++k)for(int j=js;j<=je;++j)for(int i=is;i<=ie;++i){
    double d=phydro->u(IDN,k,j,i),s=pscalars->s(0,k,j,i),f=pscalars->s(1,k,j,i);
    if(!std::isfinite(d) || !std::isfinite(phydro->u(IEN,k,j,i)) || d<=0 || s < -1e-8*d || s>f+1e-8*d)
      throw std::runtime_error("Invalid hydro state or species bound");
  }
}
