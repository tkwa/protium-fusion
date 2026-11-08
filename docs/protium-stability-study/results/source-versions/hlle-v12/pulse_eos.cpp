// Athena++ general-EOS adapter. Original implementation for this study.
#include "../eos.hpp"
#include "../../parameter_input.hpp"
#include "pulse_eos.hpp"
namespace pulse {
int nr=0,nt=0,mode=0;double gam=5./3.;
std::vector<double>lr,lt,dat;
std::atomic<unsigned long long> face_limited(0);
void load(const std::string &filename) {
  if(nr) return;
  std::ifstream f(filename,std::ios::binary);
  if(!f) throw std::runtime_error("Missing pulse EOS file");
  f.read(reinterpret_cast<char*>(&nr),4);f.read(reinterpret_cast<char*>(&nt),4);
  lr.resize(nr);lt.resize(nt);dat.resize(nr*nt*3);
  f.read(reinterpret_cast<char*>(lr.data()),8*nr);
  f.read(reinterpret_cast<char*>(lt.data()),8*nt);
  f.read(reinterpret_cast<char*>(dat.data()),8*dat.size());
  if(!f) throw std::runtime_error("Truncated pulse EOS file");
}
}
Real EquationOfState::PresFromRhoEg(Real rho,Real egas) {
  if(pulse::mode==1) return (pulse::gam-1)*egas;
  return pulse::lookup(rho,pulse::invert(rho,egas,1),0);
}
Real EquationOfState::EgasFromRhoP(Real rho,Real pres) {
  if(pulse::mode==1) return pres/(pulse::gam-1);
  return pulse::lookup(rho,pulse::invert(rho,pres,0),1);
}
Real EquationOfState::AsqFromRhoP(Real rho,Real pres) {
  if(pulse::mode==1) return pulse::gam*pres/rho;
  return pulse::lookup(rho,pulse::invert(rho,pres,0),2)*pres/rho;
}
void EquationOfState::InitEosConstants(ParameterInput *pin) {
  pulse::mode=pin->GetOrAddInteger("problem","ideal",0);
  pulse::gam=pin->GetOrAddReal("hydro","gamma",5./3.);
  if(!pulse::mode) pulse::load(pin->GetString("problem","eos_path"));
}
