#include "IsotropicLaplace.hpp"

#define INITVAL_s     DBL_MAX
#define INITVAL_sigma DBL_MAX


//--------------------------------------------------------------------------
IsotropicLaplace::
IsotropicLaplace() {
  Init();
}


//--------------------------------------------------------------------------
void IsotropicLaplace::
Init(void) {
  ReleaseAllMemory();
  s     = INITVAL_s;
  sigma = INITVAL_sigma;
}


//--------------------------------------------------------------------------
void IsotropicLaplace::
ReleaseAllMemory(void) {
}


//--------------------------------------------------------------------------
void IsotropicLaplace::
Init(double s_, double sigma_) {
  s     = s_;
  sigma = sigma_;
}


//--------------------------------------------------------------------------
IsotropicLaplace::
IsotropicLaplace(double s_, double sigma_) {
  Init(s_, sigma_);
}


//--------------------------------------------------------------------------
IsotropicLaplace::
IsotropicLaplace(const IsotropicLaplace &G) {
  Init();
  DeepCopy(G);
}


//--------------------------------------------------------------------------
IsotropicLaplace& IsotropicLaplace::
operator= (const IsotropicLaplace &G) {
  if (this != &G) {
    DeepCopy(G);
  }
  return *this;
}


//--------------------------------------------------------------------------
void IsotropicLaplace::
DeepCopy(const IsotropicLaplace &G) {
  s     = G.s;
  sigma = G.sigma;
}


//--------------------------------------------------------------------------
IsotropicLaplace::
~IsotropicLaplace() {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
double IsotropicLaplace::
GetS(void) const {
  return s;
}


//--------------------------------------------------------------------------
double IsotropicLaplace::
GetSigma(void) const {
  return sigma;
}
