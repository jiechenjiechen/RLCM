#include "IsotropicGaussian.hpp"

#define INITVAL_s     DBL_MAX
#define INITVAL_sigma DBL_MAX


//--------------------------------------------------------------------------
IsotropicGaussian::
IsotropicGaussian() {
  Init();
}


//--------------------------------------------------------------------------
void IsotropicGaussian::
Init(void) {
  ReleaseAllMemory();
  s     = INITVAL_s;
  sigma = INITVAL_sigma;
}


//--------------------------------------------------------------------------
void IsotropicGaussian::
ReleaseAllMemory(void) {
}


//--------------------------------------------------------------------------
void IsotropicGaussian::
Init(double s_, double sigma_) {
  s     = s_;
  sigma = sigma_;
}


//--------------------------------------------------------------------------
IsotropicGaussian::
IsotropicGaussian(double s_, double sigma_) {
  Init(s_, sigma_);
}


//--------------------------------------------------------------------------
IsotropicGaussian::
IsotropicGaussian(const IsotropicGaussian &G) {
  Init();
  DeepCopy(G);
}


//--------------------------------------------------------------------------
IsotropicGaussian& IsotropicGaussian::
operator= (const IsotropicGaussian &G) {
  if (this != &G) {
    DeepCopy(G);
  }
  return *this;
}


//--------------------------------------------------------------------------
void IsotropicGaussian::
DeepCopy(const IsotropicGaussian &G) {
  s     = G.s;
  sigma = G.sigma;
}


//--------------------------------------------------------------------------
IsotropicGaussian::
~IsotropicGaussian() {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
double IsotropicGaussian::
GetS(void) const {
  return s;
}


//--------------------------------------------------------------------------
double IsotropicGaussian::
GetSigma(void) const {
  return sigma;
}
