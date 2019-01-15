#include "InvMultiquadric.hpp"

#define INITVAL_s     DBL_MAX
#define INITVAL_sigma DBL_MAX


//--------------------------------------------------------------------------
InvMultiquadric::
InvMultiquadric() {
  Init();
}


//--------------------------------------------------------------------------
void InvMultiquadric::
Init(void) {
  ReleaseAllMemory();
  s     = INITVAL_s;
  sigma = INITVAL_sigma;
}


//--------------------------------------------------------------------------
void InvMultiquadric::
ReleaseAllMemory(void) {
}


//--------------------------------------------------------------------------
void InvMultiquadric::
Init(double s_, double sigma_) {
  s     = s_;
  sigma = sigma_;
}


//--------------------------------------------------------------------------
InvMultiquadric::
InvMultiquadric(double s_, double sigma_) {
  Init(s_, sigma_);
}


//--------------------------------------------------------------------------
InvMultiquadric::
InvMultiquadric(const InvMultiquadric &G) {
  Init();
  DeepCopy(G);
}


//--------------------------------------------------------------------------
InvMultiquadric& InvMultiquadric::
operator= (const InvMultiquadric &G) {
  if (this != &G) {
    DeepCopy(G);
  }
  return *this;
}


//--------------------------------------------------------------------------
void InvMultiquadric::
DeepCopy(const InvMultiquadric &G) {
  s     = G.s;
  sigma = G.sigma;
}


//--------------------------------------------------------------------------
InvMultiquadric::
~InvMultiquadric() {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
double InvMultiquadric::
GetS(void) const {
  return s;
}


//--------------------------------------------------------------------------
double InvMultiquadric::
GetSigma(void) const {
  return sigma;
}
