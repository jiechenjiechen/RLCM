#include "ProductLaplace.hpp"

#define INITVAL_s     DBL_MAX
#define INITVAL_sigma DBL_MAX


//--------------------------------------------------------------------------
ProductLaplace::
ProductLaplace() {
  Init();
}


//--------------------------------------------------------------------------
void ProductLaplace::
Init(void) {
  ReleaseAllMemory();
  s     = INITVAL_s;
  sigma = INITVAL_sigma;
}


//--------------------------------------------------------------------------
void ProductLaplace::
ReleaseAllMemory(void) {
}


//--------------------------------------------------------------------------
void ProductLaplace::
Init(double s_, double sigma_) {
  s     = s_;
  sigma = sigma_;
}


//--------------------------------------------------------------------------
ProductLaplace::
ProductLaplace(double s_, double sigma_) {
  Init(s_, sigma_);
}


//--------------------------------------------------------------------------
ProductLaplace::
ProductLaplace(const ProductLaplace &G) {
  Init();
  DeepCopy(G);
}


//--------------------------------------------------------------------------
ProductLaplace& ProductLaplace::
operator= (const ProductLaplace &G) {
  if (this != &G) {
    DeepCopy(G);
  }
  return *this;
}


//--------------------------------------------------------------------------
void ProductLaplace::
DeepCopy(const ProductLaplace &G) {
  s     = G.s;
  sigma = G.sigma;
}


//--------------------------------------------------------------------------
ProductLaplace::
~ProductLaplace() {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
double ProductLaplace::
GetS(void) const {
  return s;
}


//--------------------------------------------------------------------------
double ProductLaplace::
GetSigma(void) const {
  return sigma;
}
