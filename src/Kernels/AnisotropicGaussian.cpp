#include "AnisotropicGaussian.hpp"

#define INITVAL_s DBL_MAX
#define INITVAL_d 0


//--------------------------------------------------------------------------
AnisotropicGaussian::
AnisotropicGaussian() {
  sigma = NULL;
  Init();
}


//--------------------------------------------------------------------------
void AnisotropicGaussian::
Init(void) {
  ReleaseAllMemory();
  s = INITVAL_s;
  d = INITVAL_d;
}


//--------------------------------------------------------------------------
void AnisotropicGaussian::
ReleaseAllMemory(void) {
  Delete_1D_Array<double>(&sigma);
  d = INITVAL_d;
}


//--------------------------------------------------------------------------
void AnisotropicGaussian::
Init(double s_, INTEGER d_, const double *sigma_) {
  s = s_;
  d = d_;
  New_1D_Array<double, INTEGER>(&sigma, d);
  memcpy(sigma, sigma_, d*sizeof(double));
}


//--------------------------------------------------------------------------
AnisotropicGaussian::
AnisotropicGaussian(double s_, INTEGER d_, const double *sigma_) {
  Init(s_, d_, sigma_);
}


//--------------------------------------------------------------------------
AnisotropicGaussian::
AnisotropicGaussian(const AnisotropicGaussian &G) {
  Init();
  DeepCopy(G);
}


//--------------------------------------------------------------------------
AnisotropicGaussian& AnisotropicGaussian::
operator= (const AnisotropicGaussian &G) {
  if (this != &G) {
    DeepCopy(G);
  }
  return *this;
}


//--------------------------------------------------------------------------
void AnisotropicGaussian::
DeepCopy(const AnisotropicGaussian &G) {
  s = G.s;
  if (d != G.d) {
    ReleaseAllMemory();
  }
  if (G.sigma) {
    if (!sigma) {
      New_1D_Array<double, INTEGER>(&sigma, G.d);
    }
    memcpy(sigma, G.sigma, G.d*sizeof(double));
  }
  d = G.d;
}


//--------------------------------------------------------------------------
AnisotropicGaussian::
~AnisotropicGaussian() {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
double AnisotropicGaussian::
GetS(void) const {
  return s;
}


//--------------------------------------------------------------------------
INTEGER AnisotropicGaussian::
GetD(void) const {
  return d;
}


//--------------------------------------------------------------------------
double* AnisotropicGaussian::
GetSigma(void) const {
  return sigma;
}
