#include "IsotropicMatern.hpp"

#define INITVAL_s       DBL_MAX
#define INITVAL_nu      DBL_MAX
#define INITVAL_ell     DBL_MAX
#define INITVAL_k0      DBL_MAX
#define INITVAL_sqrt2nu DBL_MAX


//--------------------------------------------------------------------------
IsotropicMatern::
IsotropicMatern() {
  Init();
}


//--------------------------------------------------------------------------
void IsotropicMatern::
Init(void) {
  ReleaseAllMemory();
  s       = INITVAL_s;
  nu      = INITVAL_nu;
  ell     = INITVAL_ell;
  k0      = INITVAL_k0;
  sqrt2nu = INITVAL_sqrt2nu;
}


//--------------------------------------------------------------------------
void IsotropicMatern::
ReleaseAllMemory(void) {
}


//--------------------------------------------------------------------------
void IsotropicMatern::
Init(double s_, double nu_, double ell_) {
  s       = s_;
  nu      = nu_;
  ell     = ell_;
  k0      = pow(2.0, nu-1) * tgamma(nu);
  sqrt2nu = sqrt(2.0*nu);
}


//--------------------------------------------------------------------------
IsotropicMatern::
IsotropicMatern(double s_, double nu_, double ell_) {
  Init(s_, nu_, ell_);
}


//--------------------------------------------------------------------------
IsotropicMatern::
IsotropicMatern(const IsotropicMatern &G) {
  Init();
  DeepCopy(G);
}


//--------------------------------------------------------------------------
IsotropicMatern& IsotropicMatern::
operator= (const IsotropicMatern &G) {
  if (this != &G) {
    DeepCopy(G);
  }
  return *this;
}


//--------------------------------------------------------------------------
void IsotropicMatern::
DeepCopy(const IsotropicMatern &G) {
  s       = G.s;
  nu      = G.nu;
  ell     = G.ell;
  k0      = G.k0;
  sqrt2nu = G.sqrt2nu;
}


//--------------------------------------------------------------------------
IsotropicMatern::
~IsotropicMatern() {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
double IsotropicMatern::
GetS(void) const {
  return s;
}


//--------------------------------------------------------------------------
double IsotropicMatern::
GetNu(void) const {
  return nu;
}


//--------------------------------------------------------------------------
double IsotropicMatern::
GetEll(void) const {
  return ell;
}
