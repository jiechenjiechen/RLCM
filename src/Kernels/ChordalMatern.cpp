#include "ChordalMatern.hpp"

#define INITVAL_s       DBL_MAX
#define INITVAL_nu      DBL_MAX
#define INITVAL_ell     DBL_MAX
#define INITVAL_k0      DBL_MAX
#define INITVAL_sqrt2nu DBL_MAX


//--------------------------------------------------------------------------
ChordalMatern::
ChordalMatern() {
  Init();
}


//--------------------------------------------------------------------------
void ChordalMatern::
Init(void) {
  ReleaseAllMemory();
  s       = INITVAL_s;
  nu      = INITVAL_nu;
  ell     = INITVAL_ell;
  k0      = INITVAL_k0;
  sqrt2nu = INITVAL_sqrt2nu;
}


//--------------------------------------------------------------------------
void ChordalMatern::
ReleaseAllMemory(void) {
}


//--------------------------------------------------------------------------
void ChordalMatern::
Init(double s_, double nu_, double ell_) {
  s       = s_;
  nu      = nu_;
  ell     = ell_;
  k0      = pow(2.0, nu-1) * tgamma(nu);
  sqrt2nu = sqrt(2.0*nu);
}


//--------------------------------------------------------------------------
ChordalMatern::
ChordalMatern(double s_, double nu_, double ell_) {
  Init(s_, nu_, ell_);
}


//--------------------------------------------------------------------------
ChordalMatern::
ChordalMatern(const ChordalMatern &G) {
  Init();
  DeepCopy(G);
}


//--------------------------------------------------------------------------
ChordalMatern& ChordalMatern::
operator= (const ChordalMatern &G) {
  if (this != &G) {
    DeepCopy(G);
  }
  return *this;
}


//--------------------------------------------------------------------------
void ChordalMatern::
DeepCopy(const ChordalMatern &G) {
  s       = G.s;
  nu      = G.nu;
  ell     = G.ell;
  k0      = G.k0;
  sqrt2nu = G.sqrt2nu;
}


//--------------------------------------------------------------------------
ChordalMatern::
~ChordalMatern() {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
double ChordalMatern::
GetS(void) const {
  return s;
}


//--------------------------------------------------------------------------
double ChordalMatern::
GetNu(void) const {
  return nu;
}


//--------------------------------------------------------------------------
double ChordalMatern::
GetEll(void) const {
  return ell;
}
