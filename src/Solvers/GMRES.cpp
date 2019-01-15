#include "GMRES.hpp"

#define INITVAL_NormB -DBL_MAX
#define INITVAL_mIter 0


//--------------------------------------------------------------------------
GMRES::
GMRES() {
  mRes = NULL;
  Init();
}


//--------------------------------------------------------------------------
void GMRES::
Init(void) {
  ReleaseAllMemory();
  NormB = INITVAL_NormB;
  mIter = INITVAL_mIter;
}


//--------------------------------------------------------------------------
void GMRES::
ReleaseAllMemory(void) {
  x.ReleaseAllMemory();
  NormB = INITVAL_NormB;
  mIter = INITVAL_mIter;
  Delete_1D_Array<double>(&mRes);
}


//--------------------------------------------------------------------------
GMRES::
GMRES(const GMRES &G) {
  mRes = NULL;
  Init();
  DeepCopy(G);
}


//--------------------------------------------------------------------------
GMRES& GMRES::
operator= (const GMRES &G) {
  if (this != &G) {
    DeepCopy(G);
  }
  return *this;
}


//--------------------------------------------------------------------------
void GMRES::
DeepCopy(const GMRES &G) {
  ReleaseAllMemory();
  if (G.mRes) {
    New_1D_Array<double, INTEGER>(&mRes, G.mIter);
    memcpy(mRes, G.mRes, G.mIter*sizeof(double));
  }
  NormB = G.NormB;
  x     = G.x;
  mIter = G.mIter;
}


//--------------------------------------------------------------------------
GMRES::
~GMRES() {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
double GMRES::
GetNormRHS(void) const {
  return NormB;
}


//--------------------------------------------------------------------------
void GMRES::
GetSolution(DVector &Sol) const {
  Sol = x;
}


//--------------------------------------------------------------------------
const double* GMRES::
GetResHistory(INTEGER &Iter) const {
  Iter = mIter;
  return mRes;
}
