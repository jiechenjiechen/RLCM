#include "PCG.hpp"

#define INITVAL_NormB -DBL_MAX
#define INITVAL_mIter 0


//--------------------------------------------------------------------------
PCG::
PCG() {
  mRes = NULL;
  Init();
}


//--------------------------------------------------------------------------
void PCG::
Init(void) {
  ReleaseAllMemory();
  NormB = INITVAL_NormB;
  mIter = INITVAL_mIter;
}


//--------------------------------------------------------------------------
void PCG::
ReleaseAllMemory(void) {
  x.ReleaseAllMemory();
  NormB = INITVAL_NormB;
  mIter = INITVAL_mIter;
  Delete_1D_Array<double>(&mRes);
}


//--------------------------------------------------------------------------
PCG::
PCG(const PCG &G) {
  mRes = NULL;
  Init();
  DeepCopy(G);
}


//--------------------------------------------------------------------------
PCG& PCG::
operator= (const PCG &G) {
  if (this != &G) {
    DeepCopy(G);
  }
  return *this;
}


//--------------------------------------------------------------------------
void PCG::
DeepCopy(const PCG &G) {
  ReleaseAllMemory();
  if (G.mRes) {
    New_1D_Array<double, INTEGER>(&mRes, G.mIter);
    memcpy(mRes, G.mRes, G.mIter*sizeof(double));
  }
  NormB = G.NormB;
  x = G.x;
  mIter = G.mIter;
}


//--------------------------------------------------------------------------
PCG::
~PCG() {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
double PCG::
GetNormRHS(void) const {
  return NormB;
}


//--------------------------------------------------------------------------
void PCG::
GetSolution(DVector &Sol) const {
  Sol = x;
}


//--------------------------------------------------------------------------
const double* PCG::
GetResHistory(INTEGER &Iter) const {
  Iter = mIter;
  return mRes;
}
