#include "Lanczos.hpp"

#define INITVAL_NormV -DBL_MAX
#define INITVAL_mIter 0


//--------------------------------------------------------------------------
Lanczos::
Lanczos() {
  mRes = NULL;
  Init();
}


//--------------------------------------------------------------------------
void Lanczos::
Init(void) {
  ReleaseAllMemory();
  NormV = INITVAL_NormV;
  mIter = INITVAL_mIter;
}


//--------------------------------------------------------------------------
void Lanczos::
ReleaseAllMemory(void) {
  NormV = INITVAL_NormV;
  mIter = INITVAL_mIter;
  Delete_1D_Array<double>(&mRes);
  T.ReleaseAllMemory();
  V.ReleaseAllMemory();
  S.ReleaseAllMemory();
  U.ReleaseAllMemory();
}


//--------------------------------------------------------------------------
Lanczos::
Lanczos(const Lanczos &G) {
  mRes = NULL;
  Init();
  DeepCopy(G);
}


//--------------------------------------------------------------------------
Lanczos& Lanczos::
operator= (const Lanczos &G) {
  if (this != &G) {
    DeepCopy(G);
  }
  return *this;
}


//--------------------------------------------------------------------------
void Lanczos::
DeepCopy(const Lanczos &G) {
  ReleaseAllMemory();
  if (G.mRes) {
    New_1D_Array<double, INTEGER>(&mRes, G.mIter);
    memcpy(mRes, G.mRes, G.mIter*sizeof(double));
  }
  NormV = G.NormV;
  mIter = G.mIter;
  T = G.T;
  V = G.V;
  S = G.S;
  U = G.U;
}


//--------------------------------------------------------------------------
Lanczos::
~Lanczos() {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
double Lanczos::
GetNormV(void) const {
  return NormV;
}


//--------------------------------------------------------------------------
INTEGER Lanczos::
GetIter(void) const {
  return mIter;
}


//--------------------------------------------------------------------------
void Lanczos::
GetT(DMatrix &T_) const {
  T_ = T;
}


//--------------------------------------------------------------------------
void Lanczos::
GetV(DMatrix &V_) const {
  V_ = V;
}


//--------------------------------------------------------------------------
void Lanczos::
GetRitzValues(DVector &S_) const {
  S_ = S;
}


//--------------------------------------------------------------------------
void Lanczos::
GetRitzVectors(DMatrix &U_) const {
  U_ = U;
}


//--------------------------------------------------------------------------
const double* Lanczos::
GetRes(void) const {
  return mRes;
}
