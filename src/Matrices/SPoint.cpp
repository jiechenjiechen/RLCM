#include "SPoint.hpp"

#define INITVAL_d   0
#define INITVAL_nnz 0


//--------------------------------------------------------------------------
SPoint::
SPoint() {
  idx = NULL;
  x = NULL;
  d = INITVAL_d;
  nnz = INITVAL_nnz;
  Init();
}


//--------------------------------------------------------------------------
void SPoint::
Init(void) {
  Init(0, 0);
}


//--------------------------------------------------------------------------
void SPoint::
Init(INTEGER d_, INTEGER nnz_) {
  if (d_ < nnz_) {
    printf("SPoint::Init. Error: d cannot be smaller than nnz. Reset d to nnz.\n");
    d_ = nnz_;
  }
  if (d != d_ || nnz != nnz_) {
    ReleaseAllMemory();
  }
  if (!idx) {
    New_1D_Array<INTEGER, INTEGER>(&idx, nnz_);
  }
  if (!x) {
    New_1D_Array<double, INTEGER>(&x, nnz_);
  }
  d = d_;
  nnz = nnz_;
}


//--------------------------------------------------------------------------
void SPoint::
ReleaseAllMemory(void) {
  Delete_1D_Array<INTEGER>(&idx);
  Delete_1D_Array<double>(&x);
  d = INITVAL_d;
  nnz = INITVAL_nnz;
}


//--------------------------------------------------------------------------
SPoint::
SPoint(INTEGER d_, INTEGER nnz_) {
  idx = NULL;
  x = NULL;
  d = INITVAL_d;
  nnz = INITVAL_nnz;
  Init(d_, nnz_);
}


//--------------------------------------------------------------------------
SPoint::
SPoint(const SPoint& G) {
  idx = NULL;
  x = NULL;
  d = INITVAL_d;
  nnz = INITVAL_nnz;
  Init();
  DeepCopy(G);
}


//--------------------------------------------------------------------------
SPoint& SPoint::
operator= (const SPoint& G) {
  if (this != &G) {
    DeepCopy(G);
  }
  return *this;
}


//--------------------------------------------------------------------------
void SPoint::
DeepCopy(const SPoint &G) {
  if (d != G.d || nnz != G.nnz) {
    ReleaseAllMemory();
  }
  if (G.idx) {
    if (!idx) {
      New_1D_Array<INTEGER, INTEGER>(&idx, G.nnz);
    }
    memcpy(idx, G.idx, G.nnz*sizeof(INTEGER));
  }
  if (G.x) {
    if (!x) {
      New_1D_Array<double, INTEGER>(&x, G.nnz);
    }
    memcpy(x, G.x, G.nnz*sizeof(double));
  }
  d = G.d;
  nnz = G.nnz;
}


//--------------------------------------------------------------------------
SPoint::
~SPoint() {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
INTEGER SPoint::
GetD(void) const {
  return d;
}


//--------------------------------------------------------------------------
INTEGER SPoint::
GetNNZ(void) const {
  return nnz;
}


//--------------------------------------------------------------------------
INTEGER* SPoint::
GetPointerIdx(void) const {
  return idx;
}


//--------------------------------------------------------------------------
double* SPoint::
GetPointerX(void) const {
  return x;
}


//--------------------------------------------------------------------------
void SPoint::
PrintPoint(const char *name) const {
  printf("%s:\n", name);
  for (INTEGER j = 0; j < nnz; j++) {
    printf("%ld:%g ", (long)idx[j], x[j]);
  }
  printf("\n");
}


//--------------------------------------------------------------------------
double SPoint::
InProd(const DPoint &y) const {
  if (d != y.GetD()) {
    printf("SPoint::InProd. Error: Point dimensions do not match. Return NAN.\n");
    return NAN;
  }
  double *my = y.GetPointer();
  double ret = sp_d_ddot(d, nnz, idx, x, my); // not parallel
  return ret;
}


//--------------------------------------------------------------------------
double SPoint::
InProd(const SPoint &y) const {
  if (d != y.GetD()) {
    printf("SPoint::InProd. Error: Point dimensions do not match. Return NAN.\n");
    return NAN;
  }
  INTEGER nnz_y = y.GetNNZ();
  INTEGER *idx_y = y.GetPointerIdx();
  double *my = y.GetPointerX();
  double ret = sp_s_ddot(d, nnz, idx, x, nnz_y, idx_y, my); // not parallel
  return ret;
}


//--------------------------------------------------------------------------
double SPoint::
Dist1(const DPoint &y) const {
  if (d != y.GetD()) {
    printf("SPoint::Dist1. Error: Point dimensions do not match. Return NAN.\n");
    return NAN;
  }
  double *my = y.GetPointer();
  double ret = sp_d_ddist(d, nnz, idx, x, my, fabs); // not parallel
  return ret;
}


//--------------------------------------------------------------------------
double SPoint::
Dist1(const SPoint &y) const {
  if (d != y.GetD()) {
    printf("SPoint::Dist1. Error: Point dimensions do not match. Return NAN.\n");
    return NAN;
  }
  INTEGER nnz_y = y.GetNNZ();
  INTEGER *idx_y = y.GetPointerIdx();
  double *my = y.GetPointerX();
  double ret = sp_s_ddist(d, nnz, idx, x, nnz_y, idx_y, my, fabs); // not parallel
  return ret;
}


//--------------------------------------------------------------------------
double SPoint::
Dist1(const DPoint &y, double *sigma) const {
  if (d != y.GetD()) {
    printf("SPoint::Dist1. Error: Point dimensions do not match. Return NAN.\n");
    return NAN;
  }
  double *my = y.GetPointer();
  double ret = sp_d_ddists(d, nnz, idx, x, my, sigma, fabs); // not parallel
  return ret;
}


//--------------------------------------------------------------------------
double SPoint::
Dist1(const SPoint &y, double *sigma) const {
  if (d != y.GetD()) {
    printf("SPoint::Dist1. Error: Point dimensions do not match. Return NAN.\n");
    return NAN;
  }
  INTEGER nnz_y = y.GetNNZ();
  INTEGER *idx_y = y.GetPointerIdx();
  double *my = y.GetPointerX();
  double ret = sp_s_ddists(d, nnz, idx, x, nnz_y, idx_y, my, sigma, fabs); // not parallel
  return ret;
}


//--------------------------------------------------------------------------
double SPoint::
Dist2(const DPoint &y) const {
  if (d != y.GetD()) {
    printf("SPoint::Dist2. Error: Point dimensions do not match. Return NAN.\n");
    return NAN;
  }
  double *my = y.GetPointer();
  double ret = sp_d_ddist(d, nnz, idx, x, my, Square); // not parallel
  return ret;
}


//--------------------------------------------------------------------------
double SPoint::
Dist2(const SPoint &y) const {
  if (d != y.GetD()) {
    printf("SPoint::Dist2. Error: Point dimensions do not match. Return NAN.\n");
    return NAN;
  }
  INTEGER nnz_y = y.GetNNZ();
  INTEGER *idx_y = y.GetPointerIdx();
  double *my = y.GetPointerX();
  double ret = sp_s_ddist(d, nnz, idx, x, nnz_y, idx_y, my, Square); // not parallel
  return ret;
}


//--------------------------------------------------------------------------
double SPoint::
Dist2(const DPoint &y, double *sigma) const {
  if (d != y.GetD()) {
    printf("SPoint::Dist2. Error: Point dimensions do not match. Return NAN.\n");
    return NAN;
  }
  double *my = y.GetPointer();
  double ret = sp_d_ddists(d, nnz, idx, x, my, sigma, Square); // not parallel
  return ret;
}


//--------------------------------------------------------------------------
double SPoint::
Dist2(const SPoint &y, double *sigma) const {
  if (d != y.GetD()) {
    printf("SPoint::Dist2. Error: Point dimensions do not match. Return NAN.\n");
    return NAN;
  }
  INTEGER nnz_y = y.GetNNZ();
  INTEGER *idx_y = y.GetPointerIdx();
  double *my = y.GetPointerX();
  double ret = sp_s_ddists(d, nnz, idx, x, nnz_y, idx_y, my, sigma, Square); // not parallel
  return ret;
}


//--------------------------------------------------------------------------
double SPoint::
KernelFuncChi2(const DPoint &y) const {
  if (d != y.GetD()) {
    printf("DPoint::KernelFuncChi2. Error: Point dimensions do not match. Return NAN.\n");
    return NAN;
  }
  double *my = y.GetPointer();
  double ret = sp_d_dchi2(d, nnz, idx, x, my); // not parallel
  return ret;
}


//--------------------------------------------------------------------------
double SPoint::
KernelFuncChi2(const SPoint &y) const {
  if (d != y.GetD()) {
    printf("DPoint::KernelFuncChi2. Error: Point dimensions do not match. Return NAN.\n");
    return NAN;
  }
  INTEGER nnz_y = y.GetNNZ();
  INTEGER *idx_y = y.GetPointerIdx();
  double *my = y.GetPointerX();
  double ret = sp_s_dchi2(d, nnz, idx, x, nnz_y, idx_y, my); // not parallel
  return ret;
}
