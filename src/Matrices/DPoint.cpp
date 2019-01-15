#include "DPoint.hpp"

#define INITVAL_d 0


//--------------------------------------------------------------------------
DPoint::
DPoint() {
  x = NULL;
  d = INITVAL_d;
  Init();
}


//--------------------------------------------------------------------------
void DPoint::
Init(void) {
  Init(0);
}


//--------------------------------------------------------------------------
void DPoint::
Init(INTEGER d_) {
  if (d != d_) {
    ReleaseAllMemory();
  }
  if (!x) {
    New_1D_Array<double, INTEGER>(&x, d_);
  }
  else {
    memset(x, 0, d_*sizeof(double));
  }
  d = d_;
}


//--------------------------------------------------------------------------
void DPoint::
ReleaseAllMemory(void) {
  Delete_1D_Array<double>(&x);
  d = INITVAL_d;
}


//--------------------------------------------------------------------------
DPoint::
DPoint(INTEGER d_) {
  x = NULL;
  d = INITVAL_d;
  Init(d_);
}


//--------------------------------------------------------------------------
DPoint::
DPoint(const DPoint& G) {
  x = NULL;
  d = INITVAL_d;
  Init();
  DeepCopy(G);
}


//--------------------------------------------------------------------------
DPoint& DPoint::
operator= (const DPoint& G) {
  if (this != &G) {
    DeepCopy(G);
  }
  return *this;
}


//--------------------------------------------------------------------------
void DPoint::
DeepCopy(const DPoint &G) {
  if (d != G.d) {
    ReleaseAllMemory();
  }
  if (G.x) {
    if (!x) {
      New_1D_Array<double, int>(&x, G.d);
    }
    memcpy(x, G.x, G.d*sizeof(double));
  }
  d = G.d;
}


//--------------------------------------------------------------------------
DPoint::
~DPoint() {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
INTEGER DPoint::
GetD(void) const {
  return d;
}


//--------------------------------------------------------------------------
double* DPoint::
GetPointer(void) const {
  return x;
}


//--------------------------------------------------------------------------
void DPoint::
SetPoint(const double *x_, INTEGER d_) {
  Init(d_);
  memcpy(x, x_, d*sizeof(double));
}


//--------------------------------------------------------------------------
void DPoint::
SetStandardNormal(void) {
  StandardNormal(x, d);
}


//--------------------------------------------------------------------------
void DPoint::
PrintPoint(const char *name) const {
  printf("%s:\n", name);
  for (INTEGER i = 0; i < d; i++) {
    printf("%g ", x[i]);
  }
  printf("\n");
}


//--------------------------------------------------------------------------
void DPoint::
Normalize(void) {
  double len = 0.0;
  for (INTEGER i = 0; i < d; i++) {
    len += Square(x[i]);
  }
  len = sqrt(len);
  for (INTEGER i = 0; i < d; i++) {
    x[i] /= len;
  }
}


//--------------------------------------------------------------------------
double DPoint::
InProd(const DPoint &y) const {
  if (d != y.GetD()) {
    printf("DPoint::InProd. Error: Point dimensions do not match. Return NAN.\n");
    return NAN;
  }
  double *my = y.GetPointer();
  double ret = 0.0;
  for (INTEGER i = 0; i < d; i++) {
    ret += x[i]*my[i];
  }
  return ret;
}


//--------------------------------------------------------------------------
double DPoint::
Dist1(const DPoint &y) const {
  if (d != y.GetD()) {
    printf("DPoint::Dist1. Error: Point dimensions do not match. Return NAN.\n");
    return NAN;
  }
  double *my = y.GetPointer();
  double ret = 0.0;
  for (INTEGER i = 0; i < d; i++) {
    ret += fabs(x[i]-my[i]);
  }
  return ret;
}


//--------------------------------------------------------------------------
double DPoint::
Dist1(const DPoint &y, const double *sigma) const {
  if (d != y.GetD()) {
    printf("DPoint::Dist1. Error: Point dimensions do not match. Return NAN.\n");
    return NAN;
  }
  double *my = y.GetPointer();
  double ret = 0.0;
  for (INTEGER i = 0; i < d; i++) {
    ret += fabs((x[i]-my[i])/sigma[i]);
  }
  return ret;
}


//--------------------------------------------------------------------------
double DPoint::
Dist2(const DPoint &y) const {
  if (d != y.GetD()) {
    printf("DPoint::Dist2. Error: Point dimensions do not match. Return NAN.\n");
    return NAN;
  }
  double *my = y.GetPointer();
  double ret = 0.0;
  for (INTEGER i = 0; i < d; i++) {
    ret += Square(x[i]-my[i]);
  }
  return ret;
}


//--------------------------------------------------------------------------
double DPoint::
Dist2(const DPoint &y, const double *sigma) const {
  if (d != y.GetD()) {
    printf("DPoint::Dist2. Error: Point dimensions do not match. Return NAN.\n");
    return NAN;
  }
  double *my = y.GetPointer();
  double ret = 0.0;
  for (INTEGER i = 0; i < d; i++) {
    ret += Square((x[i]-my[i])/sigma[i]);
  }
  return ret;
}


//--------------------------------------------------------------------------
void DPoint::
Subtract(const DPoint &y) {
  if (d != y.GetD()) {
    printf("DPoint::Subtract. Error: Point dimensions do not match. Function call takes no effect.\n");
    return;
  }
  double *my = y.GetPointer();
  for (INTEGER i = 0; i < d; i++) {
    x[i] -= my[i];
  }
}


//--------------------------------------------------------------------------
void DPoint::
Subtract(const DPoint &y, DPoint &z) const {
  if (d != y.GetD()) {
    printf("DPoint::Subtract. Error: Point dimensions do not match. Function call takes no effect.\n");
    return;
  }
  z.Init(d);
  double *mz = z.GetPointer();
  double *my = y.GetPointer();
  for (INTEGER i = 0; i < d; i++) {
    mz[i] = x[i] - my[i];
  }
}


//--------------------------------------------------------------------------
void DPoint::
AverageWith(const DPoint &y) {
  if (d != y.GetD()) {
    printf("DPoint::AverageWith. Error: Point dimensions do not match. Function call takes no effect.\n");
    return;
  }
  double *my = y.GetPointer();
  for (INTEGER i = 0; i < d; i++) {
    x[i] = (x[i] + my[i]) / 2.0;
  }
}


//--------------------------------------------------------------------------
void DPoint::
AverageWith(const DPoint &y, DPoint &z) const {
  if (d != y.GetD()) {
    printf("DPoint::AverageWith. Error: Point dimensions do not match. Function call takes no effect.\n");
    return;
  }
  z.Init(d);
  double *mz = z.GetPointer();
  double *my = y.GetPointer();
  for (INTEGER i = 0; i < d; i++) {
    mz[i] = (x[i] + my[i]) / 2.0;
  }
}


//--------------------------------------------------------------------------
double DPoint::
KernelFuncChi2(const DPoint &y) const {
  if (d != y.GetD()) {
    printf("DPoint::KernelFuncChi2. Error: Point dimensions do not match. Return NAN.\n");
    return NAN;
  }
  double *my = y.GetPointer();
  double ret = 0.0;
  for (INTEGER i = 0; i < d; i++) {
    double den = x[i] + my[i];
    if (den != 0.0) {
      ret += 2.0 * x[i] * my[i] / den;
    }
  }
  return ret;
}
