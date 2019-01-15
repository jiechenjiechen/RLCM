#include "DVector.hpp"

#define INITVAL_N 0


//--------------------------------------------------------------------------
DVector::
DVector() {
  a = NULL;
  N = INITVAL_N;
  Init();
}


//--------------------------------------------------------------------------
void DVector::
Init(void) {
  Init(0);
}


//--------------------------------------------------------------------------
void DVector::
Init(INTEGER N_) {
  if (N != N_) {
    ReleaseAllMemory();
  }
  if (!a) {
    New_1D_Array<double, INTEGER>(&a, N_);
  }
  else {
    memset(a, 0, N_*sizeof(double));
  }
  N = N_;
}


//--------------------------------------------------------------------------
void DVector::
ReleaseAllMemory(void) {
  Delete_1D_Array<double>(&a);
  N = INITVAL_N;
}


//--------------------------------------------------------------------------
DVector::
DVector(INTEGER N_) {
  a = NULL;
  N = INITVAL_N;
  Init(N_);
}


//--------------------------------------------------------------------------
DVector::
DVector(const DVector &G) {
  a = NULL;
  N = INITVAL_N;
  Init();
  DeepCopy(G);
}


//--------------------------------------------------------------------------
DVector& DVector::
operator= (const DVector &G) {
  if (this != &G) {
    DeepCopy(G);
  }
  return *this;
}


//--------------------------------------------------------------------------
void DVector::
DeepCopy(const DVector &G) {
  if (N != G.N) {
    ReleaseAllMemory();
  }
  if (G.a) {
    if (!a) {
      New_1D_Array<double, INTEGER>(&a, G.N);
    }
    memcpy(a, G.a, G.N*sizeof(double));
  }
  N = G.N;
}


//--------------------------------------------------------------------------
DVector::
~DVector() {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
INTEGER DVector::
GetN(void) const {
  return N;
}


//--------------------------------------------------------------------------
double DVector::
GetEntry(INTEGER i) const {
  if (i < 0 || i >= N) {
    printf("DVector::GetEntry. Error: Invalid index. Return NAN.\n");
    return NAN;
  }
  return a[i];
}


//--------------------------------------------------------------------------
void DVector::
GetBlock(INTEGER RowStart, INTEGER nRow, DVector &b) const {
  if (RowStart < 0 || RowStart+nRow > N) {
    printf("DVector::GetBlock. Error: Invalid range. Function call takes no effect.\n");
    return;
  }
  b.Init(nRow);
  double *mb = b.GetPointer();
  dcopy_(&nRow, a+RowStart, &ONEi, mb, &ONEi);
}


//--------------------------------------------------------------------------
void DVector::
GetBlock(INTEGER RowStart, INTEGER nRow, double *b) const {
  if (RowStart < 0 || RowStart+nRow > N) {
    printf("DVector::GetBlock. Error: Invalid range. Function call takes no effect.\n");
    return;
  }
  if (b == NULL) {
    printf("DVector::GetBlock. Error: double *b is an empty pointer. Function call takes no effect.\n");
    return;
  }
  dcopy_(&nRow, a+RowStart, &ONEi, b, &ONEi);
}


//--------------------------------------------------------------------------
void DVector::
GetBlock(INTEGER *idx, INTEGER n, DVector &b) const {
  b.Init(n);
  double *mb = b.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < n; i++) {
    mb[i] = a[idx[i]];
  }
}


//--------------------------------------------------------------------------
double* DVector::
GetPointer(void) const {
  return a;
}


//--------------------------------------------------------------------------
void DVector::
SetEntry(INTEGER i, double b) {
  if (i < 0 || i >= N) {
    printf("DVector::SetEntry. Error: Invalid index. Function call takes no effect.\n");
    return;
  }
  a[i] = b;
}


//--------------------------------------------------------------------------
void DVector::
SetBlock(INTEGER RowStart, INTEGER nRow, const DVector &b) {
  if (RowStart < 0 || RowStart+nRow > N) {
    printf("DVector::SetBlock. Error: Invalid range. Function call takes no effect.\n");
    return;
  }
  if (nRow != b.GetN()) {
    printf("DVector::SetBlock. Error: Vector dimension not compatible. Function call takes no effect.\n");
    return;
  }
  double *mb = b.GetPointer();
  dcopy_(&nRow, mb, &ONEi, a+RowStart, &ONEi);
}


//--------------------------------------------------------------------------
void DVector::
SetBlock(INTEGER RowStart, INTEGER nRow, const double *b) {
  if (RowStart < 0 || RowStart+nRow > N) {
    printf("DVector::SetBlock. Error: Invalid range. Function call takes no effect.\n");
    return;
  }
  if (b == NULL) {
    printf("DVector::SetBlock. Error: double *b is an empty pointer. Function call takes no effect.\n");
    return;
  }
  dcopy_(&nRow, b, &ONEi, a+RowStart, &ONEi);
}


//--------------------------------------------------------------------------
void DVector::
SetConstVal(double c) {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < N; i++) {
    a[i] = c;
  }
}


//--------------------------------------------------------------------------
void DVector::
SetUniformRandom01(void) {
  UniformRandom01(a, N);
}


//--------------------------------------------------------------------------
void DVector::
SetStandardNormal(void) {
  StandardNormal(a, N);
}


//--------------------------------------------------------------------------
void DVector::
SetStudentT1(void) {
  StudentT1(a, N);
}


//--------------------------------------------------------------------------
void DVector::
SetMultivariateStudentT1(void) {
  MultivariateStudentT1(a, N);
}


//--------------------------------------------------------------------------
void DVector::
SetRandomSech(void) {
  RandomSech(a, N);
}


//--------------------------------------------------------------------------
void DVector::
Permute(const INTEGER *Perm, INTEGER N_) {
  if (N != N_) {
    printf("DVector::Permute. Error: Permutation vector does not have the right length. Function call takes no effect.\n");
    return;
  }
  if (Perm == NULL) {
    printf("DVector::Permute. Error: Perm is an empty pointer. Function call takes no effect.\n");
    return;
  }
  double *b = NULL;
  New_1D_Array<double, INTEGER>(&b, N);
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < N; i++) {
    b[i] = a[Perm[i]];
  }
  dcopy_(&N, b, &ONEi, a, &ONEi);
  Delete_1D_Array<double>(&b);
}


//--------------------------------------------------------------------------
void DVector::
iPermute(const INTEGER *iPerm, INTEGER N_) {
  if (N != N_) {
    printf("DVector::iPermute. Error: Permutation vector does not have the right length. Function call takes no effect.\n");
    return;
  }
  if (iPerm == NULL) {
    printf("DVector::iPermute. Error: iPerm is an empty pointer. Function call takes no effect.\n");
    return;
  }
  double *b = NULL;
  New_1D_Array<double, INTEGER>(&b, N);
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < N; i++) {
    b[iPerm[i]] = a[i];
  }
  dcopy_(&N, b, &ONEi, a, &ONEi);
  Delete_1D_Array<double>(&b);
}


//--------------------------------------------------------------------------
void DVector::
Sort(SortType type) {
  switch (type) {
  case ASCEND:
    qsort(a, N, sizeof(double), CompareNaturalOrderLess);
    break;
  case DESCEND:
    qsort(a, N, sizeof(double), CompareNaturalOrderGreater);
    break;
  }
}


//--------------------------------------------------------------------------
void DVector::
SortByMagnitude(SortType type) {
  switch (type) {
  case ASCEND:
    qsort(a, N, sizeof(double), CompareByMagnitudeLess);
    break;
  case DESCEND:
    qsort(a, N, sizeof(double), CompareByMagnitudeGreater);
    break;
  }
}


//--------------------------------------------------------------------------
INTEGER DVector::
FindLargerThan(double tol, INTEGER *idx) const {
  if (idx == NULL) {
    printf("DVector::FindLargerThan. Error: idx is an empty pointer. Function call takes no effect.\n");
    return 0;
  }
  INTEGER num = 0;
  INTEGER *idx_ptr = idx;
  for (INTEGER i = 0; i < N; i++) {
    if (a[i] > tol) {
      *idx_ptr++ = i;
      num++;
    }
  }
  return num;
}


//--------------------------------------------------------------------------
void DVector::
PrintVectorMatlabForm(const char *name) const {
  printf("%s = [", name);
  for (INTEGER i = 0; i < N; i++) {
    if (i != N-1) {
      printf("%g;\n", a[i]);
    }
    else {
      printf("%g];\n", a[i]);
    }
  }
  printf("\n");
}


//--------------------------------------------------------------------------
void DVector::
Negate(void) {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < N; i++) {
    a[i] = -a[i];
  }
}


//--------------------------------------------------------------------------
void DVector::
Negate(DVector &b) const {
  b.Init(N);
  double *mb = b.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < N; i++) {
    mb[i] = -a[i];
  }
}


//--------------------------------------------------------------------------
void DVector::
Add(double b) {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < N; i++) {
    a[i] += b;
  }
}


//--------------------------------------------------------------------------
void DVector::
Add(const DVector &b) {
  if (N != b.GetN()) {
    printf("DVector::Add. Error: Vector length does not match. Function call takes no effect.\n");
    return;
  }
  double *mb = b.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < N; i++) {
    a[i] += mb[i];
  }
}


//--------------------------------------------------------------------------
void DVector::
Add(double b, DVector &c) const {
  c.Init(N);
  double *mc = c.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < N; i++) {
    mc[i] = a[i] + b;
  }
}


//--------------------------------------------------------------------------
void DVector::
Add(const DVector &b, DVector &c) const {
  if (N != b.GetN()) {
    printf("DVector::Add. Error: Vector length does not match. Function call takes no effect.\n");
    return;
  }
  c.Init(N);
  double *mc = c.GetPointer();
  double *mb = b.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < N; i++) {
    mc[i] = a[i] + mb[i];
  }
}


//--------------------------------------------------------------------------
void DVector::
Subtract(double b) {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < N; i++) {
    a[i] -= b;
  }
}


//--------------------------------------------------------------------------
void DVector::
Subtract(const DVector &b) {
  if (N != b.GetN()) {
    printf("DVector::Subtract. Error: Vector length does not match. Function call takes no effect.\n");
    return;
  }
  double *mb = b.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < N; i++) {
    a[i] -= mb[i];
  }
}


//--------------------------------------------------------------------------
void DVector::
Subtract(double b, DVector &c) const {
  c.Init(N);
  double *mc = c.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < N; i++) {
    mc[i] = a[i] - b;
  }
}


//--------------------------------------------------------------------------
void DVector::
Subtract(const DVector &b, DVector &c) const {
  if (N != b.GetN()) {
    printf("DVector::Subtract. Error: Vector length does not match. Function call takes no effect.\n");
    return;
  }
  c.Init(N);
  double *mc = c.GetPointer();
  double *mb = b.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < N; i++) {
    mc[i] = a[i] - mb[i];
  }
}


//--------------------------------------------------------------------------
void DVector::
Multiply(double b) {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < N; i++) {
    a[i] *= b;
  }
}


//--------------------------------------------------------------------------
void DVector::
Multiply(const DVector &b) {
  if (N != b.GetN()) {
    printf("DVector::Multiply. Error: Vector length does not match. Function call takes no effect.\n");
    return;
  }
  double *mb = b.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < N; i++) {
    a[i] *= mb[i];
  }
}


//--------------------------------------------------------------------------
void DVector::
Multiply(double b, DVector &c) const {
  c.Init(N);
  double *mc = c.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < N; i++) {
    mc[i] = a[i] * b;
  }
}


//--------------------------------------------------------------------------
void DVector::
Multiply(const DVector &b, DVector &c) const {
  if (N != b.GetN()) {
    printf("DVector::Multiply. Error: Vector length does not match. Function call takes no effect.\n");
    return;
  }
  c.Init(N);
  double *mc = c.GetPointer();
  double *mb = b.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < N; i++) {
    mc[i] = a[i] * mb[i];
  }
}


//--------------------------------------------------------------------------
void DVector::
Divide(double b) {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < N; i++) {
    a[i] /= b;
  }
}


//--------------------------------------------------------------------------
void DVector::
Divide(const DVector &b) {
  if (N != b.GetN()) {
    printf("DVector::Divide. Error: Vector length does not match. Function call takes no effect.\n");
    return;
  }
  double *mb = b.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < N; i++) {
    a[i] /= mb[i];
  }
}


//--------------------------------------------------------------------------
void DVector::
Divide(double b, DVector &c) const {
  c.Init(N);
  double *mc = c.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < N; i++) {
    mc[i] = a[i] / b;
  }
}


//--------------------------------------------------------------------------
void DVector::
Divide(const DVector &b, DVector &c) const {
  if (N != b.GetN()) {
    printf("DVector::Divide. Error: Vector length does not match. Function call takes no effect.\n");
    return;
  }
  c.Init(N);
  double *mc = c.GetPointer();
  double *mb = b.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < N; i++) {
    mc[i] = a[i] / mb[i];
  }
}


//--------------------------------------------------------------------------
double DVector::
InProd(const DVector &b) const {
  if (N != b.GetN()) {
    printf("DVector::InProd. Error: Vector length does not match. Return NAN.\n");
    return NAN;
  }
  double *mb = b.GetPointer();
  double ret = ddot_(&N, a, &ONEi, mb, &ONEi);
  return ret;
}


//--------------------------------------------------------------------------
void DVector::
Abs(void) {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < N; i++) {
    a[i] = fabs(a[i]);
  }
}


//--------------------------------------------------------------------------
void DVector::
Abs(DVector &b) const {
  b.Init(N);
  double *mb = b.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < N; i++) {
    mb[i] = fabs(a[i]);
  }
}


//--------------------------------------------------------------------------
void DVector::
Sqrt(void) {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < N; i++) {
    a[i] = sqrt(a[i]);
  }
}


//--------------------------------------------------------------------------
void DVector::
Sqrt(DVector &b) const {
  b.Init(N);
  double *mb = b.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < N; i++) {
    mb[i] = sqrt(a[i]);
  }
}


//--------------------------------------------------------------------------
void DVector::
Square(void) {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < N; i++) {
    a[i] = a[i] * a[i];
  }
}


//--------------------------------------------------------------------------
void DVector::
Square(DVector &b) const {
  b.Init(N);
  double *mb = b.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < N; i++) {
    mb[i] = a[i] * a[i];
  }
}


//--------------------------------------------------------------------------
void DVector::
Inv(void) {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < N; i++) {
    a[i] = 1.0/a[i];
  }
}


//--------------------------------------------------------------------------
void DVector::
Inv(DVector &b) const {
  b.Init(N);
  double *mb = b.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < N; i++) {
    mb[i] = 1.0/a[i];
  }
}


//--------------------------------------------------------------------------
double DVector::
Norm2(void) const {
  double mNorm2 = dnrm2_(&N, a, &ONEi);
  return mNorm2;
}


//--------------------------------------------------------------------------
double DVector::
Norm1(void) const {
  double mNorm1 = dasum_(&N, a, &ONEi);
  return mNorm1;
}


//--------------------------------------------------------------------------
double DVector::
Min(void) const {
  INTEGER idx;
  double mMin = Min(idx);
  return mMin;
}


//--------------------------------------------------------------------------
double DVector::
Min(INTEGER &idx) const {
  double mMin = DBL_MAX;
  idx = -1;
  for (INTEGER i = 0; i < N; i++) {
    if (mMin > a[i]) {
      mMin = a[i];
      idx = i;
    }
  }
  return mMin;
}


//--------------------------------------------------------------------------
double DVector::
Max(void) const {
  INTEGER idx;
  double mMax = Max(idx);
  return mMax;
}


//--------------------------------------------------------------------------
double DVector::
Max(INTEGER &idx) const {
  double mMax = -DBL_MAX;
  idx = -1;
  for (INTEGER i = 0; i < N; i++) {
    if (mMax < a[i]) {
      mMax = a[i];
      idx = i;
    }
  }
  return mMax;
}


//--------------------------------------------------------------------------
double DVector::
Sum(void) const {
  double mSum = 0.0;
#ifdef USE_OPENMP
#pragma omp parallel for reduction(+:mSum)
#endif
  for (INTEGER i = 0; i < N; i++) {
    mSum += a[i];
  }
  return mSum;
}


//--------------------------------------------------------------------------
double DVector::
Mean(void) const {
  double mMean = Sum()/N;
  return mMean;
}
