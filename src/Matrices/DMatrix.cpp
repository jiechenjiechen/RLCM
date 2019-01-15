#include "DMatrix.hpp"

#define INITVAL_M 0
#define INITVAL_N 0


//--------------------------------------------------------------------------
DMatrix::
DMatrix() {
  A = NULL;
  IPIV = NULL;
  M = INITVAL_M;
  N = INITVAL_N;
  Init();
}


//--------------------------------------------------------------------------
void DMatrix::
Init(void) {
  Init(0);
}


//--------------------------------------------------------------------------
void DMatrix::
Init(INTEGER N_) {
  Init(N_, N_);
}


//--------------------------------------------------------------------------
void DMatrix::
Init(INTEGER M_, INTEGER N_) {
  if (M != M_ || N != N_) {
    ReleaseAllMemory();
  }
  if (!A) {
    New_1D_Array<double, INTEGER>(&A, M_*N_);
  }
  else {
    memset(A, 0, M_*N_*sizeof(double));
  }
  M = M_;
  N = N_;
  // Memory of IPIV is allocated only when needed 
  Delete_1D_Array<INTEGER>(&IPIV);
}


//--------------------------------------------------------------------------
void DMatrix::
ReleaseAllMemory(void) {
  Delete_1D_Array<double>(&A);
  Delete_1D_Array<INTEGER>(&IPIV);
  M = INITVAL_M;
  N = INITVAL_N;
}


//--------------------------------------------------------------------------
DMatrix::
DMatrix(INTEGER N_) {
  M = N_;
  N = N_;
  A = NULL;
  IPIV = NULL;
  Init(N_);
}


//--------------------------------------------------------------------------
DMatrix::
DMatrix(INTEGER M_, INTEGER N_) {
  M = M_;
  N = N_;
  A = NULL;
  IPIV = NULL;
  Init(M_, N_);
}


//--------------------------------------------------------------------------
DMatrix::
DMatrix(const DMatrix &G) {
  M = INITVAL_M;
  N = INITVAL_N;
  A = NULL;
  IPIV = NULL;
  Init();
  DeepCopy(G);
}


//--------------------------------------------------------------------------
DMatrix& DMatrix::
operator= (const DMatrix &G) {
  if (this != &G) {
    DeepCopy(G);
  }
  return *this;
}


//--------------------------------------------------------------------------
void DMatrix::
DeepCopy(const DMatrix &G) {
  if (M != G.M || N != G.N) {
    ReleaseAllMemory();
  }
  if (G.A) {
    if (!A) {
      New_1D_Array<double, INTEGER>(&A, G.M*G.N);
    }
    memcpy(A, G.A, G.M*G.N*sizeof(double));
  }
  if (G.IPIV) {
    if (!IPIV) {
      New_1D_Array<INTEGER, INTEGER>(&IPIV, G.M<G.N?G.M:G.N);
    }
    memcpy(IPIV, G.IPIV, (G.M<G.N?G.M:G.N)*sizeof(INTEGER));
  }
  else if (IPIV) {
    Delete_1D_Array<INTEGER>(&IPIV);
  }
  M = G.M;
  N = G.N;
}


//--------------------------------------------------------------------------
DMatrix::
~DMatrix() {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
INTEGER DMatrix::
GetM(void) const {
  return M;
}


//--------------------------------------------------------------------------
INTEGER DMatrix::
GetN(void) const {
  return N;
}


//--------------------------------------------------------------------------
double DMatrix::
GetEntry(INTEGER i, INTEGER j) const {
  if (i < 0 || i >= M || j < 0 || j >= N) {
    printf("DMatrix::GetEntry. Error: Invalid index. Return NAN.\n");
    return NAN;
  }
  return A[i+M*j];
}


//--------------------------------------------------------------------------
void DMatrix::
GetColumn(INTEGER i, DVector &b) const {
  if (i < 0 || i >= N) {
    printf("DMatrix::GetColumn. Error: Invalid column index. Function call takes no effect.\n");
    return;
  }
  b.Init(M);
  double *mb = b.GetPointer();
  dcopy_(&M, A+M*i, &ONEi, mb, &ONEi);
}


//--------------------------------------------------------------------------
void DMatrix::
GetColumns(INTEGER *idx, INTEGER n, DMatrix &B) const {
  B.Init(M, n);
  double *mB = B.GetPointer();
  for (INTEGER j = 0; j < n; j++) {
    dcopy_(&M, A+M*idx[j], &ONEi, mB+M*j, &ONEi);
  }
}


//--------------------------------------------------------------------------
void DMatrix::
GetRow(INTEGER i, DVector &b) const {
  if (i < 0 || i >= M) {
    printf("DMatrix::GetRow. Error: Invalid row index. Function call takes no effect.\n");
    return;
  }
  b.Init(N);
  double *mb = b.GetPointer();
  dcopy_(&N, A+i, &M, mb, &ONEi);
}


//--------------------------------------------------------------------------
void DMatrix::
GetBlock(INTEGER RowStart, INTEGER nRow, INTEGER ColStart, INTEGER nCol,
         DMatrix &B) const {
  if (RowStart < 0 || RowStart+nRow > M || ColStart < 0 || ColStart+nCol > N) {
    printf("DMatrix::GetBlock. Error: Invalid range. Function call takes no effect.\n");
    return;
  }
  B.Init(nRow, nCol);
  double *mB = B.GetPointer();
  for (INTEGER j = 0; j < nCol; j++) {
    dcopy_(&nRow, A+M*(ColStart+j)+RowStart, &ONEi, mB+nRow*j, &ONEi);
  }
}


//--------------------------------------------------------------------------
double* DMatrix::
GetPointer(void) const {
  return A;
}


//--------------------------------------------------------------------------
void DMatrix::
SetEntry(INTEGER i, INTEGER j, double b) {
  if (i < 0 || i >= M || j < 0 || j >= N) {
    printf("DMatrix::SetEntry. Error: Invalid index. Function call takes no effect.\n");
    return;
  }
  A[i+M*j] = b;
}


//--------------------------------------------------------------------------
void DMatrix::
AddToEntry(INTEGER i, INTEGER j, double b) {
  if (i < 0 || i >= M || j < 0 || j >= N) {
    printf("DMatrix::AddToEntry. Error: Invalid index. Function call takes no effect.\n");
    return;
  }
  A[i+M*j] += b;
}


//--------------------------------------------------------------------------
void DMatrix::
SetColumn(INTEGER i, const DVector &b) {
  if (b.GetN() != M) {
    printf("DMatrix::SetColumn. Error: Column length not compatible. Function call takes no effect.\n");
    return;
  }
  double *mb = b.GetPointer();
  dcopy_(&M, mb, &ONEi, A+M*i, &ONEi);
}


//--------------------------------------------------------------------------
void DMatrix::
SetRow(INTEGER i, const DVector &b) {
  if (b.GetN() != N) {
    printf("DMatrix::SetRow. Error: Row length not compatible. Function call takes no effect.\n");
    return;
  }
  double *mb = b.GetPointer();
  dcopy_(&N, mb, &ONEi, A+i, &M);
}


//--------------------------------------------------------------------------
void DMatrix::
SetBlock(INTEGER RowStart, INTEGER nRow, INTEGER ColStart, INTEGER nCol,
         const DMatrix &B) {
  if (RowStart < 0 || RowStart+nRow > M || ColStart < 0 || ColStart+nCol > N) {
    printf("DMatrix::SetBlock. Error: Invalid range. Function call takes no effect.\n");
    return;
  }
  if (nRow != B.GetM() || nCol != B.GetN()) {
    printf("DMatrix::SetBlock. Error: Matrix dimensions not compatible. Function call takes no effect.\n");
    return;
  }
  double *mB = B.GetPointer();
  for (INTEGER j = 0; j < nCol; j++) {
    dcopy_(&nRow, mB+nRow*j, &ONEi, A+M*(ColStart+j)+RowStart, &ONEi);
  }
}


//--------------------------------------------------------------------------
void DMatrix::
SetBlock(INTEGER RowStart, INTEGER nRow, INTEGER ColStart, INTEGER nCol,
         const double *B) {
  if (RowStart < 0 || RowStart+nRow > M || ColStart < 0 || ColStart+nCol > N) {
    printf("DMatrix::SetBlock. Error: Invalid range. Function call takes no effect.\n");
    return;
  }
  if (B == NULL) {
    printf("DMatrix::SetBlock. Error: double *B is an empty pointer. Function call takes no effect.\n");
    return;
  }
  for (INTEGER j = 0; j < nCol; j++) {
    dcopy_(&nRow, B+nRow*j, &ONEi, A+M*(ColStart+j)+RowStart, &ONEi);
  }
}


//--------------------------------------------------------------------------
void DMatrix::
SetIdentity(void) {
  memset(A, 0, M*N*sizeof(double));
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER j = 0; j < (M>N?N:M); j++) {
    A[M*j+j] = 1.0;
  }
}


//--------------------------------------------------------------------------
void DMatrix::
SetMultIdentity(double lambda) {
  memset(A, 0, M*N*sizeof(double));
  INTEGER minMN = M>N?N:M;
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER j = 0; j < minMN; j++) {
    A[M*j+j] = lambda;
  }
}


//--------------------------------------------------------------------------
void DMatrix::
SetConstVal(double c) {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < M*N; i++) {
    A[i] = c;
  }
}


//--------------------------------------------------------------------------
void DMatrix::
SetUniformRandom01(void) {
  UniformRandom01(A, M*N);
}


//--------------------------------------------------------------------------
void DMatrix::
SetStandardNormal(void) {
  StandardNormal(A, M*N);
}


//--------------------------------------------------------------------------
void DMatrix::
SetStudentT1(void) {
  StudentT1(A, M*N);
}


//--------------------------------------------------------------------------
void DMatrix::
SetMultivariateStudentT1(void) {
  for (INTEGER i = 0; i < N; i++) {
    DVector y;
    GetColumn(i, y);
    y.SetMultivariateStudentT1();
    SetColumn(i, y);
  }
}


//--------------------------------------------------------------------------
void DMatrix::
SetRandomSech(void) {
  RandomSech(A, M*N);
}


//--------------------------------------------------------------------------
void DMatrix::
MakeDiag(const DVector &b) {
  Init(b.GetN());
  double *mb = b.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER j = 0; j < N; j++) {
    A[N*j+j] = mb[j];
  }
}


//--------------------------------------------------------------------------
void DMatrix::
PrintMatrixMatlabForm(const char *name) const {
  double *ptr = A;
  printf("%s = [", name);
  for (INTEGER i = 0; i < M; i++) {
    ptr = A + i;
    for (INTEGER j = 0; j < N; j++) {
      printf("%g ", *ptr);
      ptr += M;
    }
    if (i != M-1) {
      printf(";\n");
    }
    else {
      printf("];\n");
    }
  }
  printf("\n");
}


//--------------------------------------------------------------------------
void DMatrix::
Transpose(void) {
  double *B = NULL;
  New_1D_Array<double, INTEGER>(&B, M*N);
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER j = 0; j < N; j++) {
    for (INTEGER i = 0; i < M; i++) {
      B[i*N+j] = A[i+j*M];
    }
  }
  INTEGER MN = M*N;
  dcopy_(&MN, B, &ONEi, A, &ONEi);
  Swap<INTEGER>(M, N);
  Delete_1D_Array<double>(&B);
}


//--------------------------------------------------------------------------
void DMatrix::
Transpose(DMatrix &B) const {
  B.Init(N,M);
  double *mB = B.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER j = 0; j < N; j++) {
    for (INTEGER i = 0; i < M; i++) {
      mB[i*N+j] = A[i+j*M];
    }
  }
}


//--------------------------------------------------------------------------
void DMatrix::
Symmetrize(void) {
  if (M != N) {
    printf("DMatrix::Symmetrize. Error: Matrix is not square. Function call takes no effect.\n");
    return;
  }
  for (INTEGER j = 0; j < N; j++) {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (INTEGER i = 0; i <= j; i++) {
      INTEGER id1 = i*N+j;
      INTEGER id2 = i+j*N;
      A[id1] = (A[id1]+A[id2])/2.0;
      A[id2] = A[id1];
    }
  }
}


//--------------------------------------------------------------------------
void DMatrix::
Symmetrize(DMatrix &B) const {
  if (M != N) {
    printf("DMatrix::Symmetrize. Error: Matrix is not square. Function call takes no effect.\n");
    return;
  }
  B.Init(M,N);
  double *mB = B.GetPointer();
  for (INTEGER j = 0; j < N; j++) {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (INTEGER i = 0; i <= j; i++) {
      INTEGER id1 = i*N+j;
      INTEGER id2 = i+j*N;
      mB[id1] = (A[id1]+A[id2])/2.0;
      mB[id2] = mB[id1];
    }
  }
}


//--------------------------------------------------------------------------
void DMatrix::
Negate(void) {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < M*N; i++) {
    A[i] = -A[i];
  }
}


//--------------------------------------------------------------------------
void DMatrix::
Negate(DMatrix &B) const {
  B.Init(M,N);
  double *mB = B.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < M*N; i++) {
    mB[i] = -A[i];
  }
}


//--------------------------------------------------------------------------
void DMatrix::
Add(const DMatrix &B) {
  if (M != B.GetM() || N != B.GetN()) {
    printf("DMatrix::Add. Error: Matrix dimensions do not match. Function call takes no effect.\n");
    return;
  }
  double *mB = B.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < M*N; i++) {
    A[i] += mB[i];
  }
}


//--------------------------------------------------------------------------
void DMatrix::
Add(const DMatrix &B, DMatrix &C) const {
  if (M != B.GetM() || N != B.GetN()) {
    printf("DMatrix::Add. Error: Matrix dimensions do not match. Function call takes no effect.\n");
    return;
  }
  C.Init(M,N);
  double *mC = C.GetPointer();
  double *mB = B.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < M*N; i++) {
    mC[i] = A[i] + mB[i];
  }
}


//--------------------------------------------------------------------------
void DMatrix::
Subtract(const DMatrix &B) {
  if (M != B.GetM() || N != B.GetN()) {
    printf("DMatrix::Subtract. Error: Matrix dimensions do not match. Function call takes no effect.\n");
    return;
  }
  double *mB = B.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < M*N; i++) {
    A[i] -= mB[i];
  }
}


//--------------------------------------------------------------------------
void DMatrix::
Subtract(const DMatrix &B, DMatrix &C) const {
  if (M != B.GetM() || N != B.GetN()) {
    printf("DMatrix::Subtract. Error: Matrix dimensions do not match. Function call takes no effect.\n");
    return;
  }
  C.Init(M,N);
  double *mC = C.GetPointer();
  double *mB = B.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < M*N; i++) {
    mC[i] = A[i] - mB[i];
  }
}


//--------------------------------------------------------------------------
void DMatrix::
Add(double b) {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < M*N; i++) {
    A[i] += b;
  }
}


//--------------------------------------------------------------------------
void DMatrix::
Add(double b, DMatrix &C) const {
  C.Init(M,N);
  double *mC = C.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < M*N; i++) {
    mC[i] = A[i] + b;
  }
}


//--------------------------------------------------------------------------
void DMatrix::
Subtract(double b) {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < M*N; i++) {
    A[i] -= b;
  }
}


//--------------------------------------------------------------------------
void DMatrix::
Subtract(double b, DMatrix &C) const {
  C.Init(M,N);
  double *mC = C.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < M*N; i++) {
    mC[i] = A[i] - b;
  }
}


//--------------------------------------------------------------------------
void DMatrix::
Multiply(const DMatrix &B) {
  if (M != B.GetM() || N != B.GetN()) {
    printf("DMatrix::Multiply. Error: Matrix dimensions do not match. Function call takes no effect.\n");
    return;
  }
  double *mB = B.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < M*N; i++) {
    A[i] *= mB[i];
  }
}


//--------------------------------------------------------------------------
void DMatrix::
Multiply(const DMatrix &B, DMatrix &C) const {
  if (M != B.GetM() || N != B.GetN()) {
    printf("DMatrix::Multiply. Error: Matrix dimensions do not match. Function call takes no effect.\n");
    return;
  }
  C.Init(M,N);
  double *mC = C.GetPointer();
  double *mB = B.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < M*N; i++) {
    mC[i] = A[i] * mB[i];
  }
}


//--------------------------------------------------------------------------
void DMatrix::
Divide(const DMatrix &B) {
  if (M != B.GetM() || N != B.GetN()) {
    printf("DMatrix::Divide. Error: Matrix dimensions do not match. Function call takes no effect.\n");
    return;
  }
  double *mB = B.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < M*N; i++) {
    A[i] /= mB[i];
  }
}


//--------------------------------------------------------------------------
void DMatrix::
Divide(const DMatrix &B, DMatrix &C) const {
  if (M != B.GetM() || N != B.GetN()) {
    printf("DMatrix::Divide. Error: Matrix dimensions do not match. Function call takes no effect.\n");
    return;
  }
  C.Init(M,N);
  double *mC = C.GetPointer();
  double *mB = B.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < M*N; i++) {
    mC[i] = A[i] / mB[i];
  }
}


//--------------------------------------------------------------------------
void DMatrix::
Multiply(double b) {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < M*N; i++) {
    A[i] *= b;
  }
}


//--------------------------------------------------------------------------
void DMatrix::
Multiply(double b, DMatrix &C) const {
  C.Init(M,N);
  double *mC = C.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < M*N; i++) {
    mC[i] = A[i] * b;
  }
}


//--------------------------------------------------------------------------
void DMatrix::
Divide(double b) {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < M*N; i++) {
    A[i] /= b;
  }
}


//--------------------------------------------------------------------------
void DMatrix::
Divide(double b, DMatrix &C) const {
  C.Init(M,N);
  double *mC = C.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < M*N; i++) {
    mC[i] = A[i] / b;
  }
}


//--------------------------------------------------------------------------
void DMatrix::
AddDiagonal(double s) {
  INTEGER minMN = M<N?M:N;
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER j = 0; j < minMN; j++) {
    A[M*j+j] += s;
  }
}


//--------------------------------------------------------------------------
void DMatrix::
AddDiagonal(double s, DMatrix &B) const {
  B = *this;
  double *mB = B.GetPointer();
  INTEGER minMN = M<N?M:N;
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER j = 0; j < minMN; j++) {
    INTEGER idx = M*j+j;
    mB[idx] = A[idx] + s;
  }
}


//--------------------------------------------------------------------------
void DMatrix::
SubtractDiagonal(double s) {
  INTEGER minMN = M<N?M:N;
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER j = 0; j < minMN; j++) {
    A[M*j+j] -= s;
  }
}


//--------------------------------------------------------------------------
void DMatrix::
SubtractDiagonal(double s, DMatrix &B) const {
  B = *this;
  double *mB = B.GetPointer();
  INTEGER minMN = M<N?M:N;
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER j = 0; j < minMN; j++) {
    INTEGER idx = M*j+j;
    mB[idx] = A[idx] - s;
  }
}


//--------------------------------------------------------------------------
void DMatrix::
Cos(void) {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER j = 0; j < M*N; j++) {
    A[j] = cos(A[j]);
  }
}


//--------------------------------------------------------------------------
void DMatrix::
Cos(DMatrix &B) const {
  B.Init(M,N);
  double *mB = B.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER j = 0; j < M*N; j++) {
    mB[j] = cos(A[j]);
  }
}


//--------------------------------------------------------------------------
void DMatrix::
Sqrt(void) {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER j = 0; j < M*N; j++) {
    A[j] = sqrt(A[j]);
  }
}


//--------------------------------------------------------------------------
void DMatrix::
Sqrt(DMatrix &B) const {
  B.Init(M,N);
  double *mB = B.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER j = 0; j < M*N; j++) {
    mB[j] = sqrt(A[j]);
  }
}


//--------------------------------------------------------------------------
void DMatrix::
Log(void) {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER j = 0; j < M*N; j++) {
    A[j] = log(A[j]);
  }
}


//--------------------------------------------------------------------------
void DMatrix::
Log(DMatrix &B) const {
  B.Init(M,N);
  double *mB = B.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER j = 0; j < M*N; j++) {
    mB[j] = log(A[j]);
  }
}


//--------------------------------------------------------------------------
void DMatrix::
Sum(DVector &b, INTEGER dim) const {
  if (dim == 1) {
    b.Init(N);
    double *mb = b.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (INTEGER j = 0; j < N; j++) {
      for (INTEGER i = 0; i < M; i++) {
        mb[j] += A[i+j*M];
      }
    }
  }
  else if (dim == 2) {
    b.Init(M);
    double *mb = b.GetPointer();
    for (INTEGER j = 0; j < N; j++) {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
      for (INTEGER i = 0; i < M; i++) {
        mb[i] += A[i+j*M];
      }
    }
  }
  else {
    printf("DMatrix::Sum. Error: Argument dim can only be 1 or 2. Function call takes no effect.\n");
    return;
  }
}


//--------------------------------------------------------------------------
void DMatrix::
Prod(DVector &b, INTEGER dim) const {
  if (dim == 1) {
    b.Init(N);
    b.SetConstVal(1.0);
    double *mb = b.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (INTEGER j = 0; j < N; j++) {
      for (INTEGER i = 0; i < M; i++) {
        mb[j] *= A[i+j*M];
      }
    }
  }
  else if (dim == 2) {
    b.Init(M);
    b.SetConstVal(1.0);
    double *mb = b.GetPointer();
    for (INTEGER j = 0; j < N; j++) {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
      for (INTEGER i = 0; i < M; i++) {
        mb[i] *= A[i+j*M];
      }
    }
  }
  else {
    printf("DMatrix::Prod. Error: Argument dim can only be 1 or 2. Function call takes no effect.\n");
    return;
  }
}


//--------------------------------------------------------------------------
void DMatrix::
OuterProduct(const DVector &b, const DVector &c) {
  Init(b.GetN(), c.GetN());
  double *mb = b.GetPointer();
  double *mc = c.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER j = 0; j < N; j++) {
    for (INTEGER i = 0; i < M; i++) {
      A[i+j*M] = mb[i] * mc[j];
    }
  }
}


//--------------------------------------------------------------------------
void DMatrix::
MatVec(const DVector &b, DVector &y, MatrixMode ModeA, MatState mState) const {

  if (mState != UNFACT) {
    printf("DMatrix::MatVec. Error: mState must be UNFACT. Function call takes no effect.\n");
    return;
  }

  if (ModeA == NORMAL) {
    y.Init(M);
  }
  else {
    y.Init(N);
  }
  DGEMV(b, y, 1.0, 0.0, ModeA);

}


//--------------------------------------------------------------------------
void DMatrix::
DGEMV(const DVector &b, DVector &y,
      double alpha, double beta, MatrixMode ModeA) const {

  if ( (ModeA == NORMAL && (N != b.GetN() || M != y.GetN())) ||
       (ModeA != NORMAL && (M != b.GetN() || N != y.GetN())) ) {
    printf("DMatrix::DGEMV. Error: Matrix/vector sizes are not compatible. Function call takes no effect.\n");
    return;
  }

  char TRANS;
  switch (ModeA) {
  case NORMAL: TRANS = 'N'; break;
  case TRANSPOSE: TRANS = 'T'; break;
  case CONJ_TRANS: TRANS = 'C'; break;
  }

  double *mb = b.GetPointer();
  double *my = y.GetPointer();

  dgemv_(&TRANS, &M, &N, &alpha, A, &M, mb, &ONEi, &beta, my, &ONEi);

}


//--------------------------------------------------------------------------
void DMatrix::
MatMat(const DMatrix &B, DMatrix &C, MatrixMode ModeA, MatrixMode ModeB) const {
  if (ModeA == NORMAL && ModeB == NORMAL) {
    C.Init(M, B.GetN());
  }
  else if (ModeA == NORMAL && ModeB != NORMAL) {
    C.Init(M, B.GetM());
  }
  else if (ModeA != NORMAL && ModeB == NORMAL) {
    C.Init(N, B.GetN());
  }
  else /*(ModeA != NORMAL && ModeB != NORMAL)*/ {
    C.Init(N, B.GetM());
  }
  DGEMM(B, C, 1.0, 0.0, ModeA, ModeB);
}


//--------------------------------------------------------------------------
void DMatrix::
DGEMM(const DMatrix &B, DMatrix &C, double alpha, double beta,
      MatrixMode ModeA, MatrixMode ModeB) const {

  if ( (ModeA == NORMAL && ModeB == NORMAL &&
        (N != B.GetM() || C.GetM() != M || C.GetN() != B.GetN())) ||
       (ModeA == NORMAL && ModeB != NORMAL &&
        (N != B.GetN() || C.GetM() != M || C.GetN() != B.GetM())) ||
       (ModeA != NORMAL && ModeB == NORMAL &&
        (M != B.GetM() || C.GetM() != N || C.GetN() != B.GetN())) ||
       (ModeA != NORMAL && ModeB != NORMAL &&
        (M != B.GetN() || C.GetM() != N || C.GetN() != B.GetM())) ) {
    printf("DMatrix::DGEMM. Error: Matrix sizes are not compatible. Function call takes no effect.\n");
    return;
  }

  INTEGER mM, mN, mK;
  char TRANSA;
  switch (ModeA) {
  case NORMAL: TRANSA = 'N'; mM = M; mK = N; break;
  case TRANSPOSE: TRANSA = 'T'; mM = N; mK = M; break;
  case CONJ_TRANS: TRANSA = 'C'; mM = N; mK = M; break;
  }
  char TRANSB;
  switch (ModeB) {
  case NORMAL: TRANSB = 'N'; mN = B.GetN(); break;
  case TRANSPOSE: TRANSB = 'T'; mN = B.GetM(); break;
  case CONJ_TRANS: TRANSB = 'C'; mN = B.GetM(); break;
  }

  INTEGER LDA = M, LDB = B.GetM();
  double *mB = B.GetPointer();
  double *mC = C.GetPointer();

  dgemm_(&TRANSA, &TRANSB, &mM, &mN, &mK, &alpha, A, &LDA, mB, &LDB, &beta, mC, &mM);

}


//--------------------------------------------------------------------------
void DMatrix::
Mldivide(const DVector &b, DVector &x,
         MatrixMode ModeA, MatrixType TypeA) {
  DMatrix B(b.GetN(),1), X(b.GetN(),1);
  B.SetColumn(0,b);
  Mldivide(B, X, ModeA, TypeA);
  X.GetColumn(0,x);
}


//--------------------------------------------------------------------------
void DMatrix::
Mldivide(const DMatrix &B, DMatrix &X,
         MatrixMode ModeA, MatrixType TypeA) {
  if ( M != N ) {
    printf("DMatrix::Mldivide. Error: Matrix is not square. Function call takes no effect.\n");
    return;
  }
  Mldivide_LinearSystemSolve(B, X, ModeA, TypeA);
}


//--------------------------------------------------------------------------
void DMatrix::
Mldivide(const DVector &b, DVector &x, double *Res) {
  DMatrix B(b.GetN(),1), X(b.GetN(),1);
  B.SetColumn(0,b);
  Mldivide(B, X, Res);
  X.GetColumn(0,x);
}


//--------------------------------------------------------------------------
void DMatrix::
Mldivide(const DMatrix &B, DMatrix &X, double *Res) {
  Mldivide_LeastSquaresSolve(B, X, Res);
}


//--------------------------------------------------------------------------
void DMatrix::
Mldivide_LinearSystemSolve(const DMatrix &B, DMatrix &X,
                           MatrixMode ModeA, MatrixType TypeA) {

  if ( N != B.GetM() ) {
    printf("DMatrix::Mldivide_LinearSystemSolve. Error: Matrix/vector sizes are not compatible. Function call takes no effect.\n");
    return;
  }

  INTEGER NRHS = B.GetN();
  double RCOND = 0.0;
  INTEGER INFO = 0;
  char EQUED_N = 'N';

  double *mB = B.GetPointer();
  X.Init(N, NRHS);
  double *mX = X.GetPointer();
  DMatrix AF = *this;
  double *mAF = AF.GetPointer();

  double *FERR = NULL;
  New_1D_Array<double, INTEGER>(&FERR, NRHS);
  double *BERR = NULL;
  New_1D_Array<double, INTEGER>(&BERR, NRHS);
  double *WORK = NULL;
  INTEGER *IWORK = NULL;
  New_1D_Array<INTEGER, INTEGER>(&IWORK, N);

  switch (TypeA) {
  case GENERAL: {

    char TRANS;
    switch (ModeA) {
    case NORMAL: TRANS = 'N'; break;
    case TRANSPOSE: TRANS = 'T'; break;
    case CONJ_TRANS: TRANS = 'C'; break;
    }

    Delete_1D_Array<INTEGER>(&IPIV);
    New_1D_Array<INTEGER, INTEGER>(&IPIV, N);

    New_1D_Array<double, INTEGER>(&WORK, 4*N);

    dgesvx_(&FACT_N, &TRANS, &N, &NRHS, A, &N, mAF, &N, IPIV, &EQUED_N, NULL, NULL, mB, &N, mX, &N, &RCOND, FERR, BERR, WORK, IWORK, &INFO);

    if (INFO < 0) {
      printf("DMatrix::Mldivide_LinearSystemSolve. DGESVX error: INFO = %ld. Computed result cannot be trusted.\n", (long)INFO);
    }
    else if (INFO > 0 && INFO <= N) {
      printf("DMatrix::Mldivide_LinearSystemSolve. DGESVX error: INFO = %ld. Matrix is singular. Computed result cannot be trusted.\n", (long)INFO);
    }
    else if (INFO == N+1) {
      printf("DMatrix::Mldivide_LinearSystemSolve. DGESVX warning: INFO = %ld. Matrix is singular to working precision.\n", (long)INFO);
    }

    break;
  }
  case SPD: {

    New_1D_Array<double, INTEGER>(&WORK, 3*N);

    dposvx_(&FACT_N, &UPLO_U, &N, &NRHS, A, &N, mAF, &N, &EQUED_N, NULL, mB, &N, mX, &N, &RCOND, FERR, BERR, WORK, IWORK, &INFO);

    if (INFO < 0) {
      printf("DMatrix::Mldivide_LinearSystemSolve. DPOSVX error: INFO = %ld. Computed result cannot be trusted.\n", (long)INFO);
    }
    else if (INFO > 0 && INFO <= N) {
      printf("DMatrix::Mldivide_LinearSystemSolve. DPOSVX error: INFO = %ld. Matrix is not positive definite. Computed result cannot be trusted.\n", (long)INFO);
    }
    else if (INFO == N+1) {
      printf("DMatrix::Mldivide_LinearSystemSolve. DPOSVX warning: INFO = %ld. Matrix is singular to working precision.\n", (long)INFO);
    }

    break;
  }

  }

  Delete_1D_Array<double>(&FERR);
  Delete_1D_Array<double>(&BERR);
  Delete_1D_Array<double>(&WORK);
  Delete_1D_Array<INTEGER>(&IWORK);

}


//--------------------------------------------------------------------------
void DMatrix::
Mldivide_LeastSquaresSolve(const DMatrix &B, DMatrix &X, double *Res) {

  if ( M != B.GetM() ) {
    printf("DMatrix::Mldivide_LeastSquaresSolve. Error: Matrix/vector sizes are not compatible. Function call takes no effect.\n");
    return;
  }

  DMatrix A_copy = *this;
  double *mA = A_copy.GetPointer();

  INTEGER NRHS = B.GetN();
  double RCOND = 0.0;
  INTEGER RANK = 0;
  INTEGER INFO = 0;
  double WORK0 = 0;
  INTEGER LWORK = -1;

  double *mB = NULL, *B_ptr = B.GetPointer();
  INTEGER maxMN = M>N?M:N;
  New_1D_Array<double, INTEGER>(&mB, maxMN*NRHS);
  for (INTEGER i = 0; i < NRHS; i++) {
    dcopy_(&M, B_ptr+M*i, &ONEi, mB+maxMN*i, &ONEi);
  }

  double *S = NULL;
  INTEGER minMN = M<N?M:N;
  New_1D_Array<double, INTEGER>(&S, minMN);

  dgelss_(&M, &N, &NRHS, mA, &M, mB, &maxMN, S, &RCOND, &RANK, &WORK0, &LWORK, &INFO);

  LWORK = (INTEGER)WORK0;
  double *WORK = NULL;
  New_1D_Array<double, INTEGER>(&WORK, LWORK);

  dgelss_(&M, &N, &NRHS, mA, &M, mB, &maxMN, S, &RCOND, &RANK, WORK, &LWORK, &INFO);

  if (INFO != 0) {
    printf("DMatrix::Mldivide_LeastSquaresSolve. DGELSS error: INFO = %ld. Computed result cannot be trusted.\n", (long)INFO);
  }

  // handle X and residual
  X.Init(N, NRHS);
  DVector x(N);
  for (INTEGER i = 0; i < NRHS; i++) {
    x.SetBlock(0, N, mB+i*maxMN);
    X.SetColumn(i, x);
    if (Res != NULL) {
      Res[i] = 0.0;
      if (M > N) {
        for (INTEGER j = N; j < maxMN; j++) {
          Res[i] += Square(mB[i*maxMN+j]);
        }
        Res[i] = sqrt(Res[i]);
      }
    }
  }

  Delete_1D_Array<double>(&mB);
  Delete_1D_Array<double>(&S);
  Delete_1D_Array<double>(&WORK);

}


//--------------------------------------------------------------------------
void DMatrix::
DGETRF(void) {

  Delete_1D_Array<INTEGER>(&IPIV);
  New_1D_Array<INTEGER, INTEGER>(&IPIV, M<N?M:N);

  INTEGER INFO = 0;

  dgetrf_(&M, &N, A, &M, IPIV, &INFO);

  if (INFO < 0) {
    printf("DMatrix::DGETRF. DGETRF error: INFO = %ld. LU factor is not computed.\n", (long)INFO);
    Delete_1D_Array<INTEGER>(&IPIV);
  }
  else if (INFO > 0) {
    printf("DMatrix::DGETRF. DGETRF error: INFO = %ld. Matrix is singular.\n", (long)INFO);
    Delete_1D_Array<INTEGER>(&IPIV);
  }

}


//--------------------------------------------------------------------------
void DMatrix::
DGETRS(const DVector &b, DVector &x, MatrixMode ModeA) const {
  DMatrix B(b.GetN(),1), X(b.GetN(),1);
  B.SetColumn(0,b);
  DGETRS(B, X, ModeA);
  X.GetColumn(0,x);
}


//--------------------------------------------------------------------------
void DMatrix::
DGETRS(const DMatrix &B, DMatrix &X, MatrixMode ModeA) const {

  if ( M != N ) {
    printf("DMatrix::DGETRS. Error: Matrix is not square. Function call takes no effect.\n");
    return;
  }
  if ( N != B.GetM() ) {
    printf("DMatrix::DGETRS. Error: Matrix/vector sizes are not compatible. Function call takes no effect.\n");
    return;
  }
  if (IPIV == NULL) {
    printf("DMatrix::DGETRS. Error: DGETRF must be called before calling this routine. Function call takes no effect.\n");
    return;
  }

  char TRANS;
  switch (ModeA) {
  case NORMAL: TRANS = 'N'; break;
  case TRANSPOSE: TRANS = 'T'; break;
  case CONJ_TRANS: TRANS = 'C'; break;
  }
  X = B;
  double *mX = X.GetPointer();
  INTEGER NRHS = B.GetN();
  INTEGER INFO = 0;

  dgetrs_(&TRANS, &N, &NRHS, A, &N, IPIV, mX, &N, &INFO);

  if (INFO < 0) {
    printf("DMatrix::DGETRS. DGETRS error: INFO = %ld. Computed result cannot be trusted.\n", (long)INFO);
  }

}


//--------------------------------------------------------------------------
double DMatrix::
DGECON(void) const {

  if (IPIV == NULL) {
    printf("DMatrix::DGECON. Error: DGETRF must be called before calling this routine. Return NAN.\n");
    return NAN;
  }

  double *WORK = NULL;
  New_1D_Array<double, INTEGER>(&WORK, 4*N);
  INTEGER *IWORK = NULL;
  New_1D_Array<INTEGER, INTEGER>(&IWORK, N);

  double RCOND = 0.0;
  INTEGER INFO = 0;

  double ANORM = dlange_(&NORM_1, &N, &N, A, &N, WORK);
  dgecon_(&NORM_1, &N, A, &N, &ANORM, &RCOND, WORK, IWORK, &INFO);

  if (INFO < 0) {
    printf("DMatrix::DGECON. DGECON error: INFO = %ld. Computed result cannot be trusted.\n", (long)INFO);
  }

  Delete_1D_Array<double>(&WORK);
  Delete_1D_Array<INTEGER>(&IWORK);

  return RCOND;

}


//--------------------------------------------------------------------------
void DMatrix::
DPOTRF(TriUPLO mTri) {

  char UPLO;
  switch (mTri) {
  case UPPER: UPLO = UPLO_U; break;
  case LOWER: UPLO = UPLO_L; break;
  }

  INTEGER INFO = 0;

  dpotrf_(&UPLO, &N, A, &N, &INFO);

  if (INFO < 0) {
    printf("DMatrix::DPOTRF. DPOTRF error: INFO = %ld. Cholesky factor is not computed.\n", (long)INFO);
  }
  else if (INFO > 0) {
    printf("DMatrix::DPOTRF. DPOTRF error: INFO = %ld. Matrix is not positive definite\n", (long)INFO);
  }

}


//--------------------------------------------------------------------------
void DMatrix::
DPOTRS(const DVector &b, DVector &x, TriUPLO mTri) const {
  DMatrix B(b.GetN(),1), X(b.GetN(),1);
  B.SetColumn(0,b);
  DPOTRS(B, X, mTri);
  X.GetColumn(0,x);
}


//--------------------------------------------------------------------------
void DMatrix::
DPOTRS(const DMatrix &B, DMatrix &X, TriUPLO mTri) const {

  if ( M != N ) {
    printf("DMatrix::DPOTRS. Error: Matrix is not square. Function call takes no effect.\n");
    return;
  }
  if ( N != B.GetM() ) {
    printf("DMatrix::DPOTRS. Error: Matrix/vector sizes are not compatible. Function call takes no effect.\n");
    return;
  }

  char UPLO;
  switch (mTri) {
  case UPPER: UPLO = UPLO_U; break;
  case LOWER: UPLO = UPLO_L; break;
  }

  X = B;
  double *mX = X.GetPointer();
  INTEGER NRHS = B.GetN();
  INTEGER INFO = 0;

  dpotrs_(&UPLO, &N, &NRHS, A, &N, mX, &N, &INFO);

  if (INFO < 0) {
    printf("DMatrix::DPOTRS. DPOTRS error: INFO = %ld. Computed result cannot be trusted.\n", (long)INFO);
  }

}


//--------------------------------------------------------------------------
void DMatrix::
DTRSV(const DVector &b, DVector &x, MatrixMode ModeA, TriUPLO mTri) const {

  if ( M != N ) {
    printf("DMatrix::DTRSV. Error: Matrix is not square. Function call takes no effect.\n");
    return;
  }
  if ( N != b.GetN() ) {
    printf("DMatrix::DTRSV. Error: Matrix/vector sizes are not compatible. Function call takes no effect.\n");
    return;
  }

  char UPLO;
  switch (mTri) {
  case UPPER: UPLO = UPLO_U; break;
  case LOWER: UPLO = UPLO_L; break;
  }

  char TRANS;
  switch (ModeA) {
  case NORMAL: TRANS = 'N'; break;
  case TRANSPOSE: TRANS = 'T'; break;
  case CONJ_TRANS: TRANS = 'C'; break;
  }

  x = b;
  double *mx = x.GetPointer();

  dtrsv_(&UPLO, &TRANS, &DIAG_N, &N, A, &N, mx, &ONEi);

}


//--------------------------------------------------------------------------
void DMatrix::
DTRSM(const DMatrix &B, DMatrix &X, MatrixMode ModeA, TriUPLO mTri) const {

  if ( M != N ) {
    printf("DMatrix::DTRSV. Error: Matrix is not square. Function call takes no effect.\n");
    return;
  }
  if ( N != B.GetM() ) {
    printf("DMatrix::DTRSV. Error: Matrix/vector sizes are not compatible. Function call takes no effect.\n");
    return;
  }

  char UPLO;
  switch (mTri) {
  case UPPER: UPLO = UPLO_U; break;
  case LOWER: UPLO = UPLO_L; break;
  }

  char TRANS;
  switch (ModeA) {
  case NORMAL: TRANS = 'N'; break;
  case TRANSPOSE: TRANS = 'T'; break;
  case CONJ_TRANS: TRANS = 'C'; break;
  }

  X = B;
  double *mX = X.GetPointer();
  INTEGER NRHS = B.GetN();

  dtrsm_(&SIDE_L, &UPLO, &TRANS, &DIAG_N, &N, &NRHS, ONE, A, &N, mX, &N);

}


//--------------------------------------------------------------------------
double DMatrix::
DPOCON(TriUPLO mTri) const {

  char UPLO;
  switch (mTri) {
  case UPPER: UPLO = UPLO_U; break;
  case LOWER: UPLO = UPLO_L; break;
  }

  double *WORK = NULL;
  New_1D_Array<double, INTEGER>(&WORK, 3*N);
  INTEGER *IWORK = NULL;
  New_1D_Array<INTEGER, INTEGER>(&IWORK, N);

  double RCOND = 0.0;
  INTEGER INFO = 0;

  double ANORM = dlange_(&NORM_1, &N, &N, A, &N, WORK);
  dpocon_(&UPLO, &N, A, &N, &ANORM, &RCOND, WORK, IWORK, &INFO);

  if (INFO < 0) {
    printf("DMatrix::DPOCON. DPOCON error: INFO = %ld. Computed result cannot be trusted.\n", (long)INFO);
  }

  Delete_1D_Array<double>(&WORK);
  Delete_1D_Array<INTEGER>(&IWORK);

  return RCOND;

}


//--------------------------------------------------------------------------
void DMatrix::
Chol(DMatrix &G, TriUPLO mTri) const {

  if (M != N) {
    printf("DMatrix::Chol. Error: Matrix is not square. Function call takes no effect.\n");
    return;
  }

  G = *this;
  double *mG = G.GetPointer();

  INTEGER INFO = 0;

  char UPLO;
  switch (mTri) {
  case UPPER: UPLO = UPLO_U; break;
  case LOWER: UPLO = UPLO_L; break;
  }

  dpotrf_(&UPLO, &N, mG, &N, &INFO);

  if (INFO != 0) {
    printf("DMatrix::Chol. DPOTRF error: INFO = %ld. Computed result cannot be trusted.\n", (long)INFO);
  }
  else {
    switch (mTri) {
    case UPPER: {
      for (INTEGER i = 0; i < N; i++) {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
        for (INTEGER j = 0; j < i; j++) {
          mG[i+j*N] = 0;
        }
      }
      break;
    }
    case LOWER: {
      for (INTEGER j = 0; j < N; j++) {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
        for (INTEGER i = 0; i < j; i++) {
          mG[i+j*N] = 0;
        }
      }
      break;
    }
    }
  }

}


//--------------------------------------------------------------------------
void DMatrix::
SymEig(DVector &lambda) const {
  DMatrix V;
  DSYEV(JOBZ_N, lambda, V);
}


//--------------------------------------------------------------------------
void DMatrix::
SymEig(DVector &lambda, DMatrix &V) const {
  DSYEV(JOBZ_V, lambda, V);
}


//--------------------------------------------------------------------------
void DMatrix::
DSYEV(char JOBZ, DVector &lambda, DMatrix &V) const {

  if (M != N) {
    printf("DMatrix::DSYEV. Error: Matrix is not square. Function call takes no effect.\n");
    return;
  }

  DMatrix A_copy = *this;
  double *mA = A_copy.GetPointer();

  double *W = NULL;
  New_1D_Array<double, INTEGER>(&W, N);

  INTEGER INFO = 0;
  double WORK0 = 0;
  INTEGER LWORK = -1;

  dsyev_(&JOBZ, &UPLO_L, &N, mA, &N, W, &WORK0, &LWORK, &INFO);

  LWORK = (INTEGER)WORK0;
  double *WORK = NULL;
  New_1D_Array<double, INTEGER>(&WORK, LWORK);

  dsyev_(&JOBZ, &UPLO_L, &N, mA, &N, W, WORK, &LWORK, &INFO);

  if (INFO != 0) {
    printf("DMatrix::DSYEV. DSYEV error: INFO = %ld. Symmetric eigen-decomposition is not computed.\n", (long)INFO);
  }
  else {
    lambda.Init(N);
    lambda.SetBlock(0, N, W);
    if (JOBZ == 'V') {
      V.Init(N);
      V.SetBlock(0, N, 0, N, mA);
    }
  }

  Delete_1D_Array<double>(&W);
  Delete_1D_Array<double>(&WORK);

}


//--------------------------------------------------------------------------
void DMatrix::
SymEigBIndef(const DMatrix &B, DVector &lambda) const {
  DMatrix V;
  DGGEV(JOBVR_N, B, lambda, V);
}


//--------------------------------------------------------------------------
void DMatrix::
SymEigBIndef(const DMatrix &B, DVector &lambda, DMatrix &V) const {
  DGGEV(JOBVR_V, B, lambda, V);
}


//--------------------------------------------------------------------------
void DMatrix::
DGGEV(char JOBVR, const DMatrix &B, DVector &lambda, DMatrix &V) const {

  if (M != N) {
    printf("DMatrix::DGGEV. Error: Matrix is not square. Function call takes no effect.\n");
    return;
  }
  if (M != B.GetM() || N != B.GetN()) {
    printf("DMatrix::DGGEV. Error: Matrices A and B do not have the same size. Function call takes no effect.\n");
    return;
  }

  DMatrix A_copy = *this;
  double *mA = A_copy.GetPointer();
  DMatrix B_copy = B;
  double *mB = B_copy.GetPointer();
  if (JOBVR == 'V') {
    V.Init(N);
  }
  double *mV = V.GetPointer();

  double *ALPHAR = NULL;
  New_1D_Array<double, INTEGER>(&ALPHAR, N);
  double *ALPHAI = NULL;
  New_1D_Array<double, INTEGER>(&ALPHAI, N);
  double *BETA = NULL;
  New_1D_Array<double, INTEGER>(&BETA, N);

  INTEGER INFO = 0;
  double WORK0 = 0;
  INTEGER LWORK = -1;

  dggev_(&JOBVL_N, &JOBVR, &N, mA, &N, mB, &N, ALPHAR, ALPHAI, BETA, NULL, &N, mV, &N, &WORK0, &LWORK, &INFO);

  LWORK = (INTEGER)WORK0;
  double *WORK = NULL;
  New_1D_Array<double, INTEGER>(&WORK, LWORK);

  dggev_(&JOBVL_N, &JOBVR, &N, mA, &N, mB, &N, ALPHAR, ALPHAI, BETA, NULL, &N, mV, &N, WORK, &LWORK, &INFO);

  if (INFO != 0) {
    printf("DMatrix::DGGEV. DGGEV error: INFO = %ld. Unsymmetric generalized eigen-decomposition is not computed.\n", (long)INFO);
  }
  else {

    // Results of DGGEV: For j = 0:N-1, the eigenvalues are
    // ALPHAR[j]/BETA[j] + ALPHAI[j]/BETA[j]*I. If ALPHAI[j] is zero,
    // this eigenvalue is real and the associated eigenvector is
    // mV(:,j). If ALPHAI[j] is nonzero, then the eigenvalues come in
    // conjugate pairs. Let's say the conjugate eigenvalue is in the
    // (j+1)-th position. Then ALPHAI[j] is positive and ALPHAI[j+1] =
    // -ALPHAI[j]. In other words, the conjugate eigenvalue is
    // ALPHAR[j]/BETA[j] - ALPHAI[j]/BETA[j]*I. In this case, the
    // associated eigenvectors also come in conjugate pairs and they
    // are mV(:,j) + mV(:,j+1)*I and mV(:,j) - mV(:,j+1)*I.
    //
    // Because the matrices A and B are symmetric, in theory all the
    // eigenvalues are real. If nonreal eigenvalues result from
    // numerical roundoff, and if the imaginary part
    // ALPHAI[j]/BETA[j]*I is small, we treat the eigenvalues to be
    // two same values ALPHAR[j]/BETA[j] and ALPHAR[j]/BETA[j] and the
    // associated eigenvectors to be mV(:,j) and mV(:,j+1).
    //
    // If we fix these eigenvectors, a better approximation of the
    // eigenvalues could be
    //
    //   ALPHAR[j]   ALPHAI[j]   dot( mV(:,j),mV(:,j+1) )
    //   --------- - --------- * ------------------------
    //    BETA[j]     BETA[j]     dot( mV(:,j),mV(:,j) )
    //
    // and
    //
    //   ALPHAR[j]   ALPHAI[j]    dot( mV(:,j),mV(:,j+1) )
    //   --------- + --------- * --------------------------
    //    BETA[j]     BETA[j]    dot( mV(:,j+1),mV(:,j+1) )
    //
    // because they minimize the residual equation. However, it is
    // perhaps not worth the effort to perform this computation and
    // hence we do not use these approximations.
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (INTEGER j = 0; j < N; j++) {
      if (ALPHAI[j]/BETA[j] > EPSx10 && ALPHAI[j]/ALPHAR[j] > EPSx10) {
#ifdef USE_OPENMP
#pragma omp critical
#endif
        printf("DMatrix::DGGEV. Warning: Spurious complex eigenvalue %g + %g * I. Treat complex part zero.\n", ALPHAR[j]/BETA[j], ALPHAI[j]/BETA[j]);
      }
      if (BETA[j] != 0.0) {
        ALPHAR[j] /= BETA[j];
      }
      else {
        ALPHAR[j] = INFINITY;
      }
    }
    lambda.Init(N);
    lambda.SetBlock(0, N, ALPHAR);

  }

  Delete_1D_Array<double>(&ALPHAR);
  Delete_1D_Array<double>(&ALPHAI);
  Delete_1D_Array<double>(&BETA);
  Delete_1D_Array<double>(&WORK);

}


//--------------------------------------------------------------------------
void DMatrix::
RealSchur(DMatrix &T, DMatrix &U, SelectType type) const {
  switch (type) {
  case LHP:
    DGEES(T, U, SelectLeftHalfPlane);
    break;
  case RHP:
    DGEES(T, U, SelectRightHalfPlane);
    break;
  }
}


//--------------------------------------------------------------------------
void DMatrix::
DGEES(DMatrix &T, DMatrix &U,
      LOGICAL (*Select)(double *WR, double *WI)) const {

  if (M != N) {
    printf("DMatrix::DGEES. Error: Matrix is not square. Function call takes no effect.\n");
    return;
  }

  T = *this;
  double *mT = T.GetPointer();
  U.Init(N);
  double *mU = U.GetPointer();

  double *WR = NULL;
  New_1D_Array<double, INTEGER>(&WR, N);
  double *WI = NULL;
  New_1D_Array<double, INTEGER>(&WI, N);
  LOGICAL *BWORK = NULL;
  New_1D_Array<LOGICAL, INTEGER>(&BWORK, N);

  INTEGER SDIM = -1;
  INTEGER INFO = 0;
  double WORK0 = 0;
  INTEGER LWORK = -1;

  dgees_(&JOBVS_V, &SORT_S, Select, &N, mT, &N, &SDIM, WR, WI, mU, &N, &WORK0, &LWORK, BWORK, &INFO);

  LWORK = (INTEGER)WORK0;
  double *WORK = NULL;
  New_1D_Array<double, INTEGER>(&WORK, LWORK);

  dgees_(&JOBVS_V, &SORT_S, Select, &N, mT, &N, &SDIM, WR, WI, mU, &N, WORK, &LWORK, BWORK, &INFO);

  if (INFO != 0) {
    printf("DMatrix::DGEES. DGEES error: INFO = %ld, N = %ld. Real Schur decomposition is not computed.\n", (long)INFO, (long)N);
  }

  Delete_1D_Array<double>(&WR);
  Delete_1D_Array<double>(&WI);
  Delete_1D_Array<double>(&WORK);
  Delete_1D_Array<LOGICAL>(&BWORK);

}


//--------------------------------------------------------------------------
void DMatrix::
Orth(DMatrix &B) const {

  if (M != N) {
    printf("DMatrix::Orth. Error: Matrix is not square. Function call takes no effect.\n");
    return;
  }

  B = *this;
  double *mB = B.GetPointer();

  double *TAU = NULL;
  New_1D_Array<double, INTEGER>(&TAU, N);

  INTEGER INFO = 0;
  double WORK0 = 0;
  INTEGER LWORK = -1;

  dgeqrf_(&N, &N, mB, &N, TAU, &WORK0, &LWORK, &INFO);

  LWORK = (INTEGER)WORK0;
  double *WORK = NULL;
  New_1D_Array<double, INTEGER>(&WORK, LWORK);

  dgeqrf_(&N, &N, mB, &N, TAU, WORK, &LWORK, &INFO);

  if (INFO != 0) {
    printf("DMatrix::Orth. DGEQRF error: INFO = %ld. QR factorziation is not computed.\n", (long)INFO);
    goto END;
  }

  dorgqr_(&N, &N, &N, mB, &N, TAU, WORK, &LWORK, &INFO);

  if (INFO != 0) {
    printf("DMatrix::Orth. DORGQR error: INFO = %ld. Orthogonal factor is not computed.\n", (long)INFO);
  }

 END:
  Delete_1D_Array<double>(&TAU);
  Delete_1D_Array<double>(&WORK);

}


//--------------------------------------------------------------------------
void DMatrix::
QRpivots(INTEGER k, INTEGER *pivots) const {

  INTEGER minMN = M<N?M:N;
  if (k > minMN) {
    printf("DMatrix::QRpivots. Error: k cannot be larger than the smallest dimension of the matrix. Function call takes no effect.\n");
    return;
  }

  DMatrix A_copy = *this;
  double *mA = A_copy.GetPointer();

  INTEGER *JPVT = NULL;
  New_1D_Array<INTEGER, INTEGER>(&JPVT, N); // Already init 0
  double *TAU = NULL;
  New_1D_Array<double, INTEGER>(&TAU, minMN);

  double WORK0 = 0;
  INTEGER LWORK = -1;
  INTEGER INFO = 0;

  dgeqp3_(&M, &N, mA, &M, JPVT, TAU, &WORK0, &LWORK, &INFO);

  LWORK = (INTEGER)WORK0;
  double *WORK = NULL;
  New_1D_Array<double, INTEGER>(&WORK, LWORK);

  dgeqp3_(&M, &N, mA, &M, JPVT, TAU, WORK, &LWORK, &INFO);

  if (INFO != 0) {
    printf("DMatrix::QRpivots. DGEQP3 error: INFO = %ld. Computed result cannot be trusted.\n", (long)INFO);
  }

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < k; i++) {
    pivots[i] = JPVT[i] - 1;
  }

  Delete_1D_Array<INTEGER>(&JPVT);
  Delete_1D_Array<double>(&TAU);
  Delete_1D_Array<double>(&WORK);

}


// TODO: So far here
//--------------------------------------------------------------------------
INTEGER DMatrix::
Rank(void) const {

  INTEGER minMN = M<N?M:N;
  double *S = NULL;
  New_1D_Array<double, INTEGER>(&S, minMN);
  DGESVD(S);

  if (S == NULL) {
    printf("DMatrix::Rank. Error: Return -1 due to DGESVD error.\n");
    return -1;
  }
  
  INTEGER maxMN = M>N?M:N;
  double smallest = (double)maxMN * S[0] * EPS;
  INTEGER mRank = 0;
  for (INTEGER i = 0; i < minMN; i++) {
    if (S[i] > smallest) {
      mRank++;
    }
  }
  Delete_1D_Array<double>(&S);
  return mRank;
}


//--------------------------------------------------------------------------
double DMatrix::
Cond2(void) const {

  INTEGER minMN = M<N?M:N;
  double *S = NULL;
  New_1D_Array<double, INTEGER>(&S, minMN);
  DGESVD(S);

  if (S == NULL) {
    printf("DMatrix::Rank. Error: Return NAN due to DGESVD error.\n");
    return NAN;
  }
  
  double mCond2 = S[0] / S[minMN-1];
  Delete_1D_Array<double>(&S);
  return mCond2;
}


//--------------------------------------------------------------------------
double DMatrix::
Norm2(void) const {

  INTEGER minMN = M<N?M:N;
  double *S = NULL;
  New_1D_Array<double, INTEGER>(&S, minMN);
  DGESVD(S);

  if (S == NULL) {
    printf("DMatrix::Rank. Error: Return NAN due to DGESVD error.\n");
    return NAN;
  }
  
  double mNorm2 = S[0];
  Delete_1D_Array<double>(&S);
  return mNorm2;
}


//--------------------------------------------------------------------------
void DMatrix::
DGESVD(double *S) const {

  DMatrix A_copy = *this;
  double *mA = A_copy.GetPointer();

  INTEGER minMN = M<N?M:N;
  INTEGER INFO = 0;
  double WORK0 = 0;
  INTEGER LWORK = -1;

  dgesvd_(&JOBU_N, &JOBVT_N, &M, &N, mA, &M, S, NULL, &M, NULL, &minMN,
          &WORK0, &LWORK, &INFO);

  LWORK = (INTEGER)WORK0;
  double *WORK = NULL;
  New_1D_Array<double, INTEGER>(&WORK, LWORK);

  dgesvd_(&JOBU_N, &JOBVT_N, &M, &N, mA, &M, S, NULL, &M, NULL, &minMN,
          WORK, &LWORK, &INFO);

  if (INFO != 0) {
    printf("DMatrix::DGESVD. DGESVD error: INFO = %ld. Computed result cannot be trusted.\n", (long)INFO);
  }

  Delete_1D_Array<double>(&WORK);

}


//--------------------------------------------------------------------------
double DMatrix::
NormF(void) const {
  INTEGER MN = M*N;
  double mNormF = dnrm2_(&MN, A, &ONEi);
  return mNormF;
}


//--------------------------------------------------------------------------
void DMatrix::
Diag(DVector &b) const {
  INTEGER minMN = M<N ? M:N;
  b.Init(minMN);
  double *mb = b.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER j = 0; j < minMN; j++) {
    mb[j] = A[M*j+j];
  }
}


//--------------------------------------------------------------------------
double DMatrix::
Tr(void) const {
  INTEGER minMN = M<N ? M:N;
  double mTr = 0.0;
#ifdef USE_OPENMP
#pragma omp parallel for reduction(+:mTr)
#endif
  for (INTEGER j = 0; j < minMN; j++) {
    mTr += A[M*j+j];
  }
  return mTr;
}


//--------------------------------------------------------------------------
LogDet DMatrix::
Det(MatrixType TypeA, MatState mState) {

  LogDet mLogDet = { 0.0, 1 };

  if ( M != N ) {
    printf("DMatrix::Det. Error: Matrix is not square. Function call takes no effect.\n");
    return mLogDet;
  }
  if ((TypeA == GENERAL && mState == CHOL_FACT) ||
      (TypeA == SPD && mState == LU_FACT)) {
    printf("DMatrix::Det. Error: Invalid combination of matrix type and matrix state. Function call takes no effect.\n");
    return mLogDet;
  }

  double *mA = NULL;
  DMatrix A_copy;
  switch (TypeA) {
  case GENERAL: {

    if (mState == UNFACT) {

      A_copy = *this;
      mA = A_copy.GetPointer();
      INTEGER INFO = 0;

      Delete_1D_Array<INTEGER>(&IPIV);
      New_1D_Array<INTEGER, INTEGER>(&IPIV, N);

      dgetrf_(&N, &N, mA, &N, IPIV, &INFO);

      if (INFO < 0) {
        printf("DMatrix:Det. DGETRF error: INFO = %ld. Computed result cannot be trusted.\n", (long)INFO);
        Delete_1D_Array<INTEGER>(&IPIV);
        return mLogDet;
      }
      else if (INFO > 0) {
        printf("DMatrix::Det. Matrix is singular and the log determinant is -INFINITY.\n");
        Delete_1D_Array<INTEGER>(&IPIV);
        mLogDet.LogAbsDet = -INFINITY;
        return mLogDet;
      }

    }
    else if (mState == LU_FACT) {
      mA = A;
    }

    for (INTEGER j = 0; j < N; j++) {
      mLogDet.LogAbsDet += log(fabs(mA[N*j+j]));
      if (mA[N*j+j] < 0) {
        mLogDet.Sign *= -1;
      }
    }
    INTEGER par = Parity(IPIV, N);
    if (par == 1) {
      mLogDet.Sign = -mLogDet.Sign;
    }

    if (mState == UNFACT) {
      Delete_1D_Array<INTEGER>(&IPIV);
    }

    break;
  }
  case SPD: {

    if (mState == UNFACT) {
      
      A_copy = *this;
      mA = A_copy.GetPointer();
      INTEGER INFO = 0;

      dpotrf_(&UPLO_L, &N, mA, &N, &INFO);

      if (INFO < 0) {
        printf("DMatrix:Det. DPOTRF error: INFO = %ld. Computed result cannot be trusted.\n", (long)INFO);
        return mLogDet;
      }
      else if (INFO > 0) {
        printf("DMatrix::Det. DPOTRF error: INFO = %ld. Matrix is not positive definite. Computed result cannot be trusted.\n", (long)INFO);
        return mLogDet;
      }

    }
    else if (mState == CHOL_FACT) {
      mA = A;
    }

    for (INTEGER j = 0; j < N; j++) {
      mLogDet.LogAbsDet += log(mA[N*j+j]);
    }
    mLogDet.LogAbsDet *= 2.0;

    break;
  }

  }

  return mLogDet;

}
