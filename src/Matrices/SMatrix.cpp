#include "SMatrix.hpp"

#define INITVAL_M   0
#define INITVAL_N   0
#define INITVAL_nnz 0


//--------------------------------------------------------------------------
SMatrix::
SMatrix() {
  start = NULL;
  idx = NULL;
  A = NULL;
  M = INITVAL_M;
  N = INITVAL_N;
  nnz = INITVAL_nnz;
  Init();
}


//--------------------------------------------------------------------------
void SMatrix::
Init(void) {
  Init(0,0,0);
}


//--------------------------------------------------------------------------
void SMatrix::
Init(INTEGER N_, INTEGER nnz_) {
  Init(N_, N_, nnz_);
}


//--------------------------------------------------------------------------
void SMatrix::
Init(INTEGER M_, INTEGER N_, INTEGER nnz_) {
  if (M != M_ || nnz != nnz_) {
    ReleaseAllMemory();
  }
  if (!start) {
    New_1D_Array<INTEGER, INTEGER>(&start, M_+1);
  }
  if (!idx) {
    New_1D_Array<INTEGER, INTEGER>(&idx, nnz_);
  }
  if (!A) {
    New_1D_Array<double, INTEGER>(&A, nnz_);
  }
  M = M_;
  N = N_;
  nnz = nnz_;
}


//--------------------------------------------------------------------------
void SMatrix::
ReleaseAllMemory(void) {
  M = INITVAL_M;
  N = INITVAL_N;
  nnz = INITVAL_nnz;
  Delete_1D_Array<INTEGER>(&start);
  Delete_1D_Array<INTEGER>(&idx);
  Delete_1D_Array<double>(&A);
}


//--------------------------------------------------------------------------
SMatrix::
SMatrix(INTEGER N_, INTEGER nnz_) {
  start = NULL;
  idx = NULL;
  A = NULL;
  M = INITVAL_M;
  N = INITVAL_N;
  nnz = INITVAL_nnz;
  Init(N_, nnz_);
}


//--------------------------------------------------------------------------
SMatrix::
SMatrix(INTEGER M_, INTEGER N_, INTEGER nnz_) {
  start = NULL;
  idx = NULL;
  A = NULL;
  M = INITVAL_M;
  N = INITVAL_N;
  nnz = INITVAL_nnz;
  Init(M_, N_, nnz_);
}


//--------------------------------------------------------------------------
SMatrix::
SMatrix(const SMatrix &G) {
  M = INITVAL_M;
  N = INITVAL_N;
  nnz = INITVAL_nnz;
  start = NULL;
  idx = NULL;
  A = NULL;
  Init();
  DeepCopy(G);
}


//--------------------------------------------------------------------------
SMatrix& SMatrix::
operator= (const SMatrix &G) {
  if (this != &G) {
    DeepCopy(G);
  }
  return *this;
}


//--------------------------------------------------------------------------
void SMatrix::
DeepCopy(const SMatrix &G) {
  if (M != G.M || nnz != G.nnz) {
    ReleaseAllMemory();
  }
  if (G.start) {
    if (!start) {
      New_1D_Array<INTEGER, INTEGER>(&start, G.M+1);
    }
    memcpy(start, G.start, (G.M+1)*sizeof(INTEGER));
  }
  if (G.idx) {
    if (!idx) {
      New_1D_Array<INTEGER, INTEGER>(&idx, G.nnz);
    }
    memcpy(idx, G.idx, G.nnz*sizeof(INTEGER));
  }
  if (G.A) {
    if (!A) {
      New_1D_Array<double, INTEGER>(&A, G.nnz);
    }
    memcpy(A, G.A, G.nnz*sizeof(double));
  }
  M = G.M;
  N = G.N;
  nnz = G.nnz;
}


//--------------------------------------------------------------------------
SMatrix::
~SMatrix() {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
INTEGER SMatrix::
GetM(void) const {
  return M;
}


//--------------------------------------------------------------------------
INTEGER SMatrix::
GetN(void) const {
  return N;
}


//--------------------------------------------------------------------------
INTEGER SMatrix::
GetNNZ(void) const {
  return nnz;
}


//--------------------------------------------------------------------------
INTEGER* SMatrix::
GetPointerStart(void) const {
  return start;
}


//--------------------------------------------------------------------------
INTEGER* SMatrix::
GetPointerIdx(void) const {
  return idx;
}


//--------------------------------------------------------------------------
double* SMatrix::
GetPointerA(void) const {
  return A;
}


//--------------------------------------------------------------------------
void SMatrix::
GetBlock(INTEGER RowStart, INTEGER nRow, INTEGER ColStart, INTEGER nCol,
         DMatrix &B) const {
  if (RowStart < 0 || RowStart+nRow > M || ColStart < 0 || ColStart+nCol > N) {
    printf("SMatrix::GetBlock. Error: Invalid range. Function call takes no effect.\n");
    return;
  }
  B.Init(nRow, nCol);
  double *mB = B.GetPointer(); // Initialized with 0
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < nRow; i++) {
    INTEGER ii = i + RowStart;
    for (INTEGER j = start[ii]; j < start[ii+1]; j++) {
      // mB[i, idx[j]-ColStart] <- A[j]
      INTEGER jj = idx[j] - ColStart;
      if (jj >= 0 && jj < nCol) {
        mB[i+jj*nRow] = A[j]; // Note that mB is in column-major order
      }
    }
  }
}


//--------------------------------------------------------------------------
void SMatrix::
GetBlock(INTEGER *IdxRow, INTEGER nRow, INTEGER *IdxCol, INTEGER nCol,
         DMatrix &B) const {
  B.Init(nRow, nCol);
  double *mB = B.GetPointer(); // Does not matter whether it's initialized 0
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < nRow; i++) {
    INTEGER ii = IdxRow[i];
    for (INTEGER j = 0; j < nCol; j++) {
      INTEGER *base = idx + start[ii];
      INTEGER *ptr = NULL;
      ptr = (INTEGER *)bsearch(&(IdxCol[j]), base, start[ii+1]-start[ii],
                               sizeof(INTEGER), CompareIntegerNaturalOrderLess);
      if (ptr != NULL) {
        // mB[i,j] <- A[loc+start[ii]].
        // Note that mB is in column-major order
        INTEGER loc = ptr - base;
        mB[i+j*nRow] = A[loc+start[ii]];
      }
      else {
        mB[i+j*nRow] = 0;
      }
    }
  }
}


//--------------------------------------------------------------------------
void SMatrix::
GetBlock(INTEGER *IdxRow, INTEGER nRow, INTEGER ColStart, INTEGER nCol,
         DMatrix &B) const {
  if (ColStart < 0 || ColStart+nCol > N) {
    printf("SMatrix::GetBlock. Error: Invalid range. Function call takes no effect.\n");
    return;
  }
  B.Init(nRow, nCol);
  double *mB = B.GetPointer(); // Initialized with 0
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < nRow; i++) {
    INTEGER ii = IdxRow[i];
    for (INTEGER j = start[ii]; j < start[ii+1]; j++) {
      // mB[i, idx[j]-ColStart] <- A[j]
      INTEGER jj = idx[j] - ColStart;
      if (jj >= 0 && jj < nCol) {
        mB[i+jj*nRow] = A[j]; // Note that mB is in column-major order
      }
    }
  }
}


//--------------------------------------------------------------------------
void SMatrix::
Divide(double b) {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < nnz; i++) {
    A[i] /= b;
  }
}


//--------------------------------------------------------------------------
void SMatrix::
Divide(double b, SMatrix &C) const {
  C = *this;
  double *mC = C.GetPointerA();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < nnz; i++) {
    mC[i] /= b;
  }
}


//--------------------------------------------------------------------------
void SMatrix::
SymmetricDivide(const DVector &b) {
  if (M != N) {
    printf("SMatrix::SymmetricDivide. Error: Matrix is not square. Function call takes no effect.\n");
    return;
  }
  if (b.GetN() != N) {
    printf("SMatrix::SymmetricDivide. Error: Vector length not compatible. Function call takes no effect.\n");
    return;
  }
  double *mb = b.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < N; i++) {
    for (INTEGER jj = start[i]; jj < start[i+1]; jj++) {
      INTEGER j = idx[jj];
      A[jj] /= (mb[i]*mb[j]);
    }
  }
}


//--------------------------------------------------------------------------
void SMatrix::
SymmetricDivide(const DVector &b, SMatrix &C) const {
  if (M != N) {
    printf("SMatrix::SymmetricDivide. Error: Matrix is not square. Function call takes no effect.\n");
    return;
  }
  if (b.GetN() != N) {
    printf("SMatrix::SymmetricDivide. Error: Vector length not compatible. Function call takes no effect.\n");
    return;
  }
  C = *this;
  double *mC = C.GetPointerA();
  double *mb = b.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < N; i++) {
    for (INTEGER jj = start[i]; jj < start[i+1]; jj++) {
      INTEGER j = idx[jj];
      mC[jj] /= (mb[i]*mb[j]);
    }
  }
}


//--------------------------------------------------------------------------
void SMatrix::
MatVec(const DVector &b, DVector &y, MatrixMode ModeA, MatState mState) const {

  if ( (ModeA == NORMAL && N != b.GetN()) ||
       (ModeA != NORMAL && M != b.GetN()) ) {
    printf("SMatrix::MatVec. Error: Matrix/vector sizes are not compatible. Function call takes no effect.\n");
    return;
  }

  if (mState != UNFACT) {
    printf("SMatrix::MatVec. Error: mState must be UNFACT. Function call takes no effect.\n");
    return;
  }

  char TRANS = 0;
  switch (ModeA) {
  case NORMAL: TRANS = 'N'; y.Init(M); break;
  case TRANSPOSE: TRANS = 'T'; y.Init(N); break;
  case CONJ_TRANS: TRANS = 'C'; y.Init(N); break;
  }

  double *mb = b.GetPointer();
  double *my = y.GetPointer();

  sp_d_dgemv(TRANS, M, N, nnz, start, idx, A, mb, my);

}


//--------------------------------------------------------------------------
void SMatrix::
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

  if ( (ModeA == NORMAL && ModeB == NORMAL &&
        (N != B.GetM() || C.GetM() != M || C.GetN() != B.GetN())) ||
       (ModeA == NORMAL && ModeB != NORMAL &&
        (N != B.GetN() || C.GetM() != M || C.GetN() != B.GetM())) ||
       (ModeA != NORMAL && ModeB == NORMAL &&
        (M != B.GetM() || C.GetM() != N || C.GetN() != B.GetN())) ||
       (ModeA != NORMAL && ModeB != NORMAL &&
        (M != B.GetN() || C.GetM() != N || C.GetN() != B.GetM())) ) {
    printf("SMatrix::MatMat. Error: Matrix sizes are not compatible. Function call takes no effect.\n");
    return;
  }

  INTEGER mM = 0, mN = 0, mK = 0;
  char TRANSA = 0;
  switch (ModeA) {
  case NORMAL: TRANSA = 'N'; mM = M; mK = N; break;
  case TRANSPOSE: TRANSA = 'T'; mM = N; mK = M; break;
  case CONJ_TRANS: TRANSA = 'C'; mM = N; mK = M; break;
  }
  char TRANSB = 0;
  switch (ModeB) {
  case NORMAL: TRANSB = 'N'; mN = B.GetN(); break;
  case TRANSPOSE: TRANSB = 'T'; mN = B.GetM(); break;
  case CONJ_TRANS: TRANSB = 'C'; mN = B.GetM(); break;
  }

  double *mB = B.GetPointer();
  double *mC = C.GetPointer();

  sp_d_dgemm(TRANSA, TRANSB, mM, mN, mK, nnz, start, idx, A, mB, mC);

}


//--------------------------------------------------------------------------
void SMatrix::
MatMat(const SMatrix &B, DMatrix &C, MatrixMode ModeA, MatrixMode ModeB) const {

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

  if ( (ModeA == NORMAL && ModeB == NORMAL &&
        (N != B.GetM() || C.GetM() != M || C.GetN() != B.GetN())) ||
       (ModeA == NORMAL && ModeB != NORMAL &&
        (N != B.GetN() || C.GetM() != M || C.GetN() != B.GetM())) ||
       (ModeA != NORMAL && ModeB == NORMAL &&
        (M != B.GetM() || C.GetM() != N || C.GetN() != B.GetN())) ||
       (ModeA != NORMAL && ModeB != NORMAL &&
        (M != B.GetN() || C.GetM() != N || C.GetN() != B.GetM())) ) {
    printf("SMatrix::MatMat. Error: Matrix sizes are not compatible. Function call takes no effect.\n");
    return;
  }

  INTEGER mM = 0, mN = 0, mK = 0;
  char TRANSA = 0;
  switch (ModeA) {
  case NORMAL: TRANSA = 'N'; mM = M; mK = N; break;
  case TRANSPOSE: TRANSA = 'T'; mM = N; mK = M; break;
  case CONJ_TRANS: TRANSA = 'C'; mM = N; mK = M; break;
  }
  char TRANSB = 0;
  switch (ModeB) {
  case NORMAL: TRANSB = 'N'; mN = B.GetN(); break;
  case TRANSPOSE: TRANSB = 'T'; mN = B.GetM(); break;
  case CONJ_TRANS: TRANSB = 'C'; mN = B.GetM(); break;
  }

  INTEGER nnzB = B.GetNNZ();
  INTEGER *startB = B.GetPointerStart();
  INTEGER *idxB = B.GetPointerIdx();
  double *mB = B.GetPointerA();
  double *mC = C.GetPointer();

  sp_s_dgemm(TRANSA, TRANSB, mM, mN, mK, nnz, start, idx, A, nnzB, startB, idxB, mB, mC);

}


//--------------------------------------------------------------------------
void SMatrix::
Diag(DVector &b) const {
  b.Init(N);
  double *mb = b.GetPointer();
  Diag(mb);
}


//--------------------------------------------------------------------------
void SMatrix::
Diag(double *b) const {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < N; i++) {
    INTEGER *loc = NULL;
    loc = (INTEGER *)bsearch(&i, idx+start[i], start[i+1]-start[i],
                             sizeof(INTEGER), CompareIntegerNaturalOrderLess);
    if (loc == NULL) {
      b[i] = 0.0;
    }
    else {
      b[i] = A[loc-idx];
    }
  }
}
