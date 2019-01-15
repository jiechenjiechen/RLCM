#include "SPointArray.hpp"

#define INITVAL_N   0
#define INITVAL_d   0
#define INITVAL_nnz 0


//--------------------------------------------------------------------------
SPointArray::
SPointArray() {
  start = NULL;
  idx = NULL;
  X = NULL;
  N = INITVAL_N;
  d = INITVAL_d;
  nnz = INITVAL_nnz;
  Init();
}


//--------------------------------------------------------------------------
void SPointArray::
Init(void) {
  Init(0,0,0);
}


//--------------------------------------------------------------------------
void SPointArray::
Init(INTEGER N_, INTEGER d_, INTEGER nnz_) {
  if (N != N_ || nnz != nnz_) {
    ReleaseAllMemory();
  }
  if (!start) {
    New_1D_Array<INTEGER, INTEGER>(&start, N_+1);
  }
  if (!idx) {
    New_1D_Array<INTEGER, INTEGER>(&idx, nnz_);
  }
  if (!X) {
    New_1D_Array<double, INTEGER>(&X, nnz_);
  }
  N = N_;
  d = d_;
  nnz = nnz_;
}


//--------------------------------------------------------------------------
void SPointArray::
ReleaseAllMemory(void) {
  N = INITVAL_N;
  d = INITVAL_d;
  nnz = INITVAL_nnz;
  Delete_1D_Array<INTEGER>(&start);
  Delete_1D_Array<INTEGER>(&idx);
  Delete_1D_Array<double>(&X);
}


//--------------------------------------------------------------------------
SPointArray::
SPointArray(INTEGER N_, INTEGER d_, INTEGER nnz_) {
  start = NULL;
  idx = NULL;
  X = NULL;
  N = INITVAL_N;
  d = INITVAL_d;
  nnz = INITVAL_nnz;
  Init(N_, d_, nnz_);
}


//--------------------------------------------------------------------------
SPointArray::
SPointArray(const SPointArray& G) {
  N = INITVAL_N;
  d = INITVAL_d;
  nnz = INITVAL_nnz;
  start = NULL;
  idx = NULL;
  X = NULL;
  Init();
  DeepCopy(G);
}


//--------------------------------------------------------------------------
SPointArray& SPointArray::
operator= (const SPointArray& G) {
  if (this != &G) {
    DeepCopy(G);
  }
  return *this;
}


//--------------------------------------------------------------------------
SPointArray& SPointArray::
operator= (const DPointArray& G) {
  ReleaseAllMemory();

  N = G.GetN();
  d = G.GetD();

  nnz = 0;
  double *GX = G.GetPointer();
  for (INTEGER i = 0; i < N*d; i++) {
    if (GX[i] != 0.0) {
      nnz++;
    }
  }
  
  New_1D_Array<INTEGER, INTEGER>(&start, N+1);
  New_1D_Array<INTEGER, INTEGER>(&idx, nnz);
  New_1D_Array<double, INTEGER>(&X, nnz);
  INTEGER k = 0;
  double val;
  for (INTEGER i = 0; i < N; i++) {
    start[i] = k;
    for (INTEGER j = 0; j < d; j++) {
      val = GX[i+j*N];
      if (val != 0.0) {
        idx[k] = j;
        X[k] = val;
        k++;
      }
    }
  }
  start[N] = k;

  return *this;
}


//--------------------------------------------------------------------------
void SPointArray::
DeepCopy(const SPointArray &G) {
  if (N != G.N || nnz != G.nnz) {
    ReleaseAllMemory();
  }
  if (G.start) {
    if (!start) {
      New_1D_Array<INTEGER, INTEGER>(&start, G.N+1);
    }
    memcpy(start, G.start, (G.N+1)*sizeof(INTEGER));
  }
  if (G.idx) {
    if (!idx) {
      New_1D_Array<INTEGER, INTEGER>(&idx, G.nnz);
    }
    memcpy(idx, G.idx, G.nnz*sizeof(INTEGER));
  }
  if (G.X) {
    if (!X) {
      New_1D_Array<double, INTEGER>(&X, G.nnz);
    }
    memcpy(X, G.X, G.nnz*sizeof(double));
  }
  N = G.N;
  d = G.d;
  nnz = G.nnz;
}


//--------------------------------------------------------------------------
SPointArray::
~SPointArray() {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
INTEGER SPointArray::
GetD(void) const {
  return d;
}


//--------------------------------------------------------------------------
INTEGER SPointArray::
GetN(void) const {
  return N;
}


//--------------------------------------------------------------------------
INTEGER SPointArray::
GetNNZ(void) const {
  return nnz;
}


//--------------------------------------------------------------------------
void SPointArray::
GetPoint(INTEGER i, SPoint &x) const {
  if (i < 0 || i >= N) {
    printf("SPointArray::GetPoint. Error: Invalid index. Function call takes no effect.\n");
    return;
  }
  INTEGER mnnz = start[i+1] - start[i];
  x.Init(d, mnnz);
  memcpy(x.GetPointerIdx(), idx+start[i], mnnz*sizeof(INTEGER));
  memcpy(x.GetPointerX(), X+start[i], mnnz*sizeof(double));
}


//--------------------------------------------------------------------------
void SPointArray::
GetPoint(INTEGER i, DPoint &x) const {
  if (i < 0 || i >= N) {
    printf("SPointArray::GetPoint. Error: Invalid index. Function call takes no effect.\n");
    return;
  }
  x.Init(d); // Initialized with 0
  double *mx = x.GetPointer();
  for (INTEGER j = start[i]; j < start[i+1]; j++) {
    mx[idx[j]] = X[j];
  }
}


//--------------------------------------------------------------------------
void SPointArray::
GetSubset(INTEGER istart, INTEGER n, SPointArray &Y) const {
  if (istart < 0 || istart+n > N) {
    printf("SPointArray::GetSubset. Error: Invalid range. Function call takes no effect.\n");
    return;
  }
  INTEGER nnz_Y = start[istart+n] - start[istart];
  Y.Init(n, d, nnz_Y);
  memcpy(Y.GetPointerIdx(), idx+start[istart], nnz_Y*sizeof(INTEGER));
  memcpy(Y.GetPointerX(), X+start[istart], nnz_Y*sizeof(double));

  INTEGER *start_Y = Y.GetPointerStart();
  memcpy(start_Y, start+istart, (n+1)*sizeof(INTEGER));
  INTEGER shift = start_Y[0];
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER j = 0; j <= n; j++) {
    start_Y[j] -= shift;
  }
}


//--------------------------------------------------------------------------
void SPointArray::
GetSubset(INTEGER istart, INTEGER n, DPointArray &Y) const {
  if (istart < 0 || istart+n > N) {
    printf("SPointArray::GetSubset. Error: Invalid range. Function call takes no effect.\n");
    return;
  }
  Y.Init(n,d); // Initialized with 0
  double *mY = Y.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < n; i++) {
    INTEGER ii = i + istart;
    for (INTEGER j = start[ii]; j < start[ii+1]; j++) {
      // mY[i, idx[j]] <- X[j]
      mY[i+idx[j]*n] = X[j]; // Note that mY is in column-major order
    }
  }
}


//--------------------------------------------------------------------------
void SPointArray::
GetSubset(INTEGER *iidx, INTEGER n, SPointArray &Y) const {
  // Get nnz for each requested point and count the total nnz
  INTEGER *len = NULL;
  New_1D_Array<INTEGER, INTEGER>(&len, n);
  INTEGER nnz_Y = 0;
#ifdef USE_OPENMP
#pragma omp parallel for reduction(+:nnz_Y)
#endif
  for (INTEGER i = 0; i < n; i++) {
    len[i] = start[iidx[i]+1] - start[iidx[i]];
    nnz_Y += len[i];
  }

  // Init Y
  Y.Init(n, d, nnz_Y);

  // Copy data
  INTEGER *start_Y = Y.GetPointerStart();
  INTEGER *idx_Y = Y.GetPointerIdx();
  double *my = Y.GetPointerX();
  start_Y[0] = 0;
  for (INTEGER i = 0; i < n; i++) {
    start_Y[i+1] = start_Y[i] + len[i]; // Cannot do OPENMP here
    memcpy(idx_Y+start_Y[i], idx+start[iidx[i]], len[i]*sizeof(INTEGER));
    memcpy(my+start_Y[i], X+start[iidx[i]], len[i]*sizeof(double));
  }

  // Clean up
  Delete_1D_Array<INTEGER>(&len);
}


//--------------------------------------------------------------------------
void SPointArray::
GetSubset(INTEGER *iidx, INTEGER n, DPointArray &Y) const {
  Y.Init(n,d); // Initialized with 0
  double *mY = Y.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < n; i++) {
    INTEGER ii = iidx[i];
    for (INTEGER j = start[ii]; j < start[ii+1]; j++) {
      // mY[i, idx[j]] <- X[j]
      mY[i+idx[j]*n] = X[j]; // Note that mY is in column-major order
    }
  }
}


//--------------------------------------------------------------------------
INTEGER* SPointArray::
GetPointerStart(void) const {
  return start;
}


//--------------------------------------------------------------------------
INTEGER* SPointArray::
GetPointerIdx(void) const {
  return idx;
}


//--------------------------------------------------------------------------
double* SPointArray::
GetPointerX(void) const {
  return X;
}


//--------------------------------------------------------------------------
void SPointArray::
AsDMatrix(DMatrix &A) const {
  A.Init(N,d); // Initialized with 0
  double *mA = A.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < N; i++) {
    for (INTEGER j = start[i]; j < start[i+1]; j++) {
      // mA[i, idx[j]] <- X[j]
      mA[i+idx[j]*N] = X[j]; // Note that mA is in column-major order
    }
  }
}


//--------------------------------------------------------------------------
void SPointArray::
PrintPointArray(const char *name) const {
  printf("%s:\n", name);
  for (INTEGER i = 0; i < N; i++) {
    for (INTEGER j = start[i]; j < start[i+1]; j++) {
      printf("%ld:%g ", (long)idx[j], X[j]);
    }
    printf("\n");
  }
  printf("\n");
}


//--------------------------------------------------------------------------
void SPointArray::
Center(DPoint &c) const {
  c.Init(d);
  double *mc = c.GetPointer();
  INTEGER j = 0;
  for (j = 0; j < nnz; j++) {
    mc[idx[j]] += X[j];
  }
  for (j = 0; j < d; j++) {
    mc[j] /= N;
  }
}


//--------------------------------------------------------------------------
double SPointArray::
StdvDist(void) const {
  DPoint c;
  Center(c);
  double dist2 = 0.0;
#ifdef USE_OPENMP
#pragma omp parallel for reduction(+:dist2)
#endif
  for (INTEGER i = 0; i < N; i++) {
    SPoint x;
    GetPoint(i, x);
    dist2 += x.Dist2(c);
  }
  double ret = sqrt(dist2/N);
  return ret;
}


//--------------------------------------------------------------------------
double SPointArray::
MaxPointNorm2(void) const {
  DVector PointNorm2(N);
  double *mPointNorm2 = PointNorm2.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < N; i++) {
    SPoint x;
    GetPoint(i, x);
    mPointNorm2[i] = x.InProd(x);
  }
  double ret = PointNorm2.Max();
  return ret;
}


//--------------------------------------------------------------------------
void SPointArray::
MatVec(const DVector &b, DVector &y, MatrixMode ModeX) const {

  if ( (ModeX == NORMAL && d != b.GetN()) ||
       (ModeX != NORMAL && N != b.GetN()) ) {
    printf("SPointArray::MatVec. Error: Matrix/vector sizes are not compatible. Function call takes no effect.\n");
    return;
  }

  char TRANS = 0;
  switch (ModeX) {
  case NORMAL: TRANS = 'N'; y.Init(N); break;
  case TRANSPOSE: TRANS = 'T'; y.Init(d); break;
  case CONJ_TRANS: TRANS = 'C'; y.Init(d); break;
  }

  double *mb = b.GetPointer();
  double *my = y.GetPointer();

  sp_d_dgemv(TRANS, N, d, nnz, start, idx, X, mb, my);

}


//--------------------------------------------------------------------------
void SPointArray::
MatVec(const DPoint &b, DVector &y, MatrixMode ModeX) const {

  if ( (ModeX == NORMAL && d != b.GetD()) ||
       (ModeX != NORMAL && N != b.GetD()) ) {
    printf("SPointArray::MatVec. Error: Matrix/vector sizes are not compatible. Function call takes no effect.\n");
    return;
  }

  char TRANS = 0;
  switch (ModeX) {
  case NORMAL: TRANS = 'N'; y.Init(N); break;
  case TRANSPOSE: TRANS = 'T'; y.Init(d); break;
  case CONJ_TRANS: TRANS = 'C'; y.Init(d); break;
  }

  double *mb = b.GetPointer();
  double *my = y.GetPointer();

  sp_d_dgemv(TRANS, N, d, nnz, start, idx, X, mb, my);

}


//--------------------------------------------------------------------------
void SPointArray::
MatMat(const DMatrix &B, DMatrix &Y, MatrixMode ModeX, MatrixMode ModeB) const {

  if ( (ModeX == NORMAL && ModeB == NORMAL && d != B.GetM()) ||
       (ModeX != NORMAL && ModeB == NORMAL && N != B.GetM()) ||
       (ModeX == NORMAL && ModeB != NORMAL && d != B.GetN()) ||
       (ModeX != NORMAL && ModeB != NORMAL && N != B.GetN()) ) {
    printf("SPointArray::MatMat. Error: Matrix sizes are not compatible. Function call takes no effect.\n");
    return;
  }

  INTEGER mM = 0, mN = 0, mK = 0;
  char TRANSX = 0;
  switch (ModeX) {
  case NORMAL: TRANSX = 'N'; mM = N; mK = d; break;
  case TRANSPOSE: TRANSX = 'T'; mM = d; mK = N; break;
  case CONJ_TRANS: TRANSX = 'C'; mM = d; mK = N; break;
  }
  char TRANSB = 0;
  switch (ModeB) {
  case NORMAL: TRANSB = 'N'; mN = B.GetN(); break;
  case TRANSPOSE: TRANSB = 'T'; mN = B.GetM(); break;
  case CONJ_TRANS: TRANSB = 'C'; mN = B.GetM(); break;
  }
  Y.Init(mM, mN);

  double *mB = B.GetPointer();
  double *mY = Y.GetPointer();

  sp_d_dgemm(TRANSX, TRANSB, mM, mN, mK, nnz, start, idx, X, mB, mY);

}


//--------------------------------------------------------------------------
void SPointArray::
MatMat(const SPointArray &B, DMatrix &Y, MatrixMode ModeX, MatrixMode ModeB) const {

  if ( (ModeX == NORMAL && ModeB == NORMAL && d != B.GetN()) ||
       (ModeX != NORMAL && ModeB == NORMAL && N != B.GetN()) ||
       (ModeX == NORMAL && ModeB != NORMAL && d != B.GetD()) ||
       (ModeX != NORMAL && ModeB != NORMAL && N != B.GetD()) ) {
    printf("SPointArray::MatMat. Error: Matrix sizes are not compatible. Function call takes no effect.\n");
    return;
  }

  INTEGER mM = 0, mN = 0, mK = 0;
  char TRANSX = 0;
  switch (ModeX) {
  case NORMAL: TRANSX = 'N'; mM = N; mK = d; break;
  case TRANSPOSE: TRANSX = 'T'; mM = d; mK = N; break;
  case CONJ_TRANS: TRANSX = 'C'; mM = d; mK = N; break;
  }
  char TRANSB = 0;
  switch (ModeB) {
  case NORMAL: TRANSB = 'N'; mN = B.GetD(); break;
  case TRANSPOSE: TRANSB = 'T'; mN = B.GetN(); break;
  case CONJ_TRANS: TRANSB = 'C'; mN = B.GetN(); break;
  }
  Y.Init(mM, mN);

  double *mY = Y.GetPointer();
  double *mB = B.GetPointerX();
  INTEGER nnz_B = B.GetNNZ();
  INTEGER *start_B = B.GetPointerStart();
  INTEGER *idx_B = B.GetPointerIdx();

  sp_s_dgemm(TRANSX, TRANSB, mM, mN, mK, nnz, start, idx, X, nnz_B, start_B, idx_B, mB, mY);

}


//--------------------------------------------------------------------------
INTEGER SPointArray::
RandomBipartition(INTEGER istart, INTEGER n, INTEGER N0,
                  INTEGER *perm, DPoint &normal, double &offset) {

  if (istart < 0 || istart+n > N) {
    printf("DPointArray::RandomPartition. Error: Invalid range. Function call takes no effect.\n");
    return 0;
  }

  // If n < 2*N0, there is no way to yield two clusters of size at least N0
  if (n < 2*N0) {
    return 0;
  }

  // Generate a random direction
  normal.Init(d);
  normal.SetStandardNormal();
  normal.Normalize();

  // Call shared subroutine Bipartition
  INTEGER m1 = Bipartition(istart, n, perm, normal, offset);
  return m1;

}


//--------------------------------------------------------------------------
INTEGER SPointArray::
PCABipartition(INTEGER istart, INTEGER n, INTEGER N0,
               INTEGER *perm, DPoint &normal, double &offset) {

  if (istart < 0 || istart+n > N) {
    printf("SPointArray::PCABipartition. Error: Invalid range. Function call takes no effect.\n");
    return 0;
  }

  // If n < 2*N0, there is no way to yield two clusters of size at least N0
  if (n < 2*N0) {
    return 0;
  }

  // Compute the principal direction. Considering that a sparse matrix
  // may have a large d, perhaps a more efficient algorithm to compute
  // the principal direction is the Lanczos method or the power
  // method, in some cases. However, due to robustness and coding
  // effort, we do not use these iterative methods here.
  //
  // Make a copy of the points
  DPointArray Y;
  GetSubset(istart, n, Y);
  // Compute the center
  DPoint c;
  Y.Center(c);
  // Center the points
  double *mc = c.GetPointer();
  double *mY = Y.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER j = 0; j < d; j++) {
    for (INTEGER i = 0; i < n; i++) {
      mY[i+j*n] -= mc[j];
    }
  }
  // Compute the dominant right singular vector of Y
  DMatrix YY, V;
  DVector lambda, v;
  normal.Init(d);
  if (n > d) {
    Y.MatMat(Y, YY, TRANSPOSE, NORMAL);
    YY.SymEig(lambda, V);
    V.GetColumn(d-1, v);
    normal.SetPoint(v.GetPointer(), d);
  }
  else {
    Y.MatMat(Y, YY, NORMAL, TRANSPOSE);
    YY.SymEig(lambda, V);
    V.GetColumn(n-1, v);
    DVector u;
    Y.MatVec(v, u, TRANSPOSE);
    double *mlambda = lambda.GetPointer();
    double sigma = sqrt(mlambda[n-1]);
    u.Divide(sigma);
    normal.SetPoint(u.GetPointer(), d);
  }

  // Call shared subroutine Bipartition
  INTEGER m1 = Bipartition(istart, n, perm, normal, offset);
  return m1;

}


//--------------------------------------------------------------------------
INTEGER SPointArray::
Bipartition(INTEGER istart, INTEGER n,
            INTEGER *perm, DPoint &normal, double &offset) {

  // Compute signed distances
  DPointArray YY;
  GetSubset(istart, n, YY);
  DVector sdist;
  YY.MatVec(normal, sdist, NORMAL);

  // Compute offset (median)
  DVector v = sdist;
  v.Sort(ASCEND);
  double *mv = v.GetPointer();
  INTEGER m1 = n/2;
  offset = (mv[m1-1]+mv[m1])/2.0;

  // Permute the points according to the signed distances.
  // p points to the start of the permutation array
  INTEGER *p = NULL;
  if (perm != NULL) {
    p = perm + istart;
  }
  // p2 is an auxilary permutation array with length only n
  INTEGER *p2 = NULL;
  New_1D_Array<INTEGER, INTEGER>(&p2, n);
  for (INTEGER i = 0; i < n; i++) {
    p2[i] = i + istart;
  }
  // pivot is the location that needs to be swapped later
  INTEGER pivot = 0;
  double *msdist = sdist.GetPointer();
  while (pivot < n && msdist[pivot] < offset) {
    pivot++;
  }
  // start j from pivot+1, swap pivot and j if needed
  for (INTEGER j = pivot+1; j < n; j++) {
    if (msdist[j] < offset) {
      // swap
      Swap<double>(msdist[pivot], msdist[j]);
      if (perm != NULL) {
        Swap<INTEGER>(p[pivot], p[j]);
      }
      Swap<INTEGER>(p2[pivot], p2[j]);
      // update pivot
      pivot++;
    }
  }
  // Permute the point array. Step 1: make temporary space
  INTEGER *start2 = NULL;
  INTEGER *idx2 = NULL;
  double *X2 = NULL;
  INTEGER nnz2 = start[istart+n]-start[istart];
  New_1D_Array<INTEGER, INTEGER>(&start2, n);
  New_1D_Array<INTEGER, INTEGER>(&idx2, nnz2);
  New_1D_Array<double, INTEGER>(&X2, nnz2);
  // Permute the point array. Step 2: permute on the temporary space
  INTEGER *idx2_ptr = idx2;
  double *X2_ptr = X2;
  start2[0] = start[istart];
  for (INTEGER i = 0; i < n; i++) {
    INTEGER j = p2[i];
    INTEGER nnz3 = start[j+1]-start[j];
    if (i+1 < n) {
      start2[i+1] = start2[i] + nnz3;
    }
    memcpy(idx2_ptr, idx+start[j], nnz3*sizeof(INTEGER));
    memcpy(X2_ptr, X+start[j], nnz3*sizeof(double));
    idx2_ptr += nnz3;
    X2_ptr += nnz3;
  }
  // Permute the point array. Step 3: copy back
  memcpy(start+istart, start2, n*sizeof(INTEGER));
  memcpy(idx+start[istart], idx2, nnz2*sizeof(INTEGER));
  memcpy(X+start[istart], X2, nnz2*sizeof(double));

  // Clean up
  Delete_1D_Array<INTEGER>(&p2);
  Delete_1D_Array<INTEGER>(&start2);
  Delete_1D_Array<INTEGER>(&idx2);
  Delete_1D_Array<double>(&X2);

  return m1;

}


//--------------------------------------------------------------------------
void SPointArray::
ComputeBBox(double **bbox, INTEGER &dim) {

  printf("SPointArray::BBoxBipartition. Error: This function is not implemented for sparse points, because the partitioning is potentially highly imbalanced, due to the excessive number of zeros. Function call takes no effect.\n");
  return;

}


//--------------------------------------------------------------------------
INTEGER SPointArray::
BBoxBipartition(INTEGER start, INTEGER n, INTEGER N0,
                INTEGER *perm, DPoint &normal, double &offset,
                const double *bbox, INTEGER &which_dim) {

  printf("SPointArray::BBoxBipartition. Error: This function is not implemented for sparse points, because the partitioning is potentially highly imbalanced, due to the excessive number of zeros. Function call takes no effect.\n");
  return -1;

}
