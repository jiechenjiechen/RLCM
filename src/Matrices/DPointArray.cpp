#include "DPointArray.hpp"

#define INITVAL_N 0
#define INITVAL_d 0


//--------------------------------------------------------------------------
DPointArray::
DPointArray() {
  X = NULL;
  N = INITVAL_N;
  d = INITVAL_d;
  Init();
}


//--------------------------------------------------------------------------
void DPointArray::
Init(void) {
  Init(0,0);
}


//--------------------------------------------------------------------------
void DPointArray::
Init(INTEGER N_, INTEGER d_) {
  if (N*d != N_*d_) {
    ReleaseAllMemory();
  }
  if (!X) {
    New_1D_Array<double, INTEGER>(&X, N_*d_);
  }
  else {
    memset(X, 0, N_*d_*sizeof(double));
  }
  N = N_;
  d = d_;
}


//--------------------------------------------------------------------------
void DPointArray::
ReleaseAllMemory(void) {
  N = INITVAL_N;
  d = INITVAL_d;
  Delete_1D_Array<double>(&X);
}


//--------------------------------------------------------------------------
DPointArray::
DPointArray(INTEGER N_, INTEGER d_) {
  X = NULL;
  N = INITVAL_N;
  d = INITVAL_d;
  Init(N_, d_);
}


//--------------------------------------------------------------------------
DPointArray::
DPointArray(const DPointArray& G) {
  X = NULL;
  N = INITVAL_N;
  d = INITVAL_d;
  Init();
  DeepCopy(G);
}


//--------------------------------------------------------------------------
DPointArray& DPointArray::
operator= (const DPointArray& G) {
  if (this != &G) {
    DeepCopy(G);
  }
  return *this;
}


//--------------------------------------------------------------------------
void DPointArray::
DeepCopy(const DPointArray &G) {
  if (N*d != G.N*G.d) {
    ReleaseAllMemory();
  }
  if (G.X) {
    if (!X) {
      New_1D_Array<double, INTEGER>(&X, G.N*G.d);
    }
    memcpy(X, G.X, G.N*G.d*sizeof(double));
  }
  N = G.N;
  d = G.d;
}


//--------------------------------------------------------------------------
DPointArray::
~DPointArray() {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
INTEGER DPointArray::
GetD(void) const {
  return d;
}


//--------------------------------------------------------------------------
INTEGER DPointArray::
GetN(void) const {
  return N;
}


//--------------------------------------------------------------------------
void DPointArray::
GetPoint(INTEGER i, DPoint &x) const {
  if (i < 0 || i >= N) {
    printf("DPointArray::GetPoint. Error: Invalid index. Function call takes no effect.\n");
    return;
  }
  x.Init(d);
  double *mx = x.GetPointer();
  for (INTEGER j = 0; j < d; j++) {
    mx[j] = X[i+j*N];
  }
}


//--------------------------------------------------------------------------
void DPointArray::
GetSubset(INTEGER start, INTEGER n, DPointArray &Y) const {
  if (start < 0 || start+n > N) {
    printf("DPointArray::GetSubset. Error: Invalid range. Function call takes no effect.\n");
    return;
  }
  Y.Init(n,d);
  double *mY = Y.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER j = 0; j < d; j++) {
    memcpy(mY+j*n, X+j*N+start, n*sizeof(double));
  }
}


//--------------------------------------------------------------------------
void DPointArray::
GetSubset(INTEGER *idx, INTEGER n, DPointArray &Y) const {
  Y.Init(n,d);
  double *mY = Y.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER j = 0; j < d; j++) {
    for (INTEGER i = 0; i < n; i++) {
      mY[i+j*n] = X[idx[i]+j*N];
    }
  }
}


//--------------------------------------------------------------------------
double* DPointArray::
GetPointer(void) const {
  return X;
}


//--------------------------------------------------------------------------
void DPointArray::
AsDMatrix(DMatrix &A) const {
  A.Init(N,d);
  double *mA = A.GetPointer();
  memcpy(mA, X, N*d*sizeof(double));
}


//--------------------------------------------------------------------------
void DPointArray::
SetUniformRandom01(void) {
  UniformRandom01(X, N*d);
}


//--------------------------------------------------------------------------
void DPointArray::
SetUniformSphere(void) {
  StandardNormal(X, N*d);
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER j = 0; j < N; j++) {
    double r = 0.0;
    for (INTEGER i = 0; i < d; i++) {
      r += Square(X[j+i*N]);
    }
    r = sqrt(r);
    for (INTEGER i = 0; i < d; i++) {
      X[j+i*N] /= r;
    }
  }
}


//--------------------------------------------------------------------------
void DPointArray::
SetRegularGrid(INTEGER d_, const INTEGER *Dim,
               const double *Lower, const double *Upper) {
  INTEGER N_ = 1;
  for (INTEGER i = 0; i < d_; i++) {
    N_ *= Dim[i];
  }
  Init(N_, d_);
  for (INTEGER i = 0; i < d_; i++) {
    INTEGER SizeThisDim = Dim[i];
    INTEGER SizePrevDim = 1;
    for (INTEGER j = 0; j < i; j++) {
      SizePrevDim *= Dim[j];
    }
    INTEGER SizePostDim = 1;
    for (INTEGER j = i+1; j < d_; j++) {
      SizePostDim *= Dim[j];
    }
    double Val = Lower[i];
    double Delta = (Upper[i]-Lower[i])/(Dim[i]-1); // May be inf, but ok
    for (INTEGER j = 0; j < SizeThisDim; j++) {
      for (INTEGER p = 0; p < SizePostDim; p++) {
        for (INTEGER q = 0; q < SizePrevDim; q++) {
          X[q+(j+p*SizeThisDim)*SizePrevDim+i*N] = Val;
        }
      }
      Val += Delta;
    }
  }
}


//--------------------------------------------------------------------------
void DPointArray::
PrintPointArray(const char *name) const {
  printf("%s:\n", name);
  for (INTEGER i = 0; i < N; i++) {
    for (INTEGER j = 0; j < d; j++) {
      printf("%g ", X[i+j*N]);
    }
    printf("\n");
  }
  printf("\n");
}


//--------------------------------------------------------------------------
void DPointArray::
Center(DPoint &c) const {
  c.Init(d);
  double *mc = c.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER j = 0; j < d; j++) {
    for (INTEGER i = 0; i < N; i++) {
      mc[j] += X[i+N*j];
    }
    mc[j] /= N;
  }
}


//--------------------------------------------------------------------------
double DPointArray::
StdvDist(void) const {
  DPoint c;
  Center(c);
  double dist2 = 0.0;
#ifdef USE_OPENMP
#pragma omp parallel for reduction(+:dist2)
#endif
  for (INTEGER i = 0; i < N; i++) {
    DPoint x;
    GetPoint(i, x);
    dist2 += x.Dist2(c);
  }
  double ret = sqrt(dist2/N);
  return ret;
}


//--------------------------------------------------------------------------
double DPointArray::
MaxPointNorm2(void) const {
  DVector PointNorm2(N);
  double *mPointNorm2 = PointNorm2.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < N; i++) {
    DPoint x;
    GetPoint(i, x);
    mPointNorm2[i] = x.InProd(x);
  }
  double ret = PointNorm2.Max();
  return ret;
}


//--------------------------------------------------------------------------
void DPointArray::
MatVec(const DVector &b, DVector &y, MatrixMode ModeX) const {

  if ( (ModeX == NORMAL && d != b.GetN()) ||
       (ModeX != NORMAL && N != b.GetN()) ) {
    printf("DPointArray::MatVec. Error: Matrix/vector sizes are not compatible. Function call takes no effect.\n");
    return;
  }

  char TRANS;
  switch (ModeX) {
  case NORMAL: TRANS = 'N'; y.Init(N); break;
  case TRANSPOSE: TRANS = 'T'; y.Init(d); break;
  case CONJ_TRANS: TRANS = 'C'; y.Init(d); break;
  }

  double *mb = b.GetPointer();
  double *my = y.GetPointer();

  dgemv_(&TRANS, &N, &d, &ONE, X, &N, mb, &ONEi, &ZERO, my, &ONEi);

}


//--------------------------------------------------------------------------
void DPointArray::
MatVec(const DPoint &b, DVector &y, MatrixMode ModeX) const {

  if ( (ModeX == NORMAL && d != b.GetD()) ||
       (ModeX != NORMAL && N != b.GetD()) ) {
    printf("DPointArray::MatVec. Error: Matrix/vector sizes are not compatible. Function call takes no effect.\n");
    return;
  }

  char TRANS;
  switch (ModeX) {
  case NORMAL: TRANS = 'N'; y.Init(N); break;
  case TRANSPOSE: TRANS = 'T'; y.Init(d); break;
  case CONJ_TRANS: TRANS = 'C'; y.Init(d); break;
  }

  double *mb = b.GetPointer();
  double *my = y.GetPointer();

  dgemv_(&TRANS, &N, &d, &ONE, X, &N, mb, &ONEi, &ZERO, my, &ONEi);

}


//--------------------------------------------------------------------------
void DPointArray::
MatMat(const DMatrix &B, DMatrix &Y, MatrixMode ModeX, MatrixMode ModeB) const {

  if ( (ModeX == NORMAL && ModeB == NORMAL && d != B.GetM()) ||
       (ModeX != NORMAL && ModeB == NORMAL && N != B.GetM()) ||
       (ModeX == NORMAL && ModeB != NORMAL && d != B.GetN()) ||
       (ModeX != NORMAL && ModeB != NORMAL && N != B.GetN()) ) {
    printf("DPointArray::MatMat. Error: Matrix sizes are not compatible. Function call takes no effect.\n");
    return;
  }

  INTEGER mM, mN, mK;
  char TRANSX;
  switch (ModeX) {
  case NORMAL: TRANSX = 'N'; mM = N; mK = d; break;
  case TRANSPOSE: TRANSX = 'T'; mM = d; mK = N; break;
  case CONJ_TRANS: TRANSX = 'C'; mM = d; mK = N; break;
  }
  char TRANSB;
  switch (ModeB) {
  case NORMAL: TRANSB = 'N'; mN = B.GetN(); break;
  case TRANSPOSE: TRANSB = 'T'; mN = B.GetM(); break;
  case CONJ_TRANS: TRANSB = 'C'; mN = B.GetM(); break;
  }
  Y.Init(mM, mN);

  INTEGER LDX = N, LDB = B.GetM();
  double *mB = B.GetPointer();
  double *mY = Y.GetPointer();

  dgemm_(&TRANSX, &TRANSB, &mM, &mN, &mK, &ONE, X, &LDX, mB, &LDB, &ZERO, mY, &mM);

}


//--------------------------------------------------------------------------
void DPointArray::
MatMat(const DPointArray &B, DMatrix &Y, MatrixMode ModeX, MatrixMode ModeB) const {

  if ( (ModeX == NORMAL && ModeB == NORMAL && d != B.GetN()) ||
       (ModeX != NORMAL && ModeB == NORMAL && N != B.GetN()) ||
       (ModeX == NORMAL && ModeB != NORMAL && d != B.GetD()) ||
       (ModeX != NORMAL && ModeB != NORMAL && N != B.GetD()) ) {
    printf("DPointArray::MatMat. Error: Matrix sizes are not compatible. Function call takes no effect.\n");
    return;
  }

  INTEGER mM = 0, mN = 0, mK = 0;
  char TRANSX;
  switch (ModeX) {
  case NORMAL: TRANSX = 'N'; mM = N; mK = d; break;
  case TRANSPOSE: TRANSX = 'T'; mM = d; mK = N; break;
  case CONJ_TRANS: TRANSX = 'C'; mM = d; mK = N; break;
  }
  char TRANSB;
  switch (ModeB) {
  case NORMAL: TRANSB = 'N'; mN = B.GetD(); break;
  case TRANSPOSE: TRANSB = 'T'; mN = B.GetN(); break;
  case CONJ_TRANS: TRANSB = 'C'; mN = B.GetN(); break;
  }
  Y.Init(mM, mN);

  INTEGER LDX = N, LDB = B.GetN();
  double *mB = B.GetPointer();
  double *mY = Y.GetPointer();

  dgemm_(&TRANSX, &TRANSB, &mM, &mN, &mK, &ONE, X, &LDX, mB, &LDB, &ZERO, mY, &mM);

}


//--------------------------------------------------------------------------
INTEGER DPointArray::
RandomBipartition(INTEGER start, INTEGER n, INTEGER N0,
                  INTEGER *perm, DPoint &normal, double &offset) {

  if (start < 0 || start + n > N) {
    printf("DPointArray::RandomBipartition. Error: Invalid range. Function call takes no effect.\n");
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
  INTEGER m1 = Bipartition(start, n, perm, normal, offset);
  return m1;

}


//--------------------------------------------------------------------------
INTEGER DPointArray::
PCABipartition(INTEGER start, INTEGER n, INTEGER N0,
               INTEGER *perm, DPoint &normal, double &offset) {

  if (start < 0 || start + n > N) {
    printf("DPointArray::PCABipartition. Error: Invalid range. Function call takes no effect.\n");
    return 0;
  }

  // If n < 2*N0, there is no way to yield two clusters of size at least N0
  if (n < 2*N0) {
    return 0;
  }

  // Compute the principal direction.
  //
  // Make a copy of the points
  DPointArray Y;
  GetSubset(start, n, Y);
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
  INTEGER m1 = Bipartition(start, n, perm, normal, offset);
  return m1;

}


//--------------------------------------------------------------------------
INTEGER DPointArray::
Bipartition(INTEGER start, INTEGER n,
            INTEGER *perm, const DPoint &normal, double &offset) {

  // Compute signed distances
  DPointArray YY;
  GetSubset(start, n, YY);
  DVector sdist;
  YY.MatVec(normal, sdist, NORMAL);

  // Compute offset (median)
  DVector v = sdist;
  v.Sort(ASCEND);
  double *mv = v.GetPointer();
  INTEGER m1 = n/2;
  offset = (mv[m1-1]+mv[m1])/2.0;

  // Permute the points according to the signed distances.
  // Y points to the start of the points
  double *Y = X + start;
  // p points to the start of the permutation array
  INTEGER *p = NULL;
  if (perm != NULL) {
    p = perm + start;
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
      for (INTEGER i = 0; i < d; i++) {
        Swap<double>(Y[pivot+N*i], Y[j+N*i]);
      }
      Swap<double>(msdist[pivot], msdist[j]);
      if (perm != NULL) {
        Swap<INTEGER>(p[pivot], p[j]);
      }
      // update pivot
      pivot++;
    }
  }

  return m1;

}


//--------------------------------------------------------------------------
void DPointArray::
ComputeBBox(double **bbox, INTEGER &dim) {

  dim = d;
  New_1D_Array<double, INTEGER>(bbox, d*2);

  double *X_ptr = X;
  for (INTEGER j = 0; j < d; j++) {
    (*bbox)[j] = DBL_MAX; // to store min
    (*bbox)[j+d] = -DBL_MAX; // to store max
    for (INTEGER i = 0; i < N; i++) {
      if ((*bbox)[j] > *X_ptr) {
        (*bbox)[j] = *X_ptr;
      }
      if ((*bbox)[j+d] < *X_ptr) {
        (*bbox)[j+d] = *X_ptr;
      }
      X_ptr++;
    }
  }

}


//--------------------------------------------------------------------------
INTEGER DPointArray::
BBoxBipartition(INTEGER start, INTEGER n, INTEGER N0,
                INTEGER *perm, DPoint &normal, double &offset,
                const double *bbox, INTEGER &which_dim) {

  if (start < 0 || start + n > N) {
    printf("DPointArray::BBoxBipartition. Error: Invalid range. Function call takes no effect.\n");
    return 0;
  }

  // If n < 2*N0, there is no way to yield two clusters of size at least N0
  if (n < 2*N0) {
    return 0;
  }

  // Find the longest dimension
  which_dim = -1;
  double len = -DBL_MAX;
  for (INTEGER i = 0; i < d; i++) {
    double cur_len = bbox[i+d] - bbox[i];
    if (len < cur_len) {
      which_dim = i;
      len = cur_len;
    }
  }

  // normal
  normal.Init(d);
  double *ptr = normal.GetPointer();
  ptr[which_dim] = 1;

  // offset
  offset = ( bbox[which_dim+d] + bbox[which_dim] ) / 2.0;

  // signed distances
  DVector sdist(n);
  double *msdist = sdist.GetPointer();
  memcpy(msdist, X+which_dim*N+start, n*sizeof(double));

  // Permute the points according to the signed distances.
  // Y points to the start of the points
  double *Y = X + start;
  // p points to the start of the permutation array
  INTEGER *p = NULL;
  if (perm != NULL) {
    p = perm + start;
  }
  // pivot is the location that needs to be swapped later
  INTEGER pivot = 0;
  while (pivot < n && msdist[pivot] < offset) {
    pivot++;
  }
  // start j from pivot+1, swap pivot and j if needed
  for (INTEGER j = pivot+1; j < n; j++) {
    if (msdist[j] < offset) {
      // swap
      for (INTEGER i = 0; i < d; i++) {
        Swap<double>(Y[pivot+N*i], Y[j+N*i]);
      }
      Swap<double>(msdist[pivot], msdist[j]);
      if (perm != NULL) {
        Swap<INTEGER>(p[pivot], p[j]);
      }
      // update pivot
      pivot++;
    }
  }

  // m1
  INTEGER m1 = pivot;
  return m1;

}
