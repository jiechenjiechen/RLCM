#include "spblas.hpp"


//--------------------------------------------------------------------------
// y = x, x sparse, y dense
void sp_n_ds2d(INTEGER N, INTEGER NNZX, INTEGER *IDXX, double *X, double *Y) {
  for (INTEGER i = 0; i < NNZX; i++) {
    Y[IDXX[i]] = X[i];
  }
}


//--------------------------------------------------------------------------
// z = dot(x,y), x sparse, y dense
double sp_d_ddot(INTEGER N, INTEGER NNZX, INTEGER *IDXX, double *X, double *Y) {
  double ret = 0.0;
  for (INTEGER i = 0; i < NNZX; i++) {
    ret += X[i] * Y[IDXX[i]];
  }
  return ret;
}


//--------------------------------------------------------------------------
// z = dot(x,y), x sparse, y sparse
double sp_s_ddot(INTEGER N, INTEGER NNZX, INTEGER *IDXX, double *X, INTEGER NNZY, INTEGER *IDXY, double *Y) {
  double ret = 0.0;
  INTEGER i = 0, j = 0;
  while (1) {
    if (i == NNZX || j == NNZY) {
      break;
    }
    INTEGER IDXXi = IDXX[i];
    INTEGER IDXYj = IDXY[j];
    if (IDXXi > IDXYj) {
      j++;
    }
    else if (IDXXi < IDXYj) {
      i++;
    }
    else /* IDXXi == IDXYj */ {
      ret += X[i] * Y[j];
      i++;
      j++;
    }
  }
  return ret;
}


//--------------------------------------------------------------------------
// z = func(x-y), x sparse, y dense, func is an elementwise function that eventually sums up
double sp_d_ddist(INTEGER N, INTEGER NNZX, INTEGER *IDXX, double *X, double *Y, double(*FUNC)(double)) {
  double ret = 0.0;
  INTEGER i = 0, j = 0, jstart = 0, jend = 0;
  for (i = 0; i < NNZX; i++) {
    if (i == 0) {
      jstart = 0;
      jend = IDXX[0];
    }
    else {
      jstart = IDXX[i-1]+1;
      jend = IDXX[i];
    }
    for (j = jstart; j < jend; j++) {
      ret += (*FUNC)( -Y[j] );
    }
    ret += (*FUNC)( X[i]-Y[IDXX[i]] );
  }
  jstart = IDXX[NNZX-1]+1;
  jend = N;
  for (j = jstart; j < jend; j++) {
    ret += (*FUNC)( -Y[j] );
  }
  return ret;
}


//--------------------------------------------------------------------------
// z = func(x-y), x sparse, y sparse, func is an elementwise function that eventually sums up
double sp_s_ddist(INTEGER N, INTEGER NNZX, INTEGER *IDXX, double *X, INTEGER NNZY, INTEGER *IDXY, double *Y, double(*FUNC)(double)) {
  double ret = 0.0;
  INTEGER i = 0, j = 0;
  while (1) {
    if (i == NNZX || j == NNZY) {
      if (i != NNZX && j == NNZY) {
        while (i < NNZX) {
          ret += (*FUNC)( X[i] );
          i++;
        }
        break;
      }
      else if (i == NNZX && j != NNZY) {
        while (j < NNZY) {
          ret += (*FUNC)( -Y[j] );
          j++;
        }
        break;
      }
      else /* i == NNZX && j == NNZY */ {
        break;
      }
    }
    INTEGER IDXXi = IDXX[i];
    INTEGER IDXYj = IDXY[j];
    if (IDXXi > IDXYj) {
      ret += (*FUNC)( -Y[j] );
      j++;
    }
    else if (IDXXi < IDXYj) {
      ret += (*FUNC)( X[i] );
      i++;
    }
    else /* IDXXi == IDXYj */ {
      ret += (*FUNC)( X[i] - Y[j] );
      i++;
      j++;
    }
  }
  return ret;
}


//--------------------------------------------------------------------------
// z = func((x-y)/s), x sparse, y dense, s dense, func is an elementwise function that eventually sums up
double sp_d_ddists(INTEGER N, INTEGER NNZX, INTEGER *IDXX, double *X, double *Y, double *S, double(*FUNC)(double)) {
  double ret = 0.0;
  INTEGER i = 0, j = 0, jstart = 0, jend = 0;
  for (i = 0; i < NNZX; i++) {
    if (i == 0) {
      jstart = 0;
      jend = IDXX[0];
    }
    else {
      jstart = IDXX[i-1]+1;
      jend = IDXX[i];
    }
    for (j = jstart; j < jend; j++) {
      ret += (*FUNC)( -Y[j] / S[j] );
    }
    ret += (*FUNC)( (X[i] - Y[IDXX[i]]) / S[IDXX[i]] );
  }
  jstart = IDXX[NNZX-1]+1;
  jend = N;
  for (j = jstart; j < jend; j++) {
    ret += (*FUNC)( -Y[j] / S[j] );
  }
  return ret;
}


//--------------------------------------------------------------------------
// z = func((x-y)/s), x sparse, y sparse, s dense, func is an elementwise function that eventually sums up
double sp_s_ddists(INTEGER N, INTEGER NNZX, INTEGER *IDXX, double *X, INTEGER NNZY, INTEGER *IDXY, double *Y, double *S, double(*FUNC)(double)) {
  double ret = 0.0;
  INTEGER i = 0, j = 0;
  while (1) {
    if (i == NNZX || j == NNZY) {
      if (i != NNZX && j == NNZY) {
        while (i < NNZX) {
          ret += (*FUNC)( X[i] / S[IDXX[i]] );
          i++;
        }
        break;
      }
      else if (i == NNZX && j != NNZY) {
        while (j < NNZY) {
          ret += (*FUNC)( -Y[j] / S[IDXY[j]] );
          j++;
        }
        break;
      }
      else /* i == NNZX && j == NNZY */ {
        break;
      }
    }
    INTEGER IDXXi = IDXX[i];
    INTEGER IDXYj = IDXY[j];
    if (IDXXi > IDXYj) {
      ret += (*FUNC)( -Y[j] / S[IDXYj] );
      j++;
    }
    else if (IDXXi < IDXYj) {
      ret += (*FUNC)( X[i] / S[IDXXi] );
      i++;
    }
    else /* IDXXi == IDXYj */ {
      ret += (*FUNC)( (X[i] - Y[j]) / S[IDXXi] );
      i++;
      j++;
    }
  }
  return ret;
}


//--------------------------------------------------------------------------
// z = sum(2*x.*y./(x+y)), used only for x,y >= 0.
// If x[i]+y[i]=0, then 2*x[i]*y[i]/(x[i]+y[i])=0.
// x sparse, y dense
double sp_d_dchi2(INTEGER N, INTEGER NNZX, INTEGER *IDXX, double *X, double *Y) {
  double ret = 0.0;
  for (INTEGER i = 0; i < NNZX; i++) {
    double denom = X[i] + Y[IDXX[i]];
    if (denom != 0.0) {
      double numer = X[i] * Y[IDXX[i]];
      ret += (numer / denom);
    }
  }
  ret *= 2.0;
  return ret;
}


//--------------------------------------------------------------------------
// z = sum(2*x.*y./(x+y)), used only for x,y >= 0.
// If x[i]+y[i]=0, then 2*x[i]*y[i]/(x[i]+y[i])=0.
// x sparse, y sparse
double sp_s_dchi2(INTEGER N, INTEGER NNZX, INTEGER *IDXX, double *X, INTEGER NNZY, INTEGER *IDXY, double *Y) {
  double ret = 0.0;
  INTEGER i = 0, j = 0;
  while (1) {
    if (i == NNZX || j == NNZY) {
      break;
    }
    INTEGER IDXXi = IDXX[i];
    INTEGER IDXYj = IDXY[j];
    if (IDXXi > IDXYj) {
      j++;
    }
    else if (IDXXi < IDXYj) {
      i++;
    }
    else /* IDXXi == IDXYj */ {
      double denom = X[i] + Y[j];
      if (denom != 0.0) {
        double numer = X[i] * Y[j];
        ret += (numer / denom);
      }
      i++;
      j++;
    }
  }
  ret *= 2.0;
  return ret;
}


//--------------------------------------------------------------------------
// y = mode(A)*x, A sparse, x dense, y dense
void sp_d_dgemv(char TRANS, INTEGER M, INTEGER N, INTEGER NNZA, INTEGER *STARTA, INTEGER *IDXA, double *A, double *X, double *Y) {
  if (TRANS == 'N') { //----------------------------------------------------
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (INTEGER i = 0; i < M; i++) {
      double Yi = 0;
      for (INTEGER j = STARTA[i]; j < STARTA[i+1]; j++) {
        Yi += A[j] * X[IDXA[j]];
      }
      Y[i] = Yi;
    }
  }
  else { //-----------------------------------------------------------------
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (INTEGER i = 0; i < N; i++) {
      Y[i] = 0;
    }
    for (INTEGER i = 0; i < M; i++) {
      double Xi = X[i];
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
      for (INTEGER j = STARTA[i]; j < STARTA[i+1]; j++) {
        Y[IDXA[j]] += A[j] * Xi;
      }
    }
  } //----------------------------------------------------------------------
}


//--------------------------------------------------------------------------
// y = mode(A)*x, A sparse, x sparse, y dense
void sp_s_dgemv(char TRANS, INTEGER M, INTEGER N, INTEGER NNZA, INTEGER *STARTA, INTEGER *IDXA, double *A, INTEGER NNZX, INTEGER *IDXX, double *X, double *Y) {
  if (TRANS == 'N') { //----------------------------------------------------
    double *DENSE_X = (double *)calloc(N,sizeof(double));
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (INTEGER i = 0; i < NNZX; i++) {
      DENSE_X[IDXX[i]] = X[i];
    }
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (INTEGER i = 0; i < M; i++) {
      double Yi = 0;
      for (INTEGER j = STARTA[i]; j < STARTA[i+1]; j++) {
        Yi += A[j] * DENSE_X[IDXA[j]];
      }
      Y[i] = Yi;
    }
    free(DENSE_X);
  }
  else { //-----------------------------------------------------------------
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (INTEGER i = 0; i < N; i++) {
      Y[i] = 0;
    }
    for (INTEGER i = 0; i < NNZX; i++) {
      INTEGER jj = IDXX[i];
      double Xi = X[i];
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
      for (INTEGER j = STARTA[jj]; j < STARTA[jj+1]; j++) {
        Y[IDXA[j]] += A[j] * Xi;
      }
    }
  } //----------------------------------------------------------------------
}


//--------------------------------------------------------------------------
// C = mode(A)*mode(B), A sparse, B dense, C dense
void sp_d_dgemm(char TRANSA, char TRANSB, INTEGER M, INTEGER N, INTEGER K, INTEGER NNZA, INTEGER *STARTA, INTEGER *IDXA, double *A, double *B, double *C) {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < N; i++) {
    memset(C+M*i, 0, M*sizeof(double));
  }
  if (TRANSA == 'N' && TRANSB == 'N') { //----------------------------
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (INTEGER l = 0; l < N; l++) {
      for (INTEGER i = 0; i < M; i++) {
        for (INTEGER j = STARTA[i]; j < STARTA[i+1]; j++) {
          INTEGER jj = IDXA[j];
          C[i+M*l] += A[j] * B[jj+K*l];
        }
      }
    }
  }
  else if (TRANSA == 'N' && TRANSB != 'N') { //-----------------------
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (INTEGER l = 0; l < N; l++) {
      for (INTEGER i = 0; i < M; i++) {
        for (INTEGER j = STARTA[i]; j < STARTA[i+1]; j++) {
          INTEGER jj = IDXA[j];
          C[i+M*l] += A[j] * B[l+N*jj];
        }
      }
    }
  }
  else if (TRANSA != 'N' && TRANSB == 'N') { //-----------------------
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (INTEGER l = 0; l < N; l++) {
      for (INTEGER i = 0; i < K; i++) {
        for (INTEGER j = STARTA[i]; j < STARTA[i+1]; j++) {
          INTEGER jj = IDXA[j];
          C[jj+M*l] += A[j] * B[i+K*l];
        }
      }
    }
  }
  else /* (TRANSA != 'N' && TRANSB != 'N') */ { //--------------------
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (INTEGER l = 0; l < N; l++) {
      for (INTEGER i = 0; i < K; i++) {
        for (INTEGER j = STARTA[i]; j < STARTA[i+1]; j++) {
          INTEGER jj = IDXA[j];
          C[jj+M*l] += A[j] * B[l+N*i];
        }
      }
    }
  } //----------------------------------------------------------------------
}


//--------------------------------------------------------------------------
// C = mode(A)*mode(B), A sparse, B sparse, C dense
void sp_s_dgemm(char TRANSA, char TRANSB, INTEGER M, INTEGER N, INTEGER K, INTEGER NNZA, INTEGER *STARTA, INTEGER *IDXA, double *A, INTEGER NNZB, INTEGER *STARTB, INTEGER *IDXB, double *B, double *C) {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < N; i++) {
    memset(C+M*i, 0, M*sizeof(double));
  }
  if (TRANSA == 'N' && TRANSB == 'N') { //----------------------------
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (INTEGER i = 0; i < M; i++) {
      for (INTEGER j = STARTA[i]; j < STARTA[i+1]; j++) {
        INTEGER jj = IDXA[j];
        double Aj = A[j];
        for (INTEGER l = STARTB[jj]; l < STARTB[jj+1]; l++) {
          INTEGER ll = IDXB[l];
          C[i+M*ll] += Aj * B[l];
        }
      }
    }
  }
  else if (TRANSA == 'N' && TRANSB != 'N') { //-----------------------
    double *DENSE_B_ROW = (double *)calloc(K,sizeof(double));
    for (INTEGER i = 0; i < N; i++) {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
      for (INTEGER l = STARTB[i]; l < STARTB[i+1]; l++) {
        DENSE_B_ROW[IDXB[l]] = B[l];
      }
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
      for (INTEGER j = 0; j < M; j++) {
        INTEGER ji = j+M*i;
        for (INTEGER l = STARTA[j]; l < STARTA[j+1]; l++) {
          INTEGER ll = IDXA[l];
          C[ji] += A[l] * DENSE_B_ROW[ll];
        }
      }
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
      for (INTEGER l = STARTB[i]; l < STARTB[i+1]; l++) {
        DENSE_B_ROW[IDXB[l]] = 0;
      }
    }
    free(DENSE_B_ROW);
  }
  else if (TRANSA != 'N' && TRANSB == 'N') { //-----------------------
    for (INTEGER i = 0; i < K; i++) {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
      for (INTEGER j = STARTA[i]; j < STARTA[i+1]; j++) {
        INTEGER jj = IDXA[j];
        for (INTEGER l = STARTB[i]; l < STARTB[i+1]; l++) {
          INTEGER ll = IDXB[l];
          C[jj+M*ll] += A[j] * B[l];
        }
      }
    }
  }
  else /* (TRANSA != 'N' && TRANSB != 'N') */ { //--------------------
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (INTEGER i = 0; i < N; i++) {
      for (INTEGER j = STARTB[i]; j < STARTB[i+1]; j++) {
        INTEGER jj = IDXB[j];
        double Bj = B[j];
        for (INTEGER l = STARTA[jj]; l < STARTA[jj+1]; l++) {
          INTEGER ll = IDXA[l];
          C[ll+M*i] += A[l] * Bj;
        }
      }
    }
  } //----------------------------------------------------------------------
}
