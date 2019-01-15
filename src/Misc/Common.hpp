// This file contains miscellaneous stuffs supporting the
// implementation of other classes.
//
// Some functions are parallelized whereas others are not.

#ifndef _COMMON_
#define _COMMON_

#include <math.h>
#include <time.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <limits.h>
#include <typeinfo>
#ifdef USE_OPENBLAS
#include <openblas/cblas.h>
#elif defined USE_OPENMP
#include <omp.h>
#endif
#ifdef USE_MKL
#include <mkl.h>
#endif

#ifdef USE_LONG
typedef long INTEGER;
#define String2Integer atol
#else
typedef int INTEGER;
#define String2Integer atoi
#endif

enum MatrixMode { NORMAL, TRANSPOSE, CONJ_TRANS }; // Used by MatVec and MatMat
enum MatrixType { GENERAL, SPD }; // Matrix type
enum SelectType { LHP, RHP }; // Left/Right-half plane; used by only RealSchur
enum SortType { ASCEND, DESCEND }; // Sort order
enum PartMethod { RAND, PCA, BBOX }; // Method for partitioning a point set
enum MatState { UNFACT, LU_FACT, CHOL_FACT }; // Sometimes factored form is used
enum TriUPLO { UPPER, LOWER }; // Upper and lower triangular

#define IF_NUMERIC_TYPE                                             \
    if (typeid(T) == typeid(int) || typeid(T) == typeid(long) ||    \
        typeid(T) == typeid(float) || typeid(T) == typeid(double))

typedef INTEGER LOGICAL;
#define FTRUE  1
#define FFALSE 0

#define EPS     2.220446049250313e-16
#define EPS_2   1.110223024625157e-16
#define SqrtEPS 1.490116119384766e-08
#define EPSx10  2.220446049250313e-15
#define PIx2    6.283185307179586

struct LogDet {
  double LogAbsDet;
  INTEGER Sign;
};

#define PREPARE_CLOCK(is_timing)                 \
    struct timespec time_start, time_end;        \
    double ELAPSED_TIME = 0.0;                   \
    bool timing = is_timing;

#define START_CLOCK                                 \
    if (timing) {                                   \
      clock_gettime(CLOCK_MONOTONIC, &time_start);  \
    }

#define END_CLOCK                                                       \
    if (timing) {                                                       \
      clock_gettime(CLOCK_MONOTONIC, &time_end);                        \
      ELAPSED_TIME = time_end.tv_sec - time_start.tv_sec;               \
      ELAPSED_TIME += (time_end.tv_nsec - time_start.tv_nsec) / 1e9;    \
    }

const INTEGER ONEi = 1;
const double ZERO = 0.0;
const double ONE = 1.0;
const double MINUS_ONE = -1.0;
const char TRANS_N = 'N';
const char TRANS_T = 'T';
const char UPLO_U = 'U';
const char UPLO_L = 'L';
const char DIAG_N = 'N';
const char FACT_N = 'N';
const char EQUED_N = 'N';
const char NORM_1 = '1';
const char JOBZ_N = 'N';
const char JOBZ_V = 'V';
const char JOBVL_N = 'N';
const char JOBVR_N = 'N';
const char JOBVR_V = 'V';
const char JOBVS_V = 'V';
const char SORT_S = 'S';
const char JOBU_N = 'N';
const char JOBVT_N = 'N';
const char SIDE_L = 'L';


//--------------------------------------------------------------------------
// x^2
double Square(double x);

// L1 difference between two arrays. Paralellized
double Diff1(const double *a, const double *b, INTEGER n);
double Diff1(const INTEGER *a, const INTEGER *b, INTEGER n);

// Endianness. Not Parallelized
bool IsBigEndian(void);
void Swap4Bytes(char *desc, char *src);
void Swap8Bytes(char *desc, char *src);
void SwapNBytes(char *desc, char *src, INTEGER N);

// Compute the parity of a permutation of order n. The permutation is
// the IPIV array output from DGETRF. It is different from the Perm
// array notation used in DVector.hpp. Specifically, Perm can be
// obtained from IPIV through the following operations:
//
// 1. Initialize Perm = 0:n-1
// 2. For i = 0:n-2, swap Perm[i] with Perm[IPIV[i]]
//
// This function is parallelized.
//
// Return value 0: even; 1: odd
INTEGER Parity(INTEGER *ipiv, INTEGER n);

// All functions related to random numbers are not parallelized.
//
// Generate [0,1] uniform random numbers
void UniformRandom01(double *a, INTEGER n);

// Generate Gaussian random numbers
void StandardNormal(double *a, INTEGER n);

// Generate student-t of degree 1 random numbers
void StudentT1(double *a, INTEGER n);

// Generate a random vector of multivariate student-t of degree 1
void MultivariateStudentT1(double *a, INTEGER n);

// Generate sech random numbers
void RandomSech(double *a, INTEGER n);

// Generate k random integers in [0,n) without repetition
void RandPerm(INTEGER n, INTEGER k, INTEGER *a);

// Needed by qsort and bsearch
int CompareNaturalOrderLess(const void *x, const void *y);
int CompareNaturalOrderGreater(const void *x, const void *y);
int CompareIntegerNaturalOrderLess(const void *x, const void *y);
int CompareIntegerNaturalOrderGreater(const void *x, const void *y);
int CompareByMagnitudeLess(const void *x, const void *y);
int CompareByMagnitudeGreater(const void *x, const void *y);

// Needed by DMatrix::DGEES()
LOGICAL SelectLeftHalfPlane(double *WR, double *WI);
LOGICAL SelectRightHalfPlane(double *WR, double *WI);

// Needed by CMatrix::Partitioning()
struct Elm1 {
  INTEGER idx;
  double val;
};
int CompareElm1Idx(const void *x, const void *y); // Increasing order

// Needed by class Node
struct Elm2 {
  INTEGER idx;
  INTEGER above;
  INTEGER below;
};
int CompareElm2Idx(const void *x, const void *y);   // Increasing order
int CompareElm2Above(const void *x, const void *y); // Increasing order
int CompareElm2Below(const void *x, const void *y); // Increasing order

// Needed by CMatrix::BuildTreeDownward2()
struct Elm3 {
  INTEGER idx;
  INTEGER dim; // Number of grid points along a direction
  double len;  // Grid length along a direction
  double inc;  // Length between two consecutive grid points
  double low;  // Lower end of the grid
};
int CompareElm3Dim(const void *x, const void *y); // Increasing order

#ifdef USE_MATERN
// Modified Bessel function of the second kind. Needed by class IsotropicMatern
double besselk(double nu, double x);
#endif

//--------------------------------------------------------------------------
// BLAS and LAPACK
#ifdef USE_MKL
#define drotg_  drotg
#define drot_   drot
#define dcopy_  dcopy
#define ddot_   ddot
#define dnrm2_  dnrm2
#define dasum_  dasum
#define dgemv_  dgemv
#define dtrsv_  dtrsv
#define dtrsm_  dtrsm
#define dgemm_  dgemm
#define dlange_ dlange
#define dgesvx_ dgesvx
#define dposvx_ dposvx
#define dgelss_ dgelss
#define dsyev_  dsyev
#define dgees_  dgees
#define dgesvd_ dgesvd
#define dggev_  dggev
#define dgetrf_ dgetrf
#define dgetrs_ dgetrs
#define dgecon_ dgecon
#define dpotrf_ dpotrf
#define dpotrs_ dpotrs
#define dpocon_ dpocon
#define dgeqp3_ dgeqp3
#define dgeqrf_ dgeqrf
#define dorgqr_ dorgqr
#endif
extern "C" {

#ifdef USE_MATERN
  // Modified Bessel function of the second kind
  void zbesk_(double *zr,    // Input, real part
              double *zi,    // Input, imaginary part
              double *fnu,   // Order
              int *kode,     // Scaling option. 1 = no scaling
              int *n,        // Computes order from fnu to fnu + n -1
              double cyr[],  // Outputs, real part, for different orders
              double cyi[],  // Outputs, imaginary part, for different orders
              int *nz,       // Number of elements in cyr/cyi set to zero
              int *ierr      // Error flag
              );
#endif

  // BLAS 1
  void drotg_(double *DA, double *DB, double *C, double *S);
  void drot_(const INTEGER *N, double *DX, const INTEGER *INCX, double *DY, const INTEGER *INCY, const double *C, const double *S);
  void dcopy_(const INTEGER *N, const double *DX, const INTEGER *INCX, double *DY, const INTEGER *INCY);
  double ddot_(const INTEGER *N, const double *DX, const INTEGER *INCX, const double *DY, const INTEGER *INCY);
  double dnrm2_(const INTEGER *N, const double *X, const INTEGER *INCX);
  double dasum_(const INTEGER *N, const double *X, const INTEGER *INCX);

  // BLAS 2
  void dgemv_(const char *TRANS, const INTEGER *M, const INTEGER *N, const double *ALPHA, const double *A, const INTEGER *LDA, const double *X, const INTEGER *INCX, const double *BETA, double *Y, const INTEGER *INCY);
  void dtrsv_(const char *UPLO, const char *TRANS, const char *DIAG, const INTEGER *N, const double *A, const INTEGER *LDA, double *X, const INTEGER *INCX);
  void dtrsm_(const char *SIDE, const char *UPLO, const char *TRANSA, const char *DIAG, const INTEGER *M, const INTEGER *N, const double &ALPHA, const double *A, const INTEGER *LDA, double *B, const INTEGER *LDB);

  // BLAS 3
  void dgemm_(const char *TRANSA, const char *TRANSB, const INTEGER *M, const INTEGER *N, const INTEGER *K, const double *ALPHA, const double *A, const INTEGER *LDA, const double *B, const INTEGER *LDB, const double *BETA, double *C, const INTEGER *LDC);

  // LAPACK matrix norm
  double dlange_(const char *NORM, const INTEGER *M, const INTEGER *N, const double *A, const INTEGER *LDA, double *WORK);

  // LAPACK driver routines: Linear equations
  void dgesvx_(const char *FACT, const char *TRANS, const INTEGER *N, const INTEGER *NRHS, double *A, const INTEGER *LDA, double *AF, const INTEGER *LDAF, INTEGER *IPIV, char *EQUED, double *R, double *C, double *B, const INTEGER *LDB, double *X, const INTEGER *LDX, double *RCOND, double *FERR, double *BERR, double *WORK, INTEGER *IWORK, INTEGER *INFO);
  void dposvx_(const char *FACT, const char *UPLO, const INTEGER *N, const INTEGER *NRHS, double *A, const INTEGER *LDA, double *AF, const INTEGER *LDAF, char *EQUED, double *S, double *B, const INTEGER *LDB, double *X, const INTEGER *LDX, double *RCOND, double *FERR, double *BERR, double *WORK, INTEGER *IWORK, INTEGER *INFO);

  // LAPACK driver routines: Linear least squares
  void dgelss_(const INTEGER *M, const INTEGER *N, const INTEGER *NRHS, double *A, const INTEGER *LDA, double *B, const INTEGER *LDB, double *S, const double *RCOND, INTEGER *RANK, double *WORK, const INTEGER *LWORK, INTEGER *INFO);

  // LAPACK driver routines: Standard symmetric eigenproblems
  void dsyev_(const char *JOBZ, const char *UPLO, const INTEGER *N, double *A, const INTEGER *LDA, double *W, double *WORK, const INTEGER *LWORK, INTEGER *INFO);

  // LAPACK driver routines: Standard nonsymmetric eigenproblems
  void dgees_(const char *JOBVS, const char *SORT, LOGICAL (*SELECT)(double *WR, double *WI), const INTEGER *N, double *A, const INTEGER *LDA, INTEGER *SDIM, double *WR, double *WI, double *VS, const INTEGER *LDVS, double *WORK, const INTEGER *LWORK, LOGICAL *BWORK, INTEGER *INFO);

  // LAPACK driver routines: Standard singular value problems
  void dgesvd_(const char *JOBU, const char *JOBVT, const INTEGER *M, const INTEGER *N, double *A, const INTEGER *LDA, double *S, double *U, const INTEGER *LDU, double *VT, const INTEGER *LDVT, double *WORK, const INTEGER *LWORK, INTEGER *INFO);

  // LAPACK driver routines: Generalized nonsymmetric eigenproblems
  void dggev_(const char *JOBVL, const char *JOBVR, const INTEGER *N, double *A, const INTEGER *LDA, double *B, const INTEGER *LDB, double *ALPHAR, double *ALPHAI, double *BETA, double *VL, const INTEGER *LDVL, double *VR, const INTEGER *LDVR, double *WORK, const INTEGER *LWORK, INTEGER *INFO);

  // LAPACK computational routines: Linear equations
  void dgetrf_(const INTEGER *M, const INTEGER *N, double *A, const INTEGER *LDA, INTEGER *IPIV, INTEGER *INFO);
  //INTEGER dgetrf_(const INTEGER *M, const INTEGER *N, double *A, const INTEGER *LDA, INTEGER *IPIV, INTEGER *INFO);
  void dgetrs_(const char *TRANS, const INTEGER *N, const INTEGER *NRHS, const double *A, const INTEGER *LDA, const INTEGER *IPIV, double *B, const INTEGER *LDB, INTEGER *INFO);
  //INTEGER dgetrs_(const char *TRANS, const INTEGER *N, const INTEGER *NRHS, const double *A, const INTEGER *LDA, const INTEGER *IPIV, double *B, const INTEGER *LDB, INTEGER *INFO);
  void dgecon_(const char *NORM, const INTEGER *N, const double *A, const INTEGER *LDA, const double *ANORM, double *RCOND, double *WORK, INTEGER *IWORK, INTEGER *INFO);
  void dpotrf_(const char *UPLO, const INTEGER *N, double *A, const INTEGER *LDA, INTEGER *INFO);
  //INTEGER dpotrf_(const char *UPLO, const INTEGER *N, double *A, const INTEGER *LDA, INTEGER *INFO);
  void dpotrs_(const char *UPLO, const INTEGER *N, const INTEGER *NRHS, const double *A, const INTEGER *LDA, double *B, const INTEGER *LDB, INTEGER *INFO);
  //INTEGER dpotrs_(const char *UPLO, const INTEGER *N, const INTEGER *NRHS, const double *A, const INTEGER *LDA, double *B, const INTEGER *LDB, INTEGER *INFO);
  void dpocon_(const char *UPLO, const INTEGER *N, const double *A, const INTEGER *LDA, const double *ANORM, double *RCOND, double *WORK, INTEGER *IWORK, INTEGER *INFO);

  // LAPACK computational routines: Orthogonal factorizations and linear least squares
  void dgeqp3_(const INTEGER *M, const INTEGER *N, double *A, const INTEGER *LDA, INTEGER *JPVT, double *TAU, double *WORK, const INTEGER *LWORK, INTEGER *INFO);
  void dgeqrf_(const INTEGER *M, const INTEGER *N, double *A, const INTEGER *LDA, double *TAU, double *WORK, const INTEGER *LWORK, INTEGER *INFO);
  void dorgqr_(const INTEGER *M, const INTEGER *N, const INTEGER *K, double *A, const INTEGER *LDA, const double *TAU, double *WORK, const INTEGER *LWORK, INTEGER *INFO);

}


//--------------------------------------------------------------------------
// Templated new and delete. For all the delete functions, the pointer
// is set to NULL. For all the new functions, the allocated memory is
// set to zero if the data type is numeric.
//
// Provides the following interfaces:
//
// 1a. Delete a 1D array
//   void Delete_1D_Array(T **a);
//
// 1b. New a 1D array of length d1
//   void New_1D_Array(T **a, T1 d1);
//
// 2a. Delete a 2D array (first dimension d1)
//   void Delete_2D_Array(T ***a, T1 d1);
//
// 2b. New a 2D array of size d1*d2
//   void New_2D_Array(T ***a, T1 d1, T2 d2);
//
// 3a. Delete a 3D array (first dimension d1, second dimension d2)
//   void Delete_3D_Array(T ****a, T1 d1, T2 d2);
//
// 3b. New a 3D array of size d1*d2*d3
//   void New_3D_Array(T ****a, T1 d1, T2 d2, T3 d3);

template <class T>
void Delete_1D_Array(T **a) {
  delete [] (*a);
  (*a) = NULL;
}

template <class T, class T1>
void New_1D_Array(T **a, T1 d1) {
  (*a) = new T [d1];
  IF_NUMERIC_TYPE{ memset((*a), 0, d1*sizeof(T)); }
}

template <class T, class T1>
void Delete_2D_Array(T ***a, T1 d1) {
  if (*a) {
    for (T1 i = 0; i < d1; i++) {
      delete [] (*a)[i];
    }
    delete [] (*a);
    (*a) = NULL;
  }
}

template <class T, class T1, class T2>
void New_2D_Array(T ***a, T1 d1, T2 d2) {
  (*a) = new T * [d1];
  for (T1 i = 0; i < d1; i++) {
    (*a)[i] = new T [d2];
    IF_NUMERIC_TYPE{ memset((*a)[i], 0, d2*sizeof(T)); }
  }
}

template <class T, class T1, class T2>
void Delete_3D_Array(T ****a, T1 d1, T2 d2) {
  if (*a) {
    for (T1 i = 0; i < d1; i++) {
      if ((*a)[i]) {
        for (T2 j = 0; j < d2; j++) {
          delete [] (*a)[i][j];
        }
        delete [] (*a)[i];
      }
    }
    delete [] (*a);
    (*a) = NULL;
  }
}

template <class T, class T1, class T2, class T3>
void New_3D_Array(T ****a, T1 d1, T2 d2, T3 d3) {
  (*a) = new T ** [d1];
  for (T1 i = 0; i < d1; i++) {
    (*a)[i] = new T * [d2];
    for (T2 j = 0; j < d2; j++) {
      (*a)[i][j] = new T [d3];
      IF_NUMERIC_TYPE{ memset((*a)[i][j], 0, d3*sizeof(T)); }
    }
  }
}


//--------------------------------------------------------------------------
// Templated swap
template <class T>
void Swap(T &a, T &b) {
  T tmp = a;
  a = b;
  b = tmp;
}


//--------------------------------------------------------------------------
// Templated max
template <class T>
T Max(T a, T b) {
  return a >= b ? a : b;
}


//--------------------------------------------------------------------------
// Templated min
template <class T>
T Min(T a, T b) {
  return a < b ? a : b;
}


#endif
