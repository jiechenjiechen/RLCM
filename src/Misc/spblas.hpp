// Sparse linear algebra for real numbers in double precision.
//
// A sparse vector x is represented by 4 quantities:
//  * n is the length;
//  * nnz is the number of nonzeros;
//  * the array idx stores the indicies of the nonzeros; and
//  * the array x stores the nonzeros.
// The indices must be in the increasing order without
// duplicates. Both arrays idx and x have a length nnz.
//
// A sparse matrix A is stored in the compressed sparse row (CSR)
// format:
//  * m is the number of rows;
//  * n is the number of columns;
//  * nnz is the number of nonzeros;
//  * the array start stores the starting location of each row;
//  * the array idx stores the column index of each nonzero; and
//  * the array A stores the nonzeros.
// The column indices for each row must be in the increasing order
// without duplicates. The array start has a length m+1, where
// start[m] is always nnz. The arrays idx and A have a length nnz.
//
// A dense matrix is stored in the column-major order.
//
// Only matrix-vector and matrix-matrix multiplications are
// parallelized; other functions ARE NOT. One reason is that in this
// particular software, the sparse vectors are used as sparse points,
// the computation of which more likely happens in the context of a
// point set. Parallelism should better be applied on iterating the
// points rather than inside the computation of each point. Also, note
// that sparse matrix-vector multiplications are memory-bound; so do
// not expect too much speed up.
//
// Naming convention:
//  * 'sp_' denotes sparse;
//  * 'd_' means the second argument is dense;
//  * 's_' means the second argument is sparse;
//  * 'n_' means no second argument; and
//  * the rest of the characters follows the naming convection of BLAS
//    if the functionality is similar to a BLAS routine, or is just a
//    novel name of the functionality.

#ifndef _SPBLAS_
#define _SPBLAS_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef USE_LONG
typedef long INTEGER;
#else
typedef int INTEGER;
#endif

// y = x, x sparse, y dense
void sp_n_ds2d(INTEGER N, INTEGER NNZX, INTEGER *IDXX, double *X, double *Y);

// z = dot(x,y), x sparse, y dense
double sp_d_ddot(INTEGER N, INTEGER NNZX, INTEGER *IDXX, double *X, double *Y);

// z = dot(x,y), x sparse, y sparse
double sp_s_ddot(INTEGER N, INTEGER NNZX, INTEGER *IDXX, double *X, INTEGER NNZY, INTEGER *IDXY, double *Y);

// z = func(x-y), x sparse, y dense, func is an elementwise function that eventually sums up
double sp_d_ddist(INTEGER N, INTEGER NNZX, INTEGER *IDXX, double *X, double *Y, double(*FUNC)(double));

// z = func(x-y), x sparse, y sparse, func is an elementwise function that eventually sums up
double sp_s_ddist(INTEGER N, INTEGER NNZX, INTEGER *IDXX, double *X, INTEGER NNZY, INTEGER *IDXY, double *Y, double(*FUNC)(double));

// z = func((x-y)/s), x sparse, y dense, s dense, func is an elementwise function that eventually sums up
double sp_d_ddists(INTEGER N, INTEGER NNZX, INTEGER *IDXX, double *X, double *Y, double *S, double(*FUNC)(double));

// z = func((x-y)/s), x sparse, y sparse, s dense, func is an elementwise function that eventually sums up
double sp_s_ddists(INTEGER N, INTEGER NNZX, INTEGER *IDXX, double *X, INTEGER NNZY, INTEGER *IDXY, double *Y, double *S, double(*FUNC)(double));

// z = sum(2*x.*y./(x+y)), used only for x,y >= 0.
// If x[i]+y[i]=0, then 2*x[i]*y[i]/(x[i]+y[i])=0.
// x sparse, y dense
double sp_d_dchi2(INTEGER N, INTEGER NNZX, INTEGER *IDXX, double *X, double *Y);

// z = sum(2*x.*y./(x+y)), used only for x,y >= 0.
// If x[i]+y[i]=0, then 2*x[i]*y[i]/(x[i]+y[i])=0.
// x sparse, y sparse
double sp_s_dchi2(INTEGER N, INTEGER NNZX, INTEGER *IDXX, double *X, INTEGER NNZY, INTEGER *IDXY, double *Y);

// [parallel] y = mode(A)*x, A sparse, x dense, y dense
void sp_d_dgemv(char TRANS, INTEGER M, INTEGER N, INTEGER NNZA, INTEGER *STARTA, INTEGER *IDXA, double *A, double *X, double *Y);

// [parallel] y = mode(A)*x, A sparse, x sparse, y dense
void sp_s_dgemv(char TRANS, INTEGER M, INTEGER N, INTEGER NNZA, INTEGER *STARTA, INTEGER *IDXA, double *A, INTEGER NNZX, INTEGER *IDXX, double *X, double *Y);

// [parallel] C = mode(A)*mode(B), A sparse, B dense, C dense
void sp_d_dgemm(char TRANSA, char TRANSB, INTEGER M, INTEGER N, INTEGER K, INTEGER NNZA, INTEGER *STARTA, INTEGER *IDXA, double *A, double *B, double *C);

// [parallel] C = mode(A)*mode(B), A sparse, B sparse, C dense
void sp_s_dgemm(char TRANSA, char TRANSB, INTEGER M, INTEGER N, INTEGER K, INTEGER NNZA, INTEGER *STARTA, INTEGER *IDXA, double *A, INTEGER NNZB, INTEGER *STARTB, INTEGER *IDXB, double *B, double *C);

#endif
