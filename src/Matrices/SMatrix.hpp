// The SMatrix class implements a sparse matrix in double
// precision. The matrix is stored in the compressed sparse row (CSR)
// format, which is in row-major order.
//
//  * M is the number of rows;
//  * N is the number of columns;
//  * nnz is the number of nonzeros;
//  * the array start stores the starting location of each row;
//  * the array idx stores the column index of each nonzero; and
//  * the array A stores the nonzeros.
//
// The column indices for each row must be in the increasing order
// without duplicates. The array start has a length M+1, where
// start[M] is always nnz. The arrays idx and A have a length nnz.
//
// Currently this class supports matrix-vector multiplications and the
// extraction of the diagonal. The former is NOT parallelized due to
// its memory-bound nature. The latter is parallelized.

#ifndef _SMATRIX_
#define _SMATRIX_

#include "DMatrix.hpp"
#include "../Misc/spblas.hpp"

class SMatrix {

public:

  SMatrix();
  SMatrix(INTEGER N_, INTEGER nnz_);
  SMatrix(INTEGER M_, INTEGER N_, INTEGER nnz_);
  SMatrix(const SMatrix &G);
  SMatrix& operator= (const SMatrix &G);
  ~SMatrix();

  void Init(void);
  void Init(INTEGER N_, INTEGER nnz_);             // Not a zero matrix
  void Init(INTEGER M_, INTEGER N_, INTEGER nnz_); // Not a zero matrix
  void ReleaseAllMemory(void);
  void DeepCopy(const SMatrix &G);

  //-------------------- Utilities --------------------
  //
  // A is the self matrix.

  // Get dimension
  INTEGER GetM(void) const;
  INTEGER GetN(void) const;

  // Get nnz
  INTEGER GetNNZ(void) const;

  // Get the pointer to start
  INTEGER* GetPointerStart(void) const;

  // Get the pointer to idx
  INTEGER* GetPointerIdx(void) const;

  // Get the pointer to A
  double* GetPointerA(void) const;

  // B = A(RowStart:RowStart+nRow-1, ColStart:ColStart+nCol-1)
  void GetBlock(INTEGER RowStart, INTEGER nRow, INTEGER ColStart, INTEGER nCol,
                DMatrix &B) const;
  // B = A(IdxRow, IdxCol)
  void GetBlock(INTEGER *IdxRow, INTEGER nRow, INTEGER *IdxCol, INTEGER nCol,
                DMatrix &B) const;
  // B = A(IdxRow, ColStart:ColStart+nCol-1)
  void GetBlock(INTEGER *IdxRow, INTEGER nRow, INTEGER ColStart, INTEGER nCol,
                DMatrix &B) const;
  // Note that sparse matrix is row-major order whereas dense matrix
  // is column-major order.

  //-------------------- Matrix computations --------------------
  //
  // A is the self matrix.

  // A = A / b  or  C = A / b  (b is scalar)
  void Divide(double b);
  void Divide(double b, SMatrix &C) const;

  // A = diag(1./b) * A * diag(1./b)  or  C = diag(1./b) * A * diag(1./b)
  void SymmetricDivide(const DVector &b);
  void SymmetricDivide(const DVector &b, SMatrix &C) const;

  // y = mode(A) * b
  //
  // NOTE: The argument mState is redundant; its appearance is only to
  // conform with the interface requirement of linear solvers. The
  // function will return error if one passes to the argument other
  // than UNFACT.
  void MatVec(const DVector &b, DVector &y, MatrixMode ModeA,
              MatState mState = UNFACT) const;

  // C = mode(A) * mode(B)
  void MatMat(const DMatrix &B, DMatrix &C,
              MatrixMode ModeA, MatrixMode ModeB) const;
  void MatMat(const SMatrix &B, DMatrix &C,
              MatrixMode ModeA, MatrixMode ModeB) const;

  // b = diag(A)
  void Diag(DVector &b) const;
  void Diag(double *b) const;

protected:

private:

  INTEGER M, N;   // Matrix dimension
  INTEGER nnz;    // Number of nonzeros
  INTEGER *start; // Start locations of a row
  INTEGER *idx;   // Column indices of the nonzeros
  double *A;      // Nonzeros

};

#endif
