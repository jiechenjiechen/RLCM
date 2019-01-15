// The DMatrix class implements a dense matrix in double precision. A
// substantial portion of the methods in this class are wrappers of
// BLAS2, BLAS3, and LAPACK routines. The interface tends to be
// coherent with Matlab constructs. The matrix data is stored by using
// column-major order.
//
// In addition, this class provides methods for building the kernel
// matrix from point sets.
//
// The implementation of this class is parallelized.

#ifndef _DMATRIX_
#define _DMATRIX_

#include "DVector.hpp"

class DMatrix {

public:

  DMatrix();
  DMatrix(INTEGER N_);
  DMatrix(INTEGER M_, INTEGER N_);
  DMatrix(const DMatrix &G);
  DMatrix& operator= (const DMatrix &G);
  ~DMatrix();

  void Init(void);
  void Init(INTEGER N_);             // Initialized as a zero matrix
  void Init(INTEGER M_, INTEGER N_); // Initialized as a zero matrix
  void ReleaseAllMemory(void);
  void DeepCopy(const DMatrix &G);

  //-------------------- Utilities --------------------
  //
  // A is the self matrix.

  // Get dimension
  INTEGER GetM(void) const;
  INTEGER GetN(void) const;

  // Get A(i,j)
  double GetEntry(INTEGER i, INTEGER j) const;
  // b = A(:,i)
  void GetColumn(INTEGER i, DVector &b) const;
  // B = A(:,idx). n is the length of idx
  void GetColumns(INTEGER *idx, INTEGER n, DMatrix &B) const;
  // b = A(i,:)
  void GetRow(INTEGER i, DVector &b) const;
  // B = A(RowStart:RowStart+nRow-1, ColStart:ColStart+nCol-1)
  void GetBlock(INTEGER RowStart, INTEGER nRow, INTEGER ColStart, INTEGER nCol,
                DMatrix &B) const;
  // Get the double* pointer
  double* GetPointer(void) const;

  // Set A(i,j) = b
  void SetEntry(INTEGER i, INTEGER j, double b);
  // A(i,j) += b
  void AddToEntry(INTEGER i, INTEGER j, double b);
  // A(:,i) = b
  void SetColumn(INTEGER i, const DVector &b);
  // A(i,:) = b
  void SetRow(INTEGER i, const DVector &b);
  // A(RowStart:RowStart+nRow-1, ColStart:ColStart+nCol-1) = B
  void SetBlock(INTEGER RowStart, INTEGER nRow, INTEGER ColStart, INTEGER nCol,
                const DMatrix &B);
  void SetBlock(INTEGER RowStart, INTEGER nRow, INTEGER ColStart, INTEGER nCol,
                const double *B);

  // A = I
  void SetIdentity(void);
  // A = lambda*I
  void SetMultIdentity(double lambda);
  // A = c
  void SetConstVal(double c);
  // A = rand()
  void SetUniformRandom01(void);
  // A = randn()
  void SetStandardNormal(void);
  // A = rtnd(1); each element is a student-t of degree 1
  void SetStudentT1(void);
  // Each column of A is a multivariate student-t of degree 1
  void SetMultivariateStudentT1(void);
  // A = random('sech')
  void SetRandomSech(void);
  // A = diag(b)
  void MakeDiag(const DVector &b);

  // Print the matrix in the Matlab form
  void PrintMatrixMatlabForm(const char *name) const;

  //-------------------- Matrix computations --------------------
  //
  // A is the self matrix.

  // A = A'  or  B = A'
  void Transpose(void);
  void Transpose(DMatrix &B) const;

  // A = (A+A')/2  or  B = (A+A')/2
  void Symmetrize(void);
  void Symmetrize(DMatrix &B) const;

  // A = -A  or  B = -A
  void Negate(void);
  void Negate(DMatrix &B) const;

  // A = A + B  or  C = A + B
  void Add(const DMatrix &B);
  void Add(const DMatrix &B, DMatrix &C) const;

  // A = A - B  or  C = A - B
  void Subtract(const DMatrix &B);
  void Subtract(const DMatrix &B, DMatrix &C) const;

  // A = A + b  or  C = A + b  (b is scalar)
  void Add(double b);
  void Add(double b, DMatrix &C) const;

  // A = A - b  or  C = A - b  (b is scalar)
  void Subtract(double b);
  void Subtract(double b, DMatrix &C) const;

  // A = A .* B  or  C = A .* B  (elementwise)
  void Multiply(const DMatrix &B);
  void Multiply(const DMatrix &B, DMatrix &C) const;

  // A = A ./ B  or  C = A ./ B  (elementwise)
  void Divide(const DMatrix &B);
  void Divide(const DMatrix &B, DMatrix &C) const;

  // A = A * b  or  C = A * b  (b is scalar)
  void Multiply(double b);
  void Multiply(double b, DMatrix &C) const;

  // A = A / b  or  C = A / b  (b is scalar)
  void Divide(double b);
  void Divide(double b, DMatrix &C) const;

  // A = A + sI  or  B = A + sI
  void AddDiagonal(double s);
  void AddDiagonal(double s, DMatrix &B) const;

  // A = A - sI  or  B = A - sI
  void SubtractDiagonal(double s);
  void SubtractDiagonal(double s, DMatrix &B) const;

  // A = cos(A)  or  B = cos(A)
  void Cos(void);
  void Cos(DMatrix &B) const;

  // A = sqrt(A)  or  B = sqrt(A)
  void Sqrt(void);
  void Sqrt(DMatrix &B) const;

  // A = log(A)  or  B = log(A)
  void Log(void);
  void Log(DMatrix &B) const;

  // b = sum(A, dim).  dim starts with 1
  void Sum(DVector &b, INTEGER dim) const;

  // b = prod(A, dim).  dim starts with 1
  void Prod(DVector &b, INTEGER dim) const;

  // A = b * c'
  void OuterProduct(const DVector &b, const DVector &c);

  // y = mode(A) * b
  //
  // NOTE: The argument mState is redundant; its appearance is only to
  // conform with the interface requirement of linear solvers. The
  // function will return error if one passes to the argument other
  // than UNFACT.
  void MatVec(const DVector &b, DVector &y, MatrixMode ModeA,
              MatState mState = UNFACT) const;

  // y = alpha * mode(A) * b + beta * y
  void DGEMV(const DVector &b, DVector &y,
             double alpha, double beta, MatrixMode ModeA) const;

  // C = mode(A) * mode(B)
  void MatMat(const DMatrix &B, DMatrix &C,
              MatrixMode ModeA, MatrixMode ModeB) const;

  // C = alpha * mode(A) * mode(B) + beta * C
  void DGEMM(const DMatrix &B, DMatrix &C, double alpha, double beta,
             MatrixMode ModeA, MatrixMode ModeB) const;

  // Linear system solve: X = mode(A)\B
  void Mldivide(const DVector &b, DVector &x,
                MatrixMode ModeA, MatrixType TypeA);
  void Mldivide(const DMatrix &B, DMatrix &X,
                MatrixMode ModeA, MatrixType TypeA);
  // Least squares solve: X = A\B
  // Res must be either NULL or an allocated array with length size(B,2).
  void Mldivide(const DVector &b, DVector &x, double *Res = NULL);
  void Mldivide(const DMatrix &B, DMatrix &X, double *Res = NULL);

  // The following four routines are related to linear system
  // solves. They can be used separately from Mldivide in places where
  // efficiency is fine grained.
  //
  // Factorize the matrix. The data array A is destroyed to hold the
  // LU factors.
  void DGETRF(void);
  //
  // Solve linear systems Ax = b by using the factorization.
  void DGETRS(const DVector &b, DVector &x, MatrixMode ModeA) const;
  void DGETRS(const DMatrix &B, DMatrix &X, MatrixMode ModeA) const;
  //
  // Estimate the reciprocal condition number (1-norm) of the original
  // matrix A.
  double DGECON(void) const;

  // The above four routines are for general matrices. The following
  // routines are for symmetric positive definite matrices instead.
  //
  // Factorize the matrix. The data array A is destroyed to hold the
  // Cholesky factors.
  void DPOTRF(TriUPLO mTri);
  //
  // Solve linear systems Ax = b by using the factorization.
  void DPOTRS(const DVector &b, DVector &x, TriUPLO mTri) const;
  void DPOTRS(const DMatrix &B, DMatrix &X, TriUPLO mTri) const;
  // On the other hand, one may also choose to solve one of the
  // triangular systems rather than both.
  void DTRSV(const DVector &b, DVector &x,
             MatrixMode ModeA, TriUPLO mTri) const;
  void DTRSM(const DMatrix &B, DMatrix &X,
             MatrixMode ModeA, TriUPLO mTri) const;
  //
  // Estimate the reciprocal condition number (1-norm) of the original
  // matrix A.
  double DPOCON(TriUPLO mTri) const;

  // A = GG'
  void Chol(DMatrix &G, TriUPLO mTri) const;

  // lambda = eig(A), A symmetric. Eigenvalues are gauranteed to be in
  // ascending order.
  void SymEig(DVector &lambda) const;
  // [V,lambda] = eig(A)
  void SymEig(DVector &lambda, DMatrix &V) const;

  // lambda = eig(A,B), A,B symmetric, but B is not definite. In such
  // a case, only nonsymmetric algorithm can be used. This function is
  // a special hack in the square root factorization of a CMatrix, in
  // the sense that the eigenvalues are known to be real.
  void SymEigBIndef(const DMatrix &B, DVector &lambda) const;
  // [V,lambda] = eig(A,B)
  void SymEigBIndef(const DMatrix &B, DVector &lambda, DMatrix &V) const;

  // [U,T] = schur(A, 'real')
  void RealSchur(DMatrix &T, DMatrix &U, SelectType type) const;

  // B = orth(A). A must be square
  void Orth(DMatrix &B) const;

  // Perform QR with column pivoting, return the first k pivots. That
  // is, pivots[0] is the index of the original matrix column that is
  // pivoted to the first column under the QR factorization
  void QRpivots(INTEGER k, INTEGER *pivots) const;

  // rank(A)
  INTEGER Rank(void) const;

  // cond(A,2)
  double Cond2(void) const;

  // norm(A,2)
  double Norm2(void) const;

  // norm(A,'fro')
  double NormF(void) const;

  // b = diag(A)
  void Diag(DVector &b) const;

  // trace(A) = sum(diag(A))
  double Tr(void) const;

  // Determinant. Only the following combinations are valid:
  //   TypeA = GENERAL, mState = UNFACT
  //   TypeA = GENERAL, mState = LU_FACT
  //   TypeA = PSD,     mState = UNFACT
  //   TypeA = PSD,     mState = CHOL_FACT
  LogDet Det(MatrixType TypeA, MatState mState);

  //-------------------- Build Kernel Matrix --------------------
  //
  // Build a kernel matrix from point set X (and Y). The class Kernel
  // must have the following methods:
  //
  //   template<class Point>
  //   double Eval(const Point &x, const Point &y) const;
  //   bool IsSymmetric(void) const;
  //
  // The class PointArray must have the following methods:
  //
  //   INTEGER GetN(void) const;
  //   void GetPoint(INTEGER i, Point &x) const;
  template<class Kernel, class Point, class PointArray>
  void BuildKernelMatrix(const Kernel &mKernel, const PointArray &X,
                         double lambda = 0.0);
  template<class Kernel, class Point, class PointArray>
  void BuildKernelMatrix(const Kernel &mKernel, const PointArray &X,
                         const PointArray &Y, double lambda = 0.0);

protected:

private:

  INTEGER M, N;  // Matrix dimension
  double *A;     // Matrix data in column major order

  INTEGER *IPIV; // Needed by DGETRF, DGETRS, and DGECON. Length = min(M,N)

  // Called by Mldivide()
  void Mldivide_LinearSystemSolve(const DMatrix &B, DMatrix &X,
                                  MatrixMode ModeA, MatrixType TypeA);
  void Mldivide_LeastSquaresSolve(const DMatrix &B, DMatrix &X,
                                  double *Res);

  // Called by SymEig()
  void DSYEV(char JOBZ, DVector &lambda, DMatrix &V) const;

  // Called by SymEigBIndef()
  void DGGEV(char JOBVR, const DMatrix &B, DVector &lambda, DMatrix &V) const;

  // Called by RealSchur()
  void DGEES(DMatrix &T, DMatrix &U,
             LOGICAL (*Select)(double *WR, double *WI)) const;

  // Called by Rank(), Cond2() and Norm2(). Singular values
  // only. Memory of S must be preallocated.
  void DGESVD(double *S) const;

};

#include "DMatrix.tpp"

#endif
