// The Lanczos class implements the Lanczos procedure for a symmetric
// matrix A:
//
//   A V = V T + ...
//
// The class supports computing Ritz values and vectors:
//
//   T = W S W',  U = V W.
//
// The user may choose whether to invoke partial reorthogonalization.
//
// For practical purposes, the Lanczos procedure is run either with a
// prescribed number of iterations, or with early stopping when the
// top k Ritz values have converged.
//
// The implementation of this class is parallelized because the
// underlying linear algebra is parallelized.

#ifndef _LANCZOS_
#define _LANCZOS_

#include "../Matrices/DMatrix.hpp"

class Lanczos {

public:

  Lanczos();
  Lanczos(const Lanczos &G);
  Lanczos& operator= (const Lanczos &G);
  ~Lanczos();

  void Init(void);
  void ReleaseAllMemory(void);
  void DeepCopy(const Lanczos &G);

  // Lanczos procedure.
  // MatrixA must have the following method:
  //   void MatVec(const DVector &b, DVector &y, MatrixMode ModeA,
  //               MatState mState);
  // Only ModeA = NORMAL is used.
  template<class MatrixA>
  void Run(MatrixA &A,         // Martix A
           MatState mStateA,   // Needed by A.MatVec()
           const DVector &v,   // Initial vector
           bool PartialReorth, // Whether to perform partial reorth
           INTEGER MaxIt,      // Maximum # of iterations
           bool EarlyStop,     // Whether do early stopping
           INTEGER k,          // Early stopping when k Ritz values converged
           double RTol         // Convergence test for Ritz values
           );

  // Get the norm of the initial vector v
  double GetNormV(void) const;

  // Get the number of iterations
  INTEGER GetIter(void) const;

  // Get tridiagonal matrix T
  void GetT(DMatrix &T_) const;

  // Get subspace V (number of columns: Iter)
  void GetV(DMatrix &V_) const;

  // Get Ritz values S. Undefined if EarlyStop is false.
  void GetRitzValues(DVector &S_) const;

  // Get Ritz vectors U (number of columns: Iter; only the top k
  // converged). Undefined if EarlyStop is false.
  void GetRitzVectors(DMatrix &U_) const;

  // Get the pointer to the array of residuals
  //
  //   Res[i] = norm( A * U_i - S_i * U_i ).
  const double* GetRes(void) const;

protected:

private:

  double NormV;
  INTEGER mIter;
  double *mRes;
  DMatrix T;
  DMatrix V;
  DVector S;
  DMatrix U;

};

#include "Lanczos.tpp"

#endif
