// The GMRES class implements the preconditioned GMRES method for
// solving a linear system of equations
//
//   Ax = b
//
// by using an initial guess x0 and a preconditioner M (that
// approximates the inverse of A).
//
// The implementation of this class is parallelized because the
// underlying linear algebra is parallelized.

#ifndef _GMRES_
#define _GMRES_

#include "../Matrices/DMatrix.hpp"

#define MAX_REFINE_STEP 5
enum SolveType { USUAL, REFINE };

class GMRES {

public:

  GMRES();
  GMRES(const GMRES &G);
  GMRES& operator= (const GMRES &G);
  ~GMRES();

  void Init(void);
  void ReleaseAllMemory(void);
  void DeepCopy(const GMRES &G);

  // Solve.
  // MatrixA and MatrixM must have the following method:
  //   void MatVec(const DVector &b, DVector &y, MatrixMode ModeA,
  //               MatState mState);
  // Only ModeA = NORMAL is used.
  //
  // Explanation: The operators A and M are used in the form of
  // matrix-vector products. In the usual case these are indeed
  // matrix-vector multiplies. However, in some scenarios
  // (particularly for the preconditioner M), the matrix-vector
  // products are indeed linear solves. In such scenarios, M has
  // already been factorized and hence mState should be either LU_FACT
  // or CHOL_FACT, depending on whether M is SPD. Then, M.MatVec()
  // knows what operation to take.
  template<class MatrixA, class MatrixM>
  void Solve(MatrixA &A,        // Martix A
             MatState mStateA,  // Whether A is factorized
             const DVector &b,  // Right-hand side b
             const DVector &x0, // Initial guess x0
             MatrixM &M,        // Preconditioner M approx inv(A)
             MatState mStateM,  // Whether M is factorized
             INTEGER m,         // Restart cycle m
             INTEGER MaxIt,     // Maximum # of iterations
             double RTol        // Relative residual tolerance
             );

  // This function refines the initial solution vector x0 =
  // M*b. Assuming that M is sufficiently close to inv(A) [such that
  // x0 is sufficiently close to inv(A)*b], the refinement solve
  // terminates when the residual > half of the previous residual, or
  // when GMRES has iterated MAX_REFINE_STEP steps.
  template<class MatrixA, class MatrixM>
  void RefineSolve(MatrixA &A,        // Martix A
                   MatState mStateA,  // Whether A is factorized
                   const DVector &b,  // Right-hand side b
                   const DVector &x0, // Initial guess x0
                   MatrixM &M,        // Preconditioner M approx inv(A)
                   MatState mStateM   // Whether M is factorized
                   );

  // Get the norm of the right hand side b.
  double GetNormRHS(void) const;

  // Get the solution vector x
  void GetSolution(DVector &Sol) const;

  // Get the pointer to the array of residual history (NOTE: not
  // relative residuals). Iter is the number of iterations. That is,
  // the residuals are stored in [0 .. Iter-1].
  const double* GetResHistory(INTEGER &Iter) const;

protected:

private:

  double NormB;
  DVector x;
  INTEGER mIter;
  double *mRes;

  // Common subroutine called by Solve() and RefinementSolve()
  template<class MatrixA, class MatrixM>
  void SolveBase(MatrixA &A,        // Martix A
                 MatState mStateA,  // Whether A is factorized
                 const DVector &b,  // Right-hand side b
                 const DVector &x0, // Initial guess x0
                 MatrixM &M,        // Preconditioner M approx inv(A)
                 MatState mStateM,  // Whether M is factorized
                 INTEGER m,         // Restart cycle m
                 INTEGER MaxIt,     // Maximum # of iterations
                 double RTol,       // Relative residual tolerance
                 SolveType Type     // Normal solve or refinement solve?
                 );

};

#include "GMRES.tpp"

#endif
