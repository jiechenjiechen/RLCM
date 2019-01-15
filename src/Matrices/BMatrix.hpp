// The BMatrix class implements a block diagonal matrix in double
// precision. The matrix is SYMMETRIC but not necessarily POSITIVE
// DEFINITE.
//
// The implementation of this class is parallelized.

#ifndef _BMATRIX_
#define _BMATRIX_

#include "Node.hpp"
#ifdef HAS_METIS
#include "SMatrix.hpp"
#include <metis.h>
#endif

class BMatrix {

public:

  BMatrix();
  BMatrix(const BMatrix &G);
  BMatrix& operator= (const BMatrix &G);
  ~BMatrix();

  void Init(void);
  void ReleaseAllMemory(void);
  void DeepCopy(const BMatrix &G);

  //-------------------- Utilities --------------------
  //
  // A is the self matrix.

  // Get matrix size
  INTEGER GetN(void) const;

  // Get root
  Node* GetRoot(void) const;

  // Get an estimation of the memory consumption of the matrix.
  INTEGER GetMemEst(void) const;

  // Inspect tree structure
  void PrintTree(void) const;

  // Convert to dense matrix B = (DMatrix)A. This routine requires
  // O(N^2) memory. Use it wisely.
  //
  // Normally, one calls this routine with mState = UNFACT. However,
  // if the matrix results from Invert(), one should do either of the
  // following:
  //
  // (i)  If Invert() is called with TypeA = GENERAL, call this routine
  //      with mState = LU_FACT.
  // (ii) If Invert() is called with TypeA = SPD, call this routine
  //      with mState = CHOL_FACT.
  //
  // The literal meaning is that if tA = inv(A), then
  // tA.ConvertToDMatrix() computes the dense matrix form of tA.
  void ConvertToDMatrix(DMatrix &B, MatState mState) const;

  //-------------------- Matrix computations --------------------
  //
  // A is the self matrix.

  // y = mode(A) * b
  //
  // Normally, one calls this routine with mState = UNFACT. However,
  // if the matrix results from Invert(), one should do either of the
  // following:
  //
  // (i)  If Invert() is called with TypeA = GENERAL, call this routine
  //      with mState = LU_FACT.
  // (ii) If Invert() is called with TypeA = SPD, call this routine
  //      with mState = CHOL_FACT.
  //
  // The literal meaning is that if tA = inv(A), then tA.MatVec()
  // computes the matrix-vector product with tA. Equivalently, this
  // means solving a linear system with A.
  void MatVec(const DVector &b, DVector &y, MatrixMode ModeA,
              MatState mState) const;

  // Y = mode(A) * B
  //
  // Normally, one calls this routine with mState = UNFACT. However,
  // if the matrix results from Invert(), one should do either of the
  // following:
  //
  // (i)  If Invert() is called with TypeA = GENERAL, call this routine
  //      with mState = LU_FACT.
  // (ii) If Invert() is called with TypeA = SPD, call this routine
  //      with mState = CHOL_FACT.
  //
  // The literal meaning is that if tA = inv(A), then tA.MatMat()
  // computes the matrix-matrix product with tA. Equivalently, this
  // means solving a linear system with A.
  void MatMat(const DMatrix &B, DMatrix &Y, MatrixMode ModeA,
              MatState mState) const;

  // tA = inv(A). The diagonal blocks of tA stores the factorizations
  // only; they are not the actual inverse. When TypeA is SPD, the
  // upper triangular factor is stored.
  void Invert(BMatrix &tA, MatrixType TypeA) const;

  // Determinant. Only the following combinations are valid:
  //   TypeA = GENERAL, mState = UNFACT
  //   TypeA = GENERAL, mState = LU_FACT
  //   TypeA = PSD,     mState = UNFACT
  //   TypeA = PSD,     mState = CHOL_FACT
  //
  // Unlike ConvertToDMatrix(), MatVec(), and MatMat(), one cannot use
  // this routine to compute the determinant of a matrix resulting
  // from Invert(). If tA = inv(A), then tA.Det() computes the
  // determinant of A, because tA contains the factored form of
  // A. tA.Det() WILL NOT COMPUTE THE DETERMINENT OF tA.
  LogDet Det(MatrixType TypeA, MatState mState) const;

  //--------------- Build the matrix from points ---------------

  // Build a block-diagonal kernel matrix from point array X. The
  // hierarchy tree is binary. The resulting matrix is SPD.
  //
  // The template classes Kernel, Point and PointArray are subject to
  // the same requirements as those for
  // DMatrix::BuildKernelMatrix. Additionally, the template class
  // PointArray must have the following method:
  //
  //   void GetSubset(INTEGER start, INTEGER n, DPointArray &Y) const;
  //   void GetSubset(INTEGER *idx, INTEGER n, DPointArray &Y) const;
  //   INTEGER RandomBipartition(INTEGER start, INTEGER n, INTEGER N0,
  //                             INTEGER *perm, DPoint &normal, double &offset);
  //   INTEGER PCABipartition(INTEGER start, INTEGER n, INTEGER N0,
  //                          INTEGER *perm, DPoint &normal, double &offset);
  //
  // For an exception, when the partitioning method mPar == BBOX
  // (which is appropriate only when the points are low-dimensional
  // and the PointArray class is DPointArray), the above functions are
  // not called. Instead, what is called is
  //
  //  INTEGER BBoxBipartition(INTEGER start, INTEGER n, INTEGER N0,
  //                          INTEGER *perm, DPoint &normal, double &offset,
  //                          const double *bbox, INTEGER &which_dim);
  //
  // The construction consists of two steps: First build the tree and
  // then build the matrix. The tree construction step is done once
  // for all, as long as the points X do not change. X will be
  // permuted. The matrix construction step may take different kernels
  // and kernel parameters as input.
  //
  // The arrays Perm and iPerm are used to record the permutation in
  // the tree construction step. Perm can be either NULL or
  // preallocated with sufficient memory. If Perm is NULL, iPerm is
  // not referenced. If Perm is preallocated, iPerm can be either NULL
  // or preallocated with sufficient memory.
  template<class Point, class PointArray>
  void BuildTree(PointArray &X,  // Points; will be permuted
                 INTEGER *Perm,  // Xnew[i] = Xold[Perm[i]]
                 INTEGER *iPerm, // Xnew[iPerm[i]] = Xold[i]
                 INTEGER N0,     // Maximum # of points per leaf node
                                 // or #levels if mPar == BBOX
                 unsigned Seed,  // Seed for RNG
                 PartMethod mPar = RAND // Partitioning method
                 );
  template<class Kernel, class Point, class PointArray>
  void BuildKernelMatrix(const Kernel &mKernel, // Kernel
                         const PointArray &X,   // Points (after permutation)
                         double lambda          // Regularization
                         );

  // Given a new set of points Y, let K be the kernel matrix between X
  // and Y; that is K = phi(X,Y). This routine computes z = K'*w for
  // any input vector w, without explicitly forming the matrix K.
  template<class Kernel, class Point, class PointArray>
  void MatVecImplicit(const PointArray &X,   // Points (after permutation)
                      const PointArray &Y,   // New points
                      const Kernel &mKernel, // Kernel
                      const DVector &w,      // The vector to multiply
                      DVector &z             // z = K'*w
                      );

  // This function is simiilar to the above one except that the
  // multiplication is with a matrix W.
  template<class Kernel, class Point, class PointArray>
  void MatMatImplicit(const PointArray &X,   // Points (after permutation)
                      const PointArray &Y,   // New points
                      const Kernel &mKernel, // Kernel
                      const DMatrix &W,      // The matrix to multiply
                      DMatrix &Z             // Z = K'*W
                      );

protected:

private:

  Node *Root;          // Tree root
  INTEGER MaxNumChild; // Maximum number of children for each node

  INTEGER N; // Matrix dimension
  INTEGER r; // Rank

  // Called by ConvertToDMatrix()
  void ConvertToDMatrixUpward(DMatrix &B, const Node *mNode,
                              MatState mState) const;

  // Called by MatVec()
  void MatVecUpward(const DVector &b, DVector &y,
                    Node *mNode, MatrixMode ModeA, MatState mState) const;

  // Called by MatMat()
  void MatMatUpward(const DMatrix &B, DMatrix &Y,
                    Node *mNode, MatrixMode ModeA, MatState mState) const;

  // Called by Invert()
  void InvertUpward(Node *mNodeTA, MatrixType TypeA) const;

  // Called by Det()
  LogDet DetUpward(Node *mNode, MatrixType TypeA, MatState mState) const;

  // Called by BuildTree()
  template<class Point, class PointArray>
  void BuildTreeDownward(Node *mNode, PointArray &X,
                         INTEGER *Perm, INTEGER N0, PartMethod mPar);

  // Called by BuildKernelMatrix()
  template<class Kernel, class Point, class PointArray>
  void BuildKernelMatrixDownward(Node *mNode, const PointArray &X,
                                 const Kernel &mKernel, double lambda);

  // Called by MatVecImplicit()
  template<class Kernel, class Point, class PointArray>
  void MatVecImplicitDownward(Node *mNode, const PointArray &X,
                              const PointArray &y, const Kernel &mKernel,
                              const DVector &w, double &z0);

  // Called by MatMatImplicit()
  template<class Kernel, class Point, class PointArray>
  void MatMatImplicitDownward(Node *mNode, const PointArray &X,
                              const PointArray &y, const Kernel &mKernel,
                              const DMatrix &W, DVector &z0);

};

#include "BMatrix.tpp"

#endif
