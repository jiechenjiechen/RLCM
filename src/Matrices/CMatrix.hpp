// The CMatrix class implements a recursively low-rank compressed
// matrix in double precision. It supports both the symmetric mode and
// unsymmetric model. When the matrix is constructed from a kernel, it
// is symmetric positive definite.
//
// The implementation of this class is parallelized.

#ifndef _CMATRIX_
#define _CMATRIX_

#include "Node.hpp"
#ifdef HAS_METIS
#include "SMatrix.hpp"
#include <metis.h>
#endif

class CMatrix {

public:

  CMatrix();
  CMatrix(const CMatrix &G);
  CMatrix& operator= (const CMatrix &G);
  ~CMatrix();

  void Init(void);
  void ReleaseAllMemory(void);
  void DeepCopy(const CMatrix &G);

  //-------------------- Utilities --------------------
  //
  // A is the self matrix.

  // Is the matrix symmetric?
  bool IsSymmetric(void) const;

  // Set the symmetry property
  void SetSymmetric(bool IsSym_);

  // Get matrix size
  INTEGER GetN(void) const;

  // Get root
  Node* GetRoot(void) const;

  // Get an estimation of the memory consumption of the matrix. The
  // estimation is based on a highly simplified model that counts only
  // A, Sigma, U, and W.
  INTEGER GetMemEst(void) const;

  // Inspect tree structure
  void PrintTree(void) const;

  // Convert to dense matrix B = (DMatrix)A. This routine requires
  // O(N^2) memory. Use it wisely.
  void ConvertToDMatrix(DMatrix &B);

  //-------------------- Matrix computations --------------------
  //
  // A is the self matrix.

  // y = mode(A) * b
  //
  // NOTE: The argument mState is redundant; its appearance is only to
  // conform with the interface requirement of linear solvers. The
  // function will return error if one passes to the argument other
  // than UNFACT.
  void MatVec(const DVector &b, DVector &y, MatrixMode ModeA,
              MatState mState = UNFACT);

  // Y = mode(A) * B
  void MatMat(const DMatrix &B, DMatrix &Y, MatrixMode ModeA);

  // tA = inv(A)
  void Invert(CMatrix &tA);

  // Determinant. The user is responsible to ensure that Invert() has
  // been called before this function is called.
  LogDet Det(void) const;

  // A = GG'
  void Sqrt(CMatrix &G) const;

  //--------------- Build the matrix from points ---------------

  // Build a compressed kernel matrix from point array X. The
  // hierarchy tree is binary. The resulting matrix is SPD.
  //
  // The template classes Kernel, Point and PointArray are subject to
  // the same requirements as those for
  // DMatrix::BuildKernelMatrix. Additionally, the template class
  // PointArray must have the following methods:
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
                 INTEGER r_,     // Rank
                 INTEGER N0,     // Maximum # of points per leaf node,
                                 // or #levels if mPar == BBOX
                 double mDiagCorrect_,  // Diagonal correction
                 unsigned Seed,         // Seed for RNG
                 PartMethod mPar = RAND // Partitioning method
                 );
  template<class Kernel, class Point, class PointArray>
  void BuildKernelMatrix(const Kernel &mKernel, // Kernel
                         const PointArray &X,   // Points (after permutation)
                         double lambda          // Regularization
                         );

  // Given a new set of points Y, let K(X,Y) be the kernel matrix
  // between X and Y. This routine computes z = K'*w for any input
  // vector w, without explicitly forming the matrix K.
  template<class Kernel, class Point, class PointArray>
  void MatVecImplicit(const PointArray &X,   // Points (after permutation)
                      const PointArray &Y,   // New points
                      const Kernel &mKernel, // Kernel
                      const DVector &w,      // The vector to multiply
                      DVector &z             // z = K(X,Y)'*w
                      );

  // This function is simiilar to the above one except that the
  // multiplication is with a matrix W.
  template<class Kernel, class Point, class PointArray>
  void MatMatImplicit(const PointArray &X,   // Points (after permutation)
                      const PointArray &Y,   // New points
                      const Kernel &mKernel, // Kernel
                      const DMatrix &W,      // The matrix to multiply
                      DMatrix &Z             // Z = K(X,Y)'*W
                      );

  // Given a new set of points Y, for each i let K(X,Y[i]) be the
  // kernel vector between X and Y[i]. This routine computes for each
  // i, z[i] = K(Y[i],X) * tK * K(X,Y[i]), without explicitly forming
  // the vector K(X,Y[i]).
  template<class Kernel, class Point, class PointArray>
  void BilinearImplicit(const PointArray &X,   // Points (after permutation)
                        const PointArray &Y,   // New points
                        const Kernel &mKernel, // Kernel
                        CMatrix &tA,           // The matrix of bilinear form
                        DVector &z             // z[i] = K(Y[i],X)*tK*K(X,Y[i])
                        );

protected:

private:

  bool IsSym;          // Is the matrix symmetric

  Node *Root;          // Tree root
  INTEGER MaxNumChild; // Maximum number of children for each node

  INTEGER N; // Matrix dimension
  INTEGER r; // Rank

  // Augmented for kernel matrices in certain applications
  double mDiagCorrect;

  // Augmented for ConvertToDMatrix()
  DMatrix Ut;
  DMatrix Vt;

  // Called by ConvertToDMatrix()
  void ConvertToDMatrixUpward(DMatrix &B, const Node *mNode);

  // Called by MatVec()
  void MatVecInitAugmentedData(Node *mNode);
  void MatVecReleaseAugmentedData(Node *mNode);
  void MatVecUpward(const DVector &b, DVector &y,
                    Node *mNode, MatrixMode ModeA);
  void MatVecDownward(const DVector &b, DVector &y,
                      Node *mNode, MatrixMode ModeA);

  // Called by MatMat()
  void MatMatInitAugmentedData(Node *mNode, INTEGER NCOL);
  void MatMatReleaseAugmentedData(Node *mNode);
  void MatMatUpward(const DMatrix &B, DMatrix &Y,
                    Node *mNode, MatrixMode ModeA);
  void MatMatDownward(const DMatrix &B, DMatrix &Y,
                      Node *mNode, MatrixMode ModeA);

  // Called by Invert()
  void InvertInitAugmentedData(Node *mNodeTA);
  void InvertReleaseAugmentedData(Node *mNodeTA);
  void InvertUpward(Node *mNodeA, Node *mNodeTA);
  void InvertDownward(Node *mNodeTA);

  // Called by Det()
  LogDet DetUpward(Node *mNode) const;

  // Called by Sqrt()
  void SqrtInitAugmentedData(Node *mNodeG) const;
  void SqrtReleaseAugmentedData(Node *mNodeG) const;
  void SqrtUpward(Node *mNodeA, Node *mNodeG) const;
  void SqrtDownward(Node *mNodeG) const;
  void SqrtSolveRiccati(const DMatrix &Lambda, const DMatrix &Xi,
                        DMatrix &D) const;

  // Called by BuildTree()
  template<class Point, class PointArray>
  void BuildTreeDownward1(Node *mNode, PointArray &X,
                          INTEGER *Perm, INTEGER N0, PartMethod mPar);
  template<class Point, class PointArray>
  void BuildTreeDownward2(Node *mNode, const PointArray &X, PartMethod mPar);
  void SamplingPivotsRegularGrid(const double *BBox, INTEGER d, DPointArray &P);

  // Called by BuildKernelMatrix()
  template<class Kernel, class Point, class PointArray>
  void BuildKernelMatrixDownward(Node *mNode, const PointArray &X,
                                 const Kernel &mKernel, double lambda);

  // Called by MatVecImplicit()
  template<class Kernel, class Point, class PointArray>
  void MatVecImplicitInitAugmentedData(Node *mNode);
  template<class Kernel, class Point, class PointArray>
  void MatVecImplicitReleaseAugmentedData(Node *mNode);
  template<class Kernel, class Point, class PointArray>
  void MatVecImplicitUpward1(Node *mNode, const DVector &w);
  template<class Kernel, class Point, class PointArray>
  void MatVecImplicitUpward2(Node *mNode, const PointArray &X,
                             const PointArray &y, const Kernel &mKernel,
                             const DVector &W, double &z0);

  // Called by MatMatImplicit()
  template<class Kernel, class Point, class PointArray>
  void MatMatImplicitInitAugmentedData(Node *mNode, const DMatrix &W);
  template<class Kernel, class Point, class PointArray>
  void MatMatImplicitReleaseAugmentedData(Node *mNode);
  template<class Kernel, class Point, class PointArray>
  void MatMatImplicitUpward1(Node *mNode, const DMatrix &W);
  template<class Kernel, class Point, class PointArray>
  void MatMatImplicitUpward2(Node *mNode, const PointArray &X,
                             const PointArray &y, const Kernel &mKernel,
                             const DMatrix &W, DVector &z0);

  // Called by BilinearImplicit()
  template<class Kernel, class Point, class PointArray>
  void BilinearImplicitInitAugmentedData(Node *mNodeA, Node *mNodeTA);
  template<class Kernel, class Point, class PointArray>
  void BilinearImplicitReleaseAugmentedData(Node *mNodeA, Node *mNodeTA);
  template<class Kernel, class Point, class PointArray>
  void BilinearImplicitUpward1(Node *mNodeA, Node *mNodeTA);
  template<class Kernel, class Point, class PointArray>
  void BilinearImplicitUpward2(Node *mNodeA, Node *mNodeTA,
                               const PointArray &X, const PointArray &y,
                               const Kernel &mKernel, double &z0);

};

#include "CMatrix.tpp"

#endif
