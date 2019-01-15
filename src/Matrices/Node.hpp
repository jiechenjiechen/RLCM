// The Node class implements the data structure for a tree node that
// is used to build the tree representation of a recursively low-rank
// compressed matrix in double precision.
//
// The implementation of this class is parallelized.

#ifndef _NODE_
#define _NODE_

#include "DMatrix.hpp"
#include "DPointArray.hpp"

class Node {

public:

  Node();
  Node(const Node &G);
  Node& operator= (const Node &G);
  ~Node();

  //-------------------- Utilities --------------------

  // Operations on this node
  void Init(void);
  void ReleaseAllMemory(void);
  void DeepCopy(const Node &G);
  void PrintNode(INTEGER MaxNumChild, INTEGER MaxN) const;

  // Get an estimation of the memory consumption of this node. The
  // estimation is based on a highly simplified model that counts only
  // A, Sigma, U, and W.
  INTEGER GetMemEstNode(void) const;

  // Get the memory consumption of A only.
  INTEGER GetMemEstNodeAonly(void) const;

  // Operations on the subtree rooted at this node. Typically called
  // by the root.
  void CopyTree(Node **mNode); // mNode points to a newly created node
  void DestroyTree(void);      // The current node is not destroyed
  void TreeReleaseAllMemory(void);
  void PrintTree(INTEGER MaxNumChild, INTEGER MaxN) const;

  // Get an estimation of the memory consumption of the subtree rooted
  // at this node. The estimation is based on a highly simplified
  // model that counts only A, Sigma, U, and W.
  INTEGER GetMemEstTree(void) const;

  // Get the memory consumption of all A's in the subtree rooted at
  // this node.
  INTEGER GetMemEstTreeAonly(void) const;

  //-------------------- Members --------------------

  // Basic tree structure
  Node *Parent;
  INTEGER NumChild;
  Node *LeftChild;
  Node *RightSibling;

  // Sizes
  INTEGER n;     // Number of indices
  INTEGER r;     // Rank
  INTEGER start; // Starting index

  // Matrix components
  DMatrix A;     // Size [n*n]
  DMatrix Sigma; // Size [r*r]
  DMatrix W;     // Size [r*r]
  DMatrix U;     // Size [n*r]

  // Augmented for unsymmetric matrix
  DMatrix Z;     // Size [r*r]
  DMatrix V;     // Size [n*r]

  // Augmented (for CMatrix::MatVec, CMatrix::MatVecImplicit, and
  //            CMatrix::BilinearImplicit)
  DVector c;     // Length [r]
  DVector d;     // Length [r]

  // Augmented (for CMatrix::MatMat and CMatrix::MatMatImplicit)
  DMatrix C;     // Size [r*m]
  DMatrix D;     // Size [r*m]

  // Augmented (for CMatrix::Invert, CMatrix::Sqrt, and
  //            CMatrix::BilinearImplicit)
  DMatrix E;     // Size [r*r]
  DMatrix Xi;    // Size [r*r]
  DMatrix Theta; // Size [r*r]
  LogDet mLogDet;

  // Augmented (for CMatrix::BuildKernelMatrix, CMatrix::MatVecImplicit, and
  //            CMatrix::MatMatImplicit)
  DPoint normal;   // Normal direction for partitioning hyperplane
  double offset;   // Inprod(Center,Normal)
  INTEGER *pivots; // Indices of the pivot points \ud{p}
  DMatrix FACT;    // Size [r*r]. Factored form of Phi_{\ud{p},\ud{p}}

  // Sometimes, it is more convenient, or even necessary, to
  // explicitly store the pivots. An example case is spatiotemporal
  // statistics (i.e., data is very low-dimensional), where the
  // partitioning uses a bounding box and the pivots, which lie on the
  // partitioning plane, are not a sampling of the original points.
  double *BBox;    // BBox[0] to BBox[Dim-1] give the lower bound,
                   // and BBox[Dim] to BBox[2*Dim-1] give the upper bound.
  INTEGER Dim;
  INTEGER PartDim; // Along which dimension is the box partitioned?
  DPointArray P;

protected:

private:

  // Called by PrintTree()
  void PrintTreeDownward(INTEGER MaxNumChild, INTEGER MaxN) const;

};

#endif
