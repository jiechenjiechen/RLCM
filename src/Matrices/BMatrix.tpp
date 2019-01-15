#ifndef _BMATRIX_TPP_
#define _BMATRIX_TPP_


//--------------------------------------------------------------------------
template<class Point, class PointArray>
void BMatrix::
BuildTree(PointArray &X,  // Points; will be permuted
          INTEGER *Perm,  // Xnew[i] = Xold[Perm[i]]
          INTEGER *iPerm, // Xnew[iPerm[i]] = Xold[i]
          INTEGER N0,     // Maximum # of points per leaf node
                          // or #levels if mPar == BBOX
          unsigned Seed,  // Seed for RNG
          PartMethod mPar // Partitioning method
          ) {

  // Seed the RNG
  srandom(Seed);

  // Destroy the old matrix
  Init();
  MaxNumChild = 2;
  N = X.GetN();
  r = N0;

  // Handle permutation
  if (Perm != NULL) {
    for (INTEGER i = 0; i < N; i++) {
      Perm[i] = i;
    }
  }

  // Create root
  Root = new Node;
  Root->r = r;
  Root->n = N;
  Root->start = 0;

  // Build a barebone tree. X will be permuted.
  BuildTreeDownward<Point, PointArray>(Root, X, Perm, N0, mPar);

  if (Root->NumChild == 0) {
    printf("BMatrix::BuildKernelMatrix. Error: Fail to partition the data set and thus cannot build the tree. Function call takes no effect.\n");
    Init();
    return;
  }

  // Compute the inverse permutation
  if (Perm != NULL && iPerm != NULL) {
    for (INTEGER i = 0; i < N; i++) {
      iPerm[Perm[i]] = i;
    }
  }

}


//--------------------------------------------------------------------------
template<class Point, class PointArray>
void BMatrix::
BuildTreeDownward(Node *mNode, PointArray &X,
                  INTEGER *Perm, INTEGER N0, PartMethod mPar) {

  // Partition the current point set
  INTEGER m1 = 0;
  switch (mPar) {
  case RAND:
    m1 = X.RandomBipartition(mNode->start, mNode->n, N0, Perm,
                             mNode->normal, mNode->offset);
    break;
  case PCA:
    m1 = X.PCABipartition(mNode->start, mNode->n, N0, Perm,
                          mNode->normal, mNode->offset);
    break;
  case BBOX: // Should never come here
    if (mNode == Root) {
      X.ComputeBBox(&(mNode->BBox), mNode->Dim);
    }
    m1 = X.BBoxBipartition(mNode->start, mNode->n, N0, Perm,
                           mNode->normal, mNode->offset,
                           mNode->BBox, mNode->PartDim);
    break;
  }

  // If partitioning yields clusters of size less than N0, quit.
  // This stopping condition applies to only mPar == RAND or PCA
  if ((mPar == RAND || mPar == PCA) && (m1 < N0)) {
    return;
  }

  // Otherwise, create children
  mNode->NumChild = 2;

  Node *NewNodeL = new Node;
  NewNodeL->Parent = mNode;
  mNode->LeftChild = NewNodeL;
  NewNodeL->r = r;
  NewNodeL->start = mNode->start;
  NewNodeL->n = m1;

  if (mPar == BBOX) { // Should never come here
    NewNodeL->Dim = X.GetD();
    New_1D_Array<double, INTEGER>(&(NewNodeL->BBox), NewNodeL->Dim*2);
    memcpy(NewNodeL->BBox, mNode->BBox, NewNodeL->Dim*2*sizeof(double));
    NewNodeL->BBox[mNode->PartDim + NewNodeL->Dim] = mNode->offset;
  }

  Node *NewNodeR = new Node;
  NewNodeR->Parent = mNode;
  NewNodeL->RightSibling = NewNodeR;
  NewNodeR->r = r;
  NewNodeR->start = mNode->start + m1;
  NewNodeR->n = mNode->n - m1;

  if (mPar == BBOX) { // Should never come here
    NewNodeR->Dim = X.GetD();
    New_1D_Array<double, INTEGER>(&(NewNodeR->BBox), NewNodeR->Dim*2);
    memcpy(NewNodeR->BBox, mNode->BBox, NewNodeR->Dim*2*sizeof(double));
    NewNodeR->BBox[mNode->PartDim] = mNode->offset;
  }


  // If no more levels are needed, quit.
  // This stopping condition applies to only mPar == BBOX
  if ((mPar == BBOX) && (N0 <= 1)) {
    return;
  }

  // Recurse on children
  switch (mPar) {
  case RAND:
  case PCA:
    BuildTreeDownward<Point, PointArray>(NewNodeL, X, Perm, N0, mPar);
    BuildTreeDownward<Point, PointArray>(NewNodeR, X, Perm, N0, mPar);
    break;
  case BBOX:
    BuildTreeDownward<Point, PointArray>(NewNodeL, X, Perm, N0-1, mPar);
    BuildTreeDownward<Point, PointArray>(NewNodeR, X, Perm, N0-1, mPar);
    break;
  }

}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void BMatrix::
BuildKernelMatrix(const Kernel &mKernel, // Kernel
                  const PointArray &X,   // Points (after permutation)
                  double lambda          // Regularization
                  ) {

  // Instantiate matrix components
  BuildKernelMatrixDownward<Kernel, Point, PointArray>
    (Root, X, mKernel, lambda);

}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void BMatrix::
BuildKernelMatrixDownward(Node *mNode, const PointArray &X,
                          const Kernel &mKernel, double lambda) {

  INTEGER i = 0;
  Node *childi = NULL;

  if (mNode->NumChild == 0) {

    // A
    PointArray Y;
    X.GetSubset(mNode->start, mNode->n, Y);
    mNode->A.BuildKernelMatrix<Kernel, Point, PointArray>(mKernel, Y);
    mNode->A.AddDiagonal(lambda);

    return;

  }
  else {

    // Recurse on children
    childi = mNode->LeftChild;
    for (i = 0; i < mNode->NumChild; i++, childi = childi->RightSibling) {
      BuildKernelMatrixDownward<Kernel, Point, PointArray>
        (childi, X, mKernel, lambda);
    }

  }

}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void BMatrix::
MatVecImplicit(const PointArray &X,   // Points (after permutation)
               const PointArray &Y,   // New points
               const Kernel &mKernel, // Kernel
               const DVector &w,      // The vector to multiply
               DVector &z             // z = B'*w
               ) {

  if (Root == NULL) {
    printf("BMatrix::MatVecImplicit. Error: Matrix is empty. Function call takes no effect.\n");
    return;
  }
  if (w.GetN() != X.GetN()) {
    printf("BMatrix::MatVecImplicit. Error: Length of w does not match number of points in X. Function call takes no effect.\n");
  }

  INTEGER n = Y.GetN();
  z.Init(n);

  // Do mat-vec for one point at a time
  for (INTEGER i = 0; i < n; i++) {
    PointArray y;
    Y.GetSubset(i, 1, y);
    double z0 = 0.0;
    MatVecImplicitDownward<Kernel, Point, PointArray>
      (Root, X, y, mKernel, w, z0);
    z.SetEntry(i, z0);
  }

}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void BMatrix::
MatVecImplicitDownward(Node *mNode, const PointArray &X,
                       const PointArray &y, const Kernel &mKernel,
                       const DVector &w, double &z0) {

  // y must be a singleton set

  Node *child = NULL;

  if (mNode->NumChild == 0) {

    // ph is the kernel matrix between the old points (corresponding
    // subset of X) and the new point (a singleton set y)
    PointArray P;
    DMatrix PH;
    DVector ph;
    X.GetSubset(mNode->start, mNode->n, P);
    PH.BuildKernelMatrix<Kernel, Point, PointArray>(mKernel, P, y);
    PH.GetColumn(0, ph);

    // Multiply ph and the corresponding block of w. The result is z0
    DVector mw;
    w.GetBlock(mNode->start, mNode->n, mw);
    z0 = ph.InProd(mw);

  }
  else {

    // Which child does y land on? The computation here is applicable
    // to only a binary tree constructed by using hyperplane
    // partitioning.
    Point y0;
    y.GetPoint(0, y0);
    double iprod = y0.InProd(mNode->normal);
    if (iprod < mNode->offset) {
      child = mNode->LeftChild;
    }
    else {
      child = mNode->LeftChild->RightSibling;
    }

    // Recurse on that child only
    MatVecImplicitDownward<Kernel, Point, PointArray>
      (child, X, y, mKernel, w, z0);

  }

}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void BMatrix::
MatMatImplicit(const PointArray &X,   // Points (after permutation)
               const PointArray &Y,   // New points
               const Kernel &mKernel, // Kernel
               const DMatrix &W,      // The matrix to multiply
               DMatrix &Z             // Z = B'*W
               ) {

  if (Root == NULL) {
    printf("BMatrix::MatMatImplicit. Error: Matrix is empty. Function call takes no effect.\n");
    return;
  }
  if (W.GetM() != X.GetN()) {
    printf("BMatrix::MatMatImplicit. Error: Size of W does not match the number of points in X. Function call takes no effect.\n");
  }

  INTEGER n = Y.GetN();
  INTEGER m = W.GetN();
  Z.Init(n, m);

  // Do mat-vec for one point at a time
  for (INTEGER i = 0; i < n; i++) {
    PointArray y;
    Y.GetSubset(i, 1, y);
    DVector z0(m);
    MatMatImplicitDownward<Kernel, Point, PointArray>
      (Root, X, y, mKernel, W, z0);
    Z.SetRow(i, z0);
  }

}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void BMatrix::
MatMatImplicitDownward(Node *mNode, const PointArray &X,
                       const PointArray &y, const Kernel &mKernel,
                       const DMatrix &W, DVector &z0) {

  // y must be a singleton set

  Node *child = NULL;
  INTEGER m = W.GetN();

  if (mNode->NumChild == 0) {

    // ph is the kernel matrix between the old points (corresponding
    // subset of X) and the new point (a singleton set y)
    PointArray P;
    DMatrix PH;
    DVector ph;
    X.GetSubset(mNode->start, mNode->n, P);
    PH.BuildKernelMatrix<Kernel, Point, PointArray>(mKernel, P, y);
    PH.GetColumn(0, ph);

    // Multiply ph and the corresponding block of W. The result is z0
    DMatrix mW;
    W.GetBlock(mNode->start, mNode->n, 0, m, mW);
    mW.MatVec(ph, z0, TRANSPOSE);

  }
  else {

    // Which child does y land on? The computation here is applicable
    // to only a binary tree constructed by using hyperplane
    // partitioning.
    Point y0;
    y.GetPoint(0, y0);
    double iprod = y0.InProd(mNode->normal);
    if (iprod < mNode->offset) {
      child = mNode->LeftChild;
    }
    else {
      child = mNode->LeftChild->RightSibling;
    }

    // Recurse on that child only
    MatMatImplicitDownward<Kernel, Point, PointArray>
      (child, X, y, mKernel, W, z0);

  }

}


#endif

