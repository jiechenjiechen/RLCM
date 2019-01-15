#ifndef _CMATRIX_TPP_
#define _CMATRIX_TPP_

//--------------------------------------------------------------------------
template<class Point, class PointArray>
void CMatrix::
BuildTree(PointArray &X,  // Points; will be permuted
          INTEGER *Perm,  // Xnew[i] = Xold[Perm[i]]
          INTEGER *iPerm, // Xnew[iPerm[i]] = Xold[i]
          INTEGER r_,     // Rank
          INTEGER N0,     // Maximum # of points per leaf node
                          // or #levels if mPar == BBOX
          double mDiagCorrect_,  // Diagonal correction
          unsigned Seed,         // Seed for RNG
          PartMethod mPar        // Partitioning method
          ) {

  if ((mPar == RAND || mPar == PCA) && (N0 < r_)) {
    printf("CMatrix::BuildTree. Error: N0 must >= r. Function call takes no effect.\n");
    return;
  }

  // Seed the RNG
  srandom(Seed);

  // Destroy the old matrix
  Init();
  IsSym = true; // Always symmetric
  MaxNumChild = 2;
  N = X.GetN();
  r = r_;
  mDiagCorrect = mDiagCorrect_;

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
  BuildTreeDownward1<Point, PointArray>(Root, X, Perm, N0, mPar);

  if (Root->NumChild == 0) {
    printf("CMatrix::BuildTree. Error: Fail to partition the data set and thus cannot build the tree. Function call takes no effect.\n");
    Init();
    return;
  }

  // Find pivots for each nonleaf node
  BuildTreeDownward2<Point, PointArray>(Root, X, mPar);

  // Compute the inverse permutation
  if (Perm != NULL && iPerm != NULL) {
    for (INTEGER i = 0; i < N; i++) {
      iPerm[Perm[i]] = i;
    }
  }

}


//--------------------------------------------------------------------------
template<class Point, class PointArray>
void CMatrix::
BuildTreeDownward1(Node *mNode, PointArray &X,
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
  case BBOX:
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

  if (mPar == BBOX) {
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

  if (mPar == BBOX) {
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
    BuildTreeDownward1<Point, PointArray>(NewNodeL, X, Perm, N0, mPar);
    BuildTreeDownward1<Point, PointArray>(NewNodeR, X, Perm, N0, mPar);
    break;
  case BBOX:
    BuildTreeDownward1<Point, PointArray>(NewNodeL, X, Perm, N0-1, mPar);
    BuildTreeDownward1<Point, PointArray>(NewNodeR, X, Perm, N0-1, mPar);
    break;
  }

}


//--------------------------------------------------------------------------
template<class Point, class PointArray>
void CMatrix::
BuildTreeDownward2(Node *mNode, const PointArray &X, PartMethod mPar) {

  INTEGER i = 0;
  Node *childi = NULL;

  // Leaf
  if (mNode->NumChild == 0) {
    return;
  }

  // Sample pivot points
  switch (mPar) {
  case RAND:
  case PCA:
    { // Sample from X
      New_1D_Array<INTEGER, INTEGER>(&mNode->pivots, r);
      RandPerm(mNode->n, r, mNode->pivots);
      for (i = 0; i < r; i++) {
        mNode->pivots[i] += mNode->start;
      }
    }
    break;

  case BBOX:
    { // Almost a regular grid
      SamplingPivotsRegularGrid(mNode->BBox, X.GetD(), mNode->P);
    }
    break;
  }

  // Recurse on children
  childi = mNode->LeftChild;
  for (i = 0; i < mNode->NumChild; i++, childi = childi->RightSibling) {
    BuildTreeDownward2<Point, PointArray>(childi, X, mPar);
  }

}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void CMatrix::
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
void CMatrix::
BuildKernelMatrixDownward(Node *mNode, const PointArray &X,
                          const Kernel &mKernel, double lambda) {

  INTEGER i = 0;
  Node *childi = NULL;
  Node *parent = mNode->Parent;

  // Leaf
  if (mNode->NumChild == 0) {

    // A
    PointArray Y;
    X.GetSubset(mNode->start, mNode->n, Y);
    mNode->A.BuildKernelMatrix<Kernel, Point, PointArray>(mKernel, Y);
    mNode->A.AddDiagonal(lambda);

    if (mNode == Root) {
      return;
    }

    // U
    PointArray P;
    if (parent->pivots != NULL) { // pivots are stored in parent->pivots
      X.GetSubset(parent->pivots, r, P);
    }
    else { // pivots are stored in parent->P
      P = parent->P;
    }
    DMatrix UU;
    UU.BuildKernelMatrix<Kernel, Point, PointArray>
      (mKernel, P, Y, mDiagCorrect*mKernel.GetS());
    parent->FACT.DPOTRS(UU, mNode->U, UPPER);
    mNode->U.Transpose();

    return;

  }

  // From here on, mNode is not a leaf and thus has children

  // Pivot points P
  PointArray P;
  if (mNode->pivots != NULL) { // pivots are stored in mNode->pivots
    X.GetSubset(mNode->pivots, r, P);
  }
  else { // pivots are stored in mNode->P
    P = mNode->P;
  }
  // Compute the kernel matrix PH of pivots
  DMatrix PH;
  PH.BuildKernelMatrix<Kernel, Point, PointArray>
    (mKernel, P, mDiagCorrect*mKernel.GetS());
  // Compute Cholesky factorization of PH
  mNode->FACT = PH;
  mNode->FACT.DPOTRF(UPPER);

  // Sigma
  mNode->Sigma = PH;

  // W
  if (mNode != Root) {
    PointArray PP;
    if (parent->pivots != NULL) { // pivots are stored in parent->pivots
      X.GetSubset(parent->pivots, r, PP);
    }
    else { // pivots are stored in parent->P
      PP = parent->P;
    }
    DMatrix HH;
    HH.BuildKernelMatrix<Kernel, Point, PointArray>
      (mKernel, PP, P, mDiagCorrect*mKernel.GetS());
    parent->FACT.DPOTRS(HH, mNode->W, UPPER);
    mNode->W.Transpose();
  }

  // Recurse on children
  childi = mNode->LeftChild;
  for (i = 0; i < mNode->NumChild; i++, childi = childi->RightSibling) {
    BuildKernelMatrixDownward<Kernel, Point, PointArray>
      (childi, X, mKernel, lambda);
  }

}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void CMatrix::
MatVecImplicit(const PointArray &X,   // Points (after permutation)
               const PointArray &Y,   // New points
               const Kernel &mKernel, // Kernel
               const DVector &w,      // The vector to multiply
               DVector &z             // z = K(X,Y)'*w
               ) {

  if (Root == NULL) {
    printf("CMatrix::MatVecImplicit. Error: Matrix is empty. Function call takes no effect.\n");
    return;
  }
  if (w.GetN() != X.GetN()) {
    printf("CMatrix::MatVecImplicit. Error: Length of w does not match the number of points in X. Function call takes no effect.\n");
  }

  INTEGER n = Y.GetN();
  z.Init(n);

  MatVecImplicitInitAugmentedData<Kernel, Point, PointArray>(Root);

  MatVecImplicitUpward1<Kernel, Point, PointArray>(Root, w);

  // Do mat-vec for one point at a time
  for (INTEGER i = 0; i < n; i++) {
    PointArray y;
    Y.GetSubset(i, 1, y);
    double z0 = 0.0;
    MatVecImplicitUpward2<Kernel, Point, PointArray>
      (Root, X, y, mKernel, w, z0);
    z.SetEntry(i, z0);
  }

  MatVecImplicitReleaseAugmentedData<Kernel, Point, PointArray>(Root);

}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void CMatrix::
MatVecImplicitInitAugmentedData(Node *mNode) {

  INTEGER i = 0;
  Node *child = NULL;

  if (mNode != Root) {
    mNode->c.Init(r);
    mNode->d.Init(r);
  }

  child = mNode->LeftChild;
  for (i = 0; i < mNode->NumChild; i++, child = child->RightSibling) {
    MatVecImplicitInitAugmentedData<Kernel, Point, PointArray>(child);
  }

}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void CMatrix::
MatVecImplicitReleaseAugmentedData(Node *mNode) {

  INTEGER i = 0;
  Node *child = NULL;

  if (mNode != Root) {
    mNode->c.ReleaseAllMemory();
    mNode->d.ReleaseAllMemory();
  }

  child = mNode->LeftChild;
  for (i = 0; i < mNode->NumChild; i++, child = child->RightSibling) {
    MatVecImplicitReleaseAugmentedData<Kernel, Point, PointArray>(child);
  }

}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void CMatrix::
MatVecImplicitUpward1(Node *mNode, const DVector &w) {

  // In the first upward phase, d is used to store the intermediate
  // mat-vec results computed bottom up, and c is used to store the
  // result from cross-multiply.

  INTEGER j = 0, k = 0;
  Node *child = NULL, *sibling = NULL;
  Node *parent = mNode->Parent;

  if (mNode->NumChild == 0) {

    // d_i = U_i' * w_i
    DVector mw;
    w.GetBlock(mNode->start, mNode->n, mw);
    mNode->U.MatVec(mw, mNode->d, TRANSPOSE);

  }
  else {

    DVector md;
    // d_i = W_i' * (sum_{j \in Ch(i)} d_j), part 0: init
    if (mNode != Root) {
      md.Init(r);
    }

    child = mNode->LeftChild;
    for (j = 0; j < mNode->NumChild; j++, child = child->RightSibling) {

      // Recurse on children
      MatVecImplicitUpward1<Kernel, Point, PointArray>(child, w);
      if (mNode == Root) { continue; }

      // d_i = W_i' * (sum_{j \in Ch(i)} d_j), part 1: summing d_j
      md.Add(child->d);

    }

    // d_i = W_i' * (sum_{j \in Ch(i)} d_j), part 2: left-mult W_i'
    if (mNode != Root) {
      mNode->W.MatVec(md, mNode->d, TRANSPOSE);
    }

  }

  if (mNode == Root) { return; }

  // c_k = Sigma_p' * d_i. Note that Sigma_p is symmetric
  DVector ck;
  parent->Sigma.MatVec(mNode->d, ck, TRANSPOSE);

  sibling = parent->LeftChild;
  for (k = 0; k < parent->NumChild; k++, sibling = sibling->RightSibling) {
    // get k (k must be unique in a binary tree)
    if (sibling == mNode) { continue; }
    sibling->c = ck;
  }

}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void CMatrix::
MatVecImplicitUpward2(Node *mNode, const PointArray &X,
                      const PointArray &y, const Kernel &mKernel,
                      const DVector &w, double &z0) {

  // In the second upward phase, the storage of d is recycled and is
  // used to store the new intermediate mat-vec results computed
  // bottom up.

  // y must be a singleton set

  Node *parent = mNode->Parent;
  Node *child = NULL;

  if (mNode->NumChild == 0) {

    // Initial d_i (at leaf)
    PointArray P;
    DMatrix PH;
    DVector ph;
    if (mNode != Root) {
      if (parent->pivots != NULL) { // pivots are stored in parent->pivots
        X.GetSubset(parent->pivots, r, P);
      }
      else { // pivots are stored in parent->P
        P = parent->P;
      }
      PH.BuildKernelMatrix<Kernel, Point, PointArray>
        (mKernel, P, y, mDiagCorrect*mKernel.GetS());
      PH.GetColumn(0, ph);
      parent->FACT.DPOTRS(ph, mNode->d, UPPER);
    }

    // ph is the kernel matrix between the old points (corresponding
    // subset of X) and the new point (a singleton set y)
    X.GetSubset(mNode->start, mNode->n, P);
    PH.BuildKernelMatrix<Kernel, Point, PointArray>(mKernel, P, y);
    PH.GetColumn(0, ph);

    // Multiply ph and the corresponding block of w. The result is
    // the starting z0. It will accumulate later.
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
    MatVecImplicitUpward2<Kernel, Point, PointArray>
      (child, X, y, mKernel, w, z0);

    // d_i = W_i' * d_j
    if (mNode != Root) {
      mNode->W.MatVec(child->d, mNode->d, TRANSPOSE);
    }

  }

  // multiply c_i with d_i and accumulate to z0
  if (mNode != Root) {
    z0 += mNode->c.InProd(mNode->d);
  }

}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void CMatrix::
MatMatImplicit(const PointArray &X,   // Points (after permutation)
               const PointArray &Y,   // New points
               const Kernel &mKernel, // Kernel
               const DMatrix &W,      // The matrix to multiply
               DMatrix &Z             // Z = K(X,Y)'*W
               ) {

  if (Root == NULL) {
    printf("CMatrix::MatMatImplicit. Error: Matrix is empty. Function call takes no effect.\n");
    return;
  }
  if (W.GetM() != X.GetN()) {
    printf("CMatrix::MatMatImplicit. Error: Size of W does not match the number of points in X. Function call takes no effect.\n");
    return;
  }

  INTEGER n = Y.GetN();
  INTEGER m = W.GetN();
  Z.Init(n,m);

  MatMatImplicitInitAugmentedData<Kernel, Point, PointArray>(Root, W);

  MatMatImplicitUpward1<Kernel, Point, PointArray>(Root, W);

  // Do mat-mat for one point at a time
  for (INTEGER i = 0; i < n; i++) {
    PointArray y;
    Y.GetSubset(i, 1, y);
    DVector z0(m);
    MatMatImplicitUpward2<Kernel, Point, PointArray>
      (Root, X, y, mKernel, W, z0);
    Z.SetRow(i, z0);
  }

  MatMatImplicitReleaseAugmentedData<Kernel, Point, PointArray>(Root);

}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void CMatrix::
MatMatImplicitInitAugmentedData(Node *mNode, const DMatrix &W) {

  INTEGER i = 0;
  Node *child = NULL;
  INTEGER m = W.GetN();

  if (mNode != Root) {
    mNode->C.Init(r, m);
    mNode->D.Init(r, m);
    mNode->d.Init(r);
  }

  child = mNode->LeftChild;
  for (i = 0; i < mNode->NumChild; i++, child = child->RightSibling) {
    MatMatImplicitInitAugmentedData<Kernel, Point, PointArray>(child, W);
  }

}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void CMatrix::
MatMatImplicitReleaseAugmentedData(Node *mNode) {

  INTEGER i = 0;
  Node *child = NULL;

  if (mNode != Root) {
    mNode->C.ReleaseAllMemory();
    mNode->D.ReleaseAllMemory();
    mNode->d.ReleaseAllMemory();
  }

  child = mNode->LeftChild;
  for (i = 0; i < mNode->NumChild; i++, child = child->RightSibling) {
    MatMatImplicitReleaseAugmentedData<Kernel, Point, PointArray>(child);
  }

}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void CMatrix::
MatMatImplicitUpward1(Node *mNode, const DMatrix &W) {

  // In the first upward phase, D is used to store the intermediate
  // mat-mat results computed bottom up, and C is used to store the
  // result from cross-multiply.

  INTEGER j = 0, k = 0;
  Node *child = NULL, *sibling = NULL;
  Node *parent = mNode->Parent;
  INTEGER m = W.GetN();

  if (mNode->NumChild == 0) {

    // D_i = U_i' * W_i
    DMatrix mW;
    W.GetBlock(mNode->start, mNode->n, 0, m, mW);
    mNode->U.MatMat(mW, mNode->D, TRANSPOSE, NORMAL);

  }
  else {

    DMatrix mD;
    // D_i = W_i' * (sum_{j \in Ch(i)} D_j), part 0: init
    if (mNode != Root) {
      mD.Init(r, m);
    }

    child = mNode->LeftChild;
    for (j = 0; j < mNode->NumChild; j++, child = child->RightSibling) {

      // Recurse on children
      MatMatImplicitUpward1<Kernel, Point, PointArray>(child, W);
      if (mNode == Root) { continue; }

      // D_i = W_i' * (sum_{j \in Ch(i)} D_j), part 1: summing D_j
      mD.Add(child->D);

    }

    // D_i = W_i' * (sum_{j \in Ch(i)} D_j), part 2: left-mult W_i'
    if (mNode != Root) {
      mNode->W.MatMat(mD, mNode->D, TRANSPOSE, NORMAL);
    }

  }

  if (mNode == Root) { return; }

  // C_k = Sigma_p' * D_i. Note that Sigma_p is symmetric
  DMatrix Ck;
  parent->Sigma.MatMat(mNode->D, Ck, TRANSPOSE, NORMAL);

  sibling = parent->LeftChild;
  for (k = 0; k < parent->NumChild; k++, sibling = sibling->RightSibling) {
    // get k (k must be unique in a binary tree)
    if (sibling == mNode) { continue; }
    sibling->C = Ck;
  }

}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void CMatrix::
MatMatImplicitUpward2(Node *mNode, const PointArray &X,
                      const PointArray &y, const Kernel &mKernel,
                      const DMatrix &W, DVector &z0) {

  // In the second upward phase, the storage of d is used to store the
  // new intermediate mat-mat results computed bottom up.

  // y must be a singleton set

  Node *parent = mNode->Parent;
  Node *child = NULL;
  INTEGER m = W.GetN();

  if (mNode->NumChild == 0) {

    // Initial d_i (at leaf)
    PointArray P;
    DMatrix PH;
    DVector ph;
    if (mNode != Root) {
      if (parent->pivots != NULL) { // pivots are stored in parent->pivots
        X.GetSubset(parent->pivots, r, P);
      }
      else { // pivots are stored in parent->P
        P = parent->P;
      }
      PH.BuildKernelMatrix<Kernel, Point, PointArray>
        (mKernel, P, y, mDiagCorrect*mKernel.GetS());
      PH.GetColumn(0, ph);
      parent->FACT.DPOTRS(ph, mNode->d, UPPER);
    }

    // ph is the kernel matrix between the old points (corresponding
    // subset of X) and the new point (a singleton set y)
    X.GetSubset(mNode->start, mNode->n, P);
    PH.BuildKernelMatrix<Kernel, Point, PointArray>(mKernel, P, y);
    PH.GetColumn(0, ph);

    // Multiply ph and the corresponding block of W. The result is
    // the starting z0. It will accumulate later.
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
    MatMatImplicitUpward2<Kernel, Point, PointArray>
      (child, X, y, mKernel, W, z0);

    // d_i = W_i' * d_j
    if (mNode != Root) {
      mNode->W.MatVec(child->d, mNode->d, TRANSPOSE);
    }

  }

  // multiply C_i with d_i and accumulate to z0
  if (mNode != Root) {
    mNode->C.DGEMV(mNode->d, z0, 1.0, 1.0, TRANSPOSE);
  }

}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void CMatrix::
BilinearImplicit(const PointArray &X,   // Points (after permutation)
                 const PointArray &Y,   // New points
                 const Kernel &mKernel, // Kernel
                 CMatrix &tA,           // The matrix of bilinear form
                 DVector &z             // z[i] = K(Y[i],X)*tK*K(X,Y[i])
                 ) {

  if (Root == NULL) {
    printf("CMatrix::BilinearImplicit. Error: Matrix is empty. Function call takes no effect.\n");
    return;
  }

  INTEGER n = Y.GetN();
  z.Init(n);
  Node *RootTA = tA.GetRoot();

  BilinearImplicitInitAugmentedData<Kernel, Point, PointArray>(Root, RootTA);

  BilinearImplicitUpward1<Kernel, Point, PointArray>(Root, RootTA);

  // One point at a time
  for (INTEGER i = 0; i < n; i++) {
    PointArray y;
    Y.GetSubset(i, 1, y);
    double z0 = 0.0;
    BilinearImplicitUpward2<Kernel, Point, PointArray>
      (Root, RootTA, X, y, mKernel, z0);
    z.SetEntry(i, z0);
  }

  BilinearImplicitReleaseAugmentedData<Kernel, Point, PointArray>(Root, RootTA);

}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void CMatrix::
BilinearImplicitInitAugmentedData(Node *mNodeA, Node *mNodeTA) {

  INTEGER i = 0;
  Node *childA = NULL, *childTA = NULL;

  if (mNodeA != Root) {
    mNodeA->c.Init(r);
    mNodeA->d.Init(r);
    mNodeA->Theta.Init(r);
    mNodeA->Xi.Init(r);
    mNodeTA->Theta.Init(r);
    mNodeTA->Xi.Init(r);
  }

  childA = mNodeA->LeftChild;
  childTA = mNodeTA->LeftChild;
  for (i = 0; i < mNodeA->NumChild; i++, childA = childA->RightSibling,
         childTA = childTA->RightSibling) {
    BilinearImplicitInitAugmentedData<Kernel, Point, PointArray>
      (childA, childTA);
  }

}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void CMatrix::
BilinearImplicitReleaseAugmentedData(Node *mNodeA, Node *mNodeTA) {

  INTEGER i = 0;
  Node *childA = NULL, *childTA = NULL;

  if (mNodeA != Root) {
    mNodeA->c.ReleaseAllMemory();
    mNodeA->d.ReleaseAllMemory();
    mNodeA->Theta.ReleaseAllMemory();
    mNodeA->Xi.ReleaseAllMemory();
    mNodeTA->Theta.ReleaseAllMemory();
    mNodeTA->Xi.ReleaseAllMemory();
  }

  childA = mNodeA->LeftChild;
  childTA = mNodeTA->LeftChild;
  for (i = 0; i < mNodeA->NumChild; i++, childA = childA->RightSibling,
         childTA = childTA->RightSibling) {
    BilinearImplicitReleaseAugmentedData<Kernel, Point, PointArray>
      (childA, childTA);
  }

}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void CMatrix::
BilinearImplicitUpward1(Node *mNodeA, Node *mNodeTA) {

  INTEGER j = 0, k = 0;
  Node *childA = NULL, *childAj = NULL, *childAk = NULL, *childTA = NULL;
  Node *parentA = mNodeA->Parent;

  if (mNodeA->NumChild == 0 && parentA != NULL) {

    // Theta_i = tU_i' * U_i
    mNodeTA->U.MatMat(mNodeA->U, mNodeA->Theta, TRANSPOSE, NORMAL);

    // tTheta_i = Theta_i * Sigma_p
    mNodeA->Theta.MatMat(parentA->Sigma, mNodeTA->Theta, NORMAL, NORMAL);

    // Xi_i = U_i' * tA_{ii} * U_i
    DMatrix mUT;
    mNodeA->U.MatMat(mNodeTA->A, mUT, TRANSPOSE, NORMAL);
    mUT.MatMat(mNodeA->U, mNodeA->Xi, NORMAL, NORMAL);

    // tXi_i = Sigma_p' * Xi_i * Sigma_p
    DMatrix mSX;
    parentA->Sigma.MatMat(mNodeA->Xi, mSX, TRANSPOSE, NORMAL);
    mSX.MatMat(parentA->Sigma, mNodeTA->Xi, NORMAL, NORMAL);

    return;

  }

  childA = mNodeA->LeftChild;
  childTA = mNodeTA->LeftChild;
  for (j = 0; j < mNodeA->NumChild; j++, childA = childA->RightSibling,
         childTA = childTA->RightSibling) {

    // Recurse on children
    BilinearImplicitUpward1<Kernel, Point, PointArray>(childA, childTA);

  }

  if (parentA == NULL) {
    return;
  }

  // Theta_i = tW_i' * (sum_{j \in Ch(i)} Theta_j) * W_i
  DMatrix mST;
  mST.Init(r);
  childA = mNodeA->LeftChild;
  for (j = 0; j < mNodeA->NumChild; j++, childA = childA->RightSibling) {
    mST.Add(childA->Theta);
  }
  DMatrix mWST;
  mNodeTA->W.MatMat(mST, mWST, TRANSPOSE, NORMAL);
  mWST.MatMat(mNodeA->W, mNodeA->Theta, NORMAL, NORMAL);

  // tTheta_i = Theta_i * Sigma_p
  mNodeA->Theta.MatMat(parentA->Sigma, mNodeTA->Theta, NORMAL, NORMAL);

  // Xi_i = W_i' * [ sum_{j \in Ch(i)} Xi_j
  //    + sum_{j,k \in Chi(i), j \ne k} Theta_j' * tSigma_i * Theta_k ] * W_i
  DMatrix mXi2;
  mXi2.Init(r);
  childAj = mNodeA->LeftChild;
  for (j = 0; j < mNodeA->NumChild; j++, childAj = childAj->RightSibling) {
    childAk = mNodeA->LeftChild;
    for (k = 0; k < mNodeA->NumChild; k++, childAk = childAk->RightSibling) {

      if (j == k) {
        mXi2.Add(childAj->Xi);
      }
      else {
        DMatrix mTS, mTST;
        childAj->Theta.MatMat(mNodeTA->Sigma, mTS, TRANSPOSE, NORMAL);
        mTS.MatMat(childAk->Theta, mTST, NORMAL, NORMAL);
        mXi2.Add(mTST);
      }

    }
  }
  DMatrix mWX;
  mNodeA->W.MatMat(mXi2, mWX, TRANSPOSE, NORMAL);
  mWX.MatMat(mNodeA->W, mNodeA->Xi, NORMAL, NORMAL);

  // tXi_i = Sigma_p' * Xi_i * Sigma_p
  DMatrix mSX;
  parentA->Sigma.MatMat(mNodeA->Xi, mSX, TRANSPOSE, NORMAL);
  mSX.MatMat(parentA->Sigma, mNodeTA->Xi, NORMAL, NORMAL);

}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void CMatrix::
BilinearImplicitUpward2(Node *mNodeA, Node *mNodeTA,
                        const PointArray &X, const PointArray &y,
                        const Kernel &mKernel, double &z0) {

  // y must be a singleton set

  INTEGER k = 0;
  Node *parentA = mNodeA->Parent, *parentTA = mNodeTA->Parent;
  Node *childA = NULL, *childTA = NULL, *siblingA = NULL, *siblingTA = NULL;

  if (mNodeA->NumChild == 0) {

    // d_i = K(landmark,lanmark) \ k(lankmark,x)
    //                                   ||
    //                                   ph
    PointArray P;
    DMatrix PH;
    DVector ph;
    if (mNodeA != Root) {
      if (parentA->pivots != NULL) { // pivots are stored in parentA->pivots
        X.GetSubset(parentA->pivots, r, P);
      }
      else { // pivots are stored in parentA->P
        P = parentA->P;
      }
      PH.BuildKernelMatrix<Kernel, Point, PointArray>
        (mKernel, P, y, mDiagCorrect*mKernel.GetS());
      PH.GetColumn(0, ph);
      parentA->FACT.DPOTRS(ph, mNodeA->d, UPPER);
    }

    // Reuse PH and ph: ph is the kernel matrix between the old points
    // (corresponding subset of X) and the new point (a singleton set
    // y)
    X.GetSubset(mNodeA->start, mNodeA->n, P);
    PH.BuildKernelMatrix<Kernel, Point, PointArray>(mKernel, P, y);
    PH.GetColumn(0, ph);

    // c_i = tU_i' * k(X_i,x)
    //                  ||
    //                  ph
    if (mNodeA != Root) {
      mNodeTA->U.MatVec(ph, mNodeA->c, TRANSPOSE);
    }

    // Initial z0 = k(x,X_i) * tA_{ii} * k(X_i,x)
    DVector phh;
    mNodeTA->A.MatVec(ph, phh, NORMAL);
    z0 = ph.InProd(phh);

  }
  else {

    // Which child does y land on? The computation here is applicable
    // to only a binary tree constructed by using hyperplane
    // partitioning.
    Point y0;
    y.GetPoint(0, y0);
    double iprod = y0.InProd(mNodeA->normal);
    if (iprod < mNodeA->offset) {
      childA = mNodeA->LeftChild;
      childTA = mNodeTA->LeftChild;
    }
    else {
      childA = mNodeA->LeftChild->RightSibling;
      childTA = mNodeTA->LeftChild->RightSibling;
    }

    // Recurse on that child only
    BilinearImplicitUpward2<Kernel, Point, PointArray>
      (childA, childTA, X, y, mKernel, z0);

    // d_i = W_i' * d_j
    if (mNodeA != Root) {
      mNodeA->W.MatVec(childA->d, mNodeA->d, TRANSPOSE);
    }

  }

  if (mNodeA != Root) {

    // c_p = tW_p' * (sum_{j \in Ch(p)} c_j), part 1: summing c_j for itself
    DVector mc;
    if (parentA != Root) {
      mc = mNodeA->c;
    }

    siblingA = parentA->LeftChild;
    siblingTA = parentTA->LeftChild;
    for (k = 0; k < parentA->NumChild; k++, siblingA = siblingA->RightSibling,
           siblingTA = siblingTA->RightSibling) {
      if (siblingA == mNodeA) { continue; }

      // c_k = tTheta_k * d_i
      siblingTA->Theta.MatVec(mNodeA->d, siblingA->c, NORMAL);

      // z0 += d_i' * tXi_k * d_i + 2 * c_k' * tSigma_p * c_i
      DVector mXd;
      siblingTA->Xi.MatVec(mNodeA->d, mXd, NORMAL);
      z0 += mNodeA->d.InProd(mXd);
      DVector mSc;
      parentTA->Sigma.MatVec(mNodeA->c, mSc, NORMAL);
      z0 += 2.0 * siblingA->c.InProd(mSc);

      // c_p = tW_p' * (sum_{j \in Ch(p)} c_j), part 2: summing c_j for siblings
      if (parentA != Root) {
        mc.Add(siblingA->c);
      }
    }

    // c_p = tW_p' * (sum_{j \in Ch(p)} c_j), part 3: left-mult tW_p'
    if (parentA != Root) {
      parentTA->W.MatVec(mc, parentA->c, TRANSPOSE);
    }

  }

}


#endif

