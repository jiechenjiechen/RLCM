#include "BMatrix.hpp"

#define INITVAL_MaxNumChild 0
#define INITVAL_N           0
#define INITVAL_r           0


//--------------------------------------------------------------------------
BMatrix::
BMatrix() {
  Root = NULL;
  Init();
}


//--------------------------------------------------------------------------
void BMatrix::
Init(void) {
  ReleaseAllMemory();
  MaxNumChild = INITVAL_MaxNumChild;
  N           = INITVAL_N;
  r           = INITVAL_r;
}


//--------------------------------------------------------------------------
void BMatrix::
ReleaseAllMemory(void) {
  if (Root) {
    Root->DestroyTree();
    delete Root;
    Root = NULL;
  }
}


//--------------------------------------------------------------------------
BMatrix::
BMatrix(const BMatrix &G) {
  Root = NULL;
  Init();
  DeepCopy(G);
}


//--------------------------------------------------------------------------
BMatrix& BMatrix::
operator= (const BMatrix &G) {
  if (this != &G) {
    DeepCopy(G);
  }
  return *this;
}


//--------------------------------------------------------------------------
void BMatrix::
DeepCopy(const BMatrix &G) {
  Node *GRoot = G.GetRoot();
  if (Root) {
    ReleaseAllMemory();
  }
  if (GRoot) {
    GRoot->CopyTree(&Root);
  }
  MaxNumChild = G.MaxNumChild;
  N           = G.N;
  r           = G.r;
}


//--------------------------------------------------------------------------
BMatrix::
~BMatrix() {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
INTEGER BMatrix::
GetN(void) const {
  return N;
}


//--------------------------------------------------------------------------
Node* BMatrix::
GetRoot(void) const {
  return Root;
}


//--------------------------------------------------------------------------
INTEGER BMatrix::
GetMemEst(void) const {
  if (Root == NULL) {
    return 0;
  }
  INTEGER mem = Root->GetMemEstTreeAonly();
  return mem;
}


//--------------------------------------------------------------------------
void BMatrix::
PrintTree(void) const {
  if (Root == NULL) {
    printf("Tree is empty.\n");
    return;
  }
  Root->PrintTree(MaxNumChild, N);
}


//--------------------------------------------------------------------------
void BMatrix::
ConvertToDMatrix(DMatrix &B, MatState mState) const {
  if (Root == NULL) {
    printf("BMatrix::ConvertToDMatrix. Error: Matrix is empty. Function call takes no effect.\n");
    return;
  }
  B.Init(N, N);
  ConvertToDMatrixUpward(B, Root, mState);
}


//--------------------------------------------------------------------------
void BMatrix::
ConvertToDMatrixUpward(DMatrix &B, const Node *mNode, MatState mState) const {

  INTEGER i = 0;
  Node *childi = NULL;

  // If leaf, fill diagonal block
  if (mNode->NumChild == 0) {
    switch (mState) {
    case UNFACT: {
      B.SetBlock(mNode->start, mNode->n, mNode->start, mNode->n, mNode->A);
      break;
    }
    case LU_FACT: {
      DMatrix Eye(mNode->n);
      Eye.SetIdentity();
      DMatrix invA;
      mNode->A.DGETRS(Eye, invA, NORMAL);
      B.SetBlock(mNode->start, mNode->n, mNode->start, mNode->n, invA);
      break;
    }
    case CHOL_FACT: {
      DMatrix Eye(mNode->n);
      Eye.SetIdentity();
      DMatrix invG;
      mNode->A.DTRSM(Eye, invG, TRANSPOSE, UPPER);
      DMatrix invA;
      invG.MatMat(invG, invA, TRANSPOSE, NORMAL);
      B.SetBlock(mNode->start, mNode->n, mNode->start, mNode->n, invA);
      break;
    }
    }
    return;
  }

  // Recurse on children
  childi = mNode->LeftChild;
  for (i = 0; i < mNode->NumChild; i++, childi = childi->RightSibling) {
    ConvertToDMatrixUpward(B, childi, mState);
  }

}


//--------------------------------------------------------------------------
void BMatrix::
MatVec(const DVector &b, DVector &y, MatrixMode ModeA, MatState mState) const {

  if (Root == NULL) {
    printf("BMatrix::MatVec. Error: Matrix is empty. Function call takes no effect.\n");
    return;
  }
  else if (N != b.GetN()) {
    printf("BMatrix::MatVec. Error: Matrix/vector sizes are not compatible. Function call takes no effect.\n");
    return;
  }

  y.Init(N);
  MatVecUpward(b, y, Root, ModeA, mState);

}


//--------------------------------------------------------------------------
void BMatrix::
MatVecUpward(const DVector &b, DVector &y, Node *mNode, MatrixMode ModeA,
             MatState mState) const {

  INTEGER j = 0;
  Node *child = NULL;

  if (mNode->NumChild == 0) {

    // y_i = A_{ii} * b_i
    DVector mb;
    b.GetBlock(mNode->start, mNode->n, mb);
    DVector my;
    switch (mState) {
    case UNFACT: mNode->A.MatVec(mb, my, ModeA); break;
    case LU_FACT: mNode->A.DGETRS(mb, my, ModeA); break;
    case CHOL_FACT: mNode->A.DPOTRS(mb, my, UPPER); break;
    }
    y.SetBlock(mNode->start, mNode->n, my);

  }
  else {

    // Recurse on children
    child = mNode->LeftChild;
    for (j = 0; j < mNode->NumChild; j++, child = child->RightSibling) {
      MatVecUpward(b, y, child, ModeA, mState);
    }

  }

}


//--------------------------------------------------------------------------
void BMatrix::
MatMat(const DMatrix &B, DMatrix &Y, MatrixMode ModeA, MatState mState) const {

  if (Root == NULL) {
    printf("BMatrix::MatMat. Error: Matrix is empty. Function call takes no effect.\n");
    return;
  }
  else if (N != B.GetM()) {
    printf("BMatrix::MatMat. Error: Matrix sizes are not compatible. Function call takes no effect.\n");
    return;
  }

  INTEGER NCOL = B.GetN();
  Y.Init(N, NCOL);
  MatMatUpward(B, Y, Root, ModeA, mState);

}


//--------------------------------------------------------------------------
void BMatrix::
MatMatUpward(const DMatrix &B, DMatrix &Y, Node *mNode, MatrixMode ModeA,
             MatState mState) const {

  INTEGER j = 0;
  Node *child = NULL;
  INTEGER NCOL = B.GetN();

  if (mNode->NumChild == 0) {

    // Y_i = A_{ii} * B_i
    DMatrix mB;
    B.GetBlock(mNode->start, mNode->n, 0, NCOL, mB);
    DMatrix mY;
    switch (mState) {
    case UNFACT: mNode->A.MatMat(mB, mY, ModeA, NORMAL); break;
    case LU_FACT: mNode->A.DGETRS(mB, mY, ModeA); break;
    case CHOL_FACT: mNode->A.DPOTRS(mB, mY, UPPER); break;
    }
    Y.SetBlock(mNode->start, mNode->n, 0, NCOL, mY);

  }
  else {

    // Recurse on children
    child = mNode->LeftChild;
    for (j = 0; j < mNode->NumChild; j++, child = child->RightSibling) {
      MatMatUpward(B, Y, child, ModeA, mState);
    }

  }

}


//--------------------------------------------------------------------------
void BMatrix::
Invert(BMatrix &tA, MatrixType TypeA) const {

  if (Root == NULL) {
    printf("BMatrix::Invert. Error: Matrix is empty. Function call takes no effect.\n");
    return;
  }

  tA = *this; // Init
  Node *RootTA = tA.GetRoot();
  InvertUpward(RootTA, TypeA);

}


//--------------------------------------------------------------------------
void BMatrix::
InvertUpward(Node *mNodeTA, MatrixType TypeA) const {

  INTEGER j = 0;
  Node *childjTA = NULL;

  if (mNodeTA->NumChild == 0) {

    // Factorization
    switch (TypeA) {
    case GENERAL: mNodeTA->A.DGETRF(); break;
    case SPD: mNodeTA->A.DPOTRF(UPPER); break;
    }

  }
  else {

    // Recurse on children
    childjTA = mNodeTA->LeftChild;
    for (j = 0; j < mNodeTA->NumChild; j++, childjTA = childjTA->RightSibling) {
      InvertUpward(childjTA, TypeA);
    }

  }

}


//--------------------------------------------------------------------------
LogDet BMatrix::
Det(MatrixType TypeA, MatState mState) const {

  LogDet mLogDet = { 0.0, 1 };

  if (Root == NULL) {
    printf("BMatrix::Det. Error: Matrix is empty. Function call takes no effect.\n");
    return mLogDet;
  }
  if ((TypeA == GENERAL && mState == CHOL_FACT) ||
      (TypeA == SPD && mState == LU_FACT)) {
    printf("BMatrix::Det. Error: Invalid combination of matrix type and matrix state. Function call takes no effect.\n");
    return mLogDet;
  }

  mLogDet = DetUpward(Root, TypeA, mState);

  return mLogDet;

}


//--------------------------------------------------------------------------
LogDet BMatrix::
DetUpward(Node *mNode, MatrixType TypeA, MatState mState) const {

  INTEGER j = 0;
  Node *child = NULL;
  LogDet mLogDet = { 0.0, 1 };

  if (mNode->NumChild == 0) {

    mLogDet = mNode->A.Det(TypeA, mState);

  }
  else {

    // Recurse on children
    child = mNode->LeftChild;
    for (j = 0; j < mNode->NumChild; j++, child = child->RightSibling) {
      LogDet mLogDet2 = DetUpward(child, TypeA, mState);
      mLogDet.LogAbsDet += mLogDet2.LogAbsDet;
      mLogDet.Sign *= mLogDet2.Sign;
    }

  }

  return mLogDet;

}
