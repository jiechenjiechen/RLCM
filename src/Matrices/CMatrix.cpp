#include "CMatrix.hpp"

#define INITVAL_IsSym        false
#define INITVAL_MaxNumChild  0
#define INITVAL_N            0
#define INITVAL_r            0
#define INITVAL_mDiagCorrect 0.0


//--------------------------------------------------------------------------
CMatrix::
CMatrix() {
  Root = NULL;
  Init();
}


//--------------------------------------------------------------------------
void CMatrix::
Init(void) {
  ReleaseAllMemory();
  IsSym        = INITVAL_IsSym;
  MaxNumChild  = INITVAL_MaxNumChild;
  N            = INITVAL_N;
  r            = INITVAL_r;
  mDiagCorrect = INITVAL_mDiagCorrect;
}


//--------------------------------------------------------------------------
void CMatrix::
ReleaseAllMemory(void) {
  if (Root) {
    Root->DestroyTree();
    delete Root;
    Root = NULL;
  }
  Ut.ReleaseAllMemory();
  Vt.ReleaseAllMemory();
}


//--------------------------------------------------------------------------
CMatrix::
CMatrix(const CMatrix &G) {
  Root = NULL;
  Init();
  DeepCopy(G);
}


//--------------------------------------------------------------------------
CMatrix& CMatrix::
operator= (const CMatrix &G) {
  if (this != &G) {
    DeepCopy(G);
  }
  return *this;
}


//--------------------------------------------------------------------------
void CMatrix::
DeepCopy(const CMatrix &G) {
  Node *GRoot = G.GetRoot();
  if (Root) {
    ReleaseAllMemory();
  }
  if (GRoot) {
    GRoot->CopyTree(&Root);
  }
  IsSym        = G.IsSym;
  MaxNumChild  = G.MaxNumChild;
  N            = G.N;
  r            = G.r;
  Ut           = G.Ut;
  Vt           = G.Vt;
  mDiagCorrect = G.mDiagCorrect;
}


//--------------------------------------------------------------------------
CMatrix::
~CMatrix() {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
bool CMatrix::
IsSymmetric(void) const {
  return IsSym;
}


//--------------------------------------------------------------------------
void CMatrix::
SetSymmetric(bool IsSym_) {
  IsSym = IsSym_;
}


//--------------------------------------------------------------------------
INTEGER CMatrix::
GetN(void) const {
  return N;
}


//--------------------------------------------------------------------------
Node* CMatrix::
GetRoot(void) const {
  return Root;
}


//--------------------------------------------------------------------------
INTEGER CMatrix::
GetMemEst(void) const {
  if (Root == NULL) {
    return 0;
  }
  INTEGER mem = Root->GetMemEstTree();
  return mem;
}


//--------------------------------------------------------------------------
void CMatrix::
PrintTree(void) const {
  if (Root == NULL) {
    printf("Tree is empty.\n");
    return;
  }
  Root->PrintTree(MaxNumChild, N);
}


//--------------------------------------------------------------------------
void CMatrix::
ConvertToDMatrix(DMatrix &B) {
  if (Root == NULL) {
    printf("CMatrix::ConvertToDMatrix. Error: Matrix is empty. Function call takes no effect.\n");
    return;
  }
  Ut.Init(N,r);
  if (IsSym == false) { // If matrix unsymmetric, also need Vt
    Vt.Init(N,r);
  }
  B.Init(N,N);
  ConvertToDMatrixUpward(B, Root);
  Ut.ReleaseAllMemory();
  if (IsSym == false) {
    Vt.ReleaseAllMemory();
  }
}


//--------------------------------------------------------------------------
void CMatrix::
ConvertToDMatrixUpward(DMatrix &B, const Node *mNode) {

  INTEGER i = 0, j = 0;
  Node *childi = NULL, *childj = NULL;

  // If leaf, fill diagonal block and initialize Ut
  if (mNode->NumChild == 0) {
    B.SetBlock(mNode->start, mNode->n, mNode->start, mNode->n, mNode->A);
    Ut.SetBlock(mNode->start, mNode->n, 0, r, mNode->U);
    if (IsSym == false) { // If matrix unsymmetric, also need Vt
      Vt.SetBlock(mNode->start, mNode->n, 0, r, mNode->V);
    }
    return;
  }

  // Recurse on children
  childi = mNode->LeftChild;
  for (i = 0; i < mNode->NumChild; i++, childi = childi->RightSibling) {
    ConvertToDMatrixUpward(B, childi);
  }

  // Fill off-diagonal blocks
  childi = mNode->LeftChild;
  for (i = 0; i < mNode->NumChild; i++, childi = childi->RightSibling) {
    childj = mNode->LeftChild;
    for (j = 0; j < mNode->NumChild; j++, childj = childj->RightSibling) {

      if (childi == childj) { continue; }
      if (i > j && IsSym == true) { continue; }

      // Extract the corresponding parts of Ut
      DMatrix mU, mV;
      Ut.GetBlock(childi->start, childi->n, 0, r, mU);
      if (IsSym == true) {
        Ut.GetBlock(childj->start, childj->n, 0, r, mV);
      }
      else {
        Vt.GetBlock(childj->start, childj->n, 0, r, mV);
      }

      // U * Sigma * V' (V could be either U or V, depending on
      // whether the matrix is symmetric or not. See above)
      DMatrix mUS, mUSV;
      mU.MatMat(mNode->Sigma, mUS, NORMAL, NORMAL);
      mUS.MatMat(mV, mUSV, NORMAL, TRANSPOSE);

      // Set the part to B
      B.SetBlock(childi->start, childi->n, childj->start, childj->n, mUSV);

      // Set the symmetric part (only when the matrix is symmetric. If
      // not, program will not reach here.)
      if (IsSym == true) {
        DMatrix mUSVT;
        mUSV.Transpose(mUSVT);
        B.SetBlock(childj->start, childj->n, childi->start, childi->n, mUSVT);
      }

    }
  }

  // Transform Ut
  if (mNode == Root) { return; }

  childi = mNode->LeftChild;
  for (i = 0; i < mNode->NumChild; i++, childi = childi->RightSibling) {
    
    // Extract the corresponding parts of Ut
    DMatrix mU;
    Ut.GetBlock(childi->start, childi->n, 0, r, mU);

    // U * W
    DMatrix mU2;
    mU.MatMat(mNode->W, mU2, NORMAL, NORMAL);

    // Set the corresponding parts of Ut
    Ut.SetBlock(childi->start, childi->n, 0, r, mU2);

    // If matrix unsymmetric, also need Vt
    if (IsSym == false) {

      // Extract the corresponding parts of Vt
      DMatrix mV;
      Vt.GetBlock(childi->start, childi->n, 0, r, mV);

      // V * Z
      DMatrix mV2;
      mV.MatMat(mNode->Z, mV2, NORMAL, NORMAL);

      // Set the corresponding parts of Vt
      Vt.SetBlock(childi->start, childi->n, 0, r, mV2);

    }

  }

}


//--------------------------------------------------------------------------
void CMatrix::
MatVec(const DVector &b, DVector &y, MatrixMode ModeA, MatState mState) {

  if (Root == NULL) {
    printf("CMatrix::MatVec. Error: Matrix is empty. Function call takes no effect.\n");
    return;
  }
  else if (N != b.GetN()) {
    printf("CMatrix::MatVec. Error: Matrix/vector sizes are not compatible. Function call takes no effect.\n");
    return;
  }

  if (mState != UNFACT) {
    printf("CMatrix::MatVec. Error: mState must be UNFACT. Function call takes no effect.\n");
    return;
  }

  y.Init(N);
  MatVecInitAugmentedData(Root);
  MatVecUpward(b, y, Root, ModeA);
  MatVecDownward(b, y, Root, ModeA);
  MatVecReleaseAugmentedData(Root);

}


//--------------------------------------------------------------------------
void CMatrix::
MatVecInitAugmentedData(Node *mNode) {

  INTEGER i = 0;
  Node *child = NULL;

  if (mNode != Root) {
    mNode->c.Init(r);
    mNode->d.Init(r);
  }

  child = mNode->LeftChild;
  for (i = 0; i < mNode->NumChild; i++, child = child->RightSibling) {
    MatVecInitAugmentedData(child);
  }

}


//--------------------------------------------------------------------------
void CMatrix::
MatVecReleaseAugmentedData(Node *mNode) {

  INTEGER i = 0;
  Node *child = NULL;

  if (mNode != Root) {
    mNode->c.ReleaseAllMemory();
    mNode->d.ReleaseAllMemory();
  }

  child = mNode->LeftChild;
  for (i = 0; i < mNode->NumChild; i++, child = child->RightSibling) {
    MatVecReleaseAugmentedData(child);
  }

}


//--------------------------------------------------------------------------
void CMatrix::
MatVecUpward(const DVector &b, DVector &y, Node *mNode, MatrixMode ModeA) {

  INTEGER j = 0, k = 0;
  Node *child = NULL, *sibling = NULL;
  Node *parent = mNode->Parent;

  if (mNode->NumChild == 0) {

    // c_i = U_i' * b_i
    DVector mb;
    b.GetBlock(mNode->start, mNode->n, mb);
    if (IsSym == true) {
      mNode->U.MatVec(mb, mNode->c, TRANSPOSE);
    }
    else { // If matrix unsymmetric, should use V_i instead of U_i
      mNode->V.MatVec(mb, mNode->c, TRANSPOSE);
    }

    // y_i = A_{ii} * b_i
    DVector my;
    mNode->A.MatVec(mb, my, ModeA);
    y.SetBlock(mNode->start, mNode->n, my);

  }
  else {

    DVector mc;
    // c_i = W_i' * (sum_{j \in Ch(i)} c_j), part 0: init
    if (mNode != Root) {
      mc.Init(r);
    }

    child = mNode->LeftChild;
    for (j = 0; j < mNode->NumChild; j++, child = child->RightSibling) {

      // Recurse on children
      MatVecUpward(b, y, child, ModeA);
      if (mNode == Root) { continue; }

      // c_i = W_i' * (sum_{j \in Ch(i)} c_j), part 1: summing c_j
      mc.Add(child->c);

    }

    // c_i = W_i' * (sum_{j \in Ch(i)} c_j), part 2: left-mult W_i'
    if (mNode != Root) {
      if (IsSym == true) {
        mNode->W.MatVec(mc, mNode->c, TRANSPOSE);
      }
      else { // If matrix unsymmetric, should use Z instead of W
        mNode->Z.MatVec(mc, mNode->c, TRANSPOSE);
      }
    }

  }

  if (mNode == Root) { return; }

  // d_k = d_k + Sigma_p * c_i
  DVector dk;
  parent->Sigma.MatVec(mNode->c, dk, ModeA);

  sibling = parent->LeftChild;
  for (k = 0; k < parent->NumChild; k++, sibling = sibling->RightSibling) {
    if (sibling == mNode) { continue; }
    sibling->d.Add(dk);
  }

}


//--------------------------------------------------------------------------
void CMatrix::
MatVecDownward(const DVector &b, DVector &y, Node *mNode, MatrixMode ModeA) {

  INTEGER j = 0;
  Node *child = NULL;

  if (mNode->NumChild == 0) {

    // y_i = y_i + U_i * d_i
    DVector my;
    y.GetBlock(mNode->start, mNode->n, my);
    mNode->U.DGEMV(mNode->d, my, 1.0, 1.0, NORMAL);
    y.SetBlock(mNode->start, mNode->n, my);
    return;

  }

  // d_j = d_j + W_i * d_i
  DVector dj;
  if (mNode != Root) {
    mNode->W.MatVec(mNode->d, dj, NORMAL);
  }

  child = mNode->LeftChild;
  for (j = 0; j < mNode->NumChild; j++, child = child->RightSibling) {

    // d_j = d_j + W_i * d_i
    if (mNode != Root) {
      child->d.Add(dj);
    }

    // Recurse on children
    MatVecDownward(b, y, child, ModeA);

  }

}


//--------------------------------------------------------------------------
void CMatrix::
MatMat(const DMatrix &B, DMatrix &Y, MatrixMode ModeA) {

  if (Root == NULL) {
    printf("CMatrix::MatMat. Error: Matrix is empty. Function call takes no effect.\n");
    return;
  }
  else if (N != B.GetM()) {
    printf("CMatrix::MatMat. Error: Matrix sizes are not compatible. Function call takes no effect.\n");
    return;
  }

  INTEGER NCOL = B.GetN();
  Y.Init(N, NCOL);
  MatMatInitAugmentedData(Root, NCOL);
  MatMatUpward(B, Y, Root, ModeA);
  MatMatDownward(B, Y, Root, ModeA);
  MatMatReleaseAugmentedData(Root);

}


//--------------------------------------------------------------------------
void CMatrix::
MatMatInitAugmentedData(Node *mNode, INTEGER NCOL) {

  INTEGER i = 0;
  Node *child = NULL;

  if (mNode != Root) {
    mNode->C.Init(r, NCOL);
    mNode->D.Init(r, NCOL);
  }

  child = mNode->LeftChild;
  for (i = 0; i < mNode->NumChild; i++, child = child->RightSibling) {
    MatMatInitAugmentedData(child, NCOL);
  }

}


//--------------------------------------------------------------------------
void CMatrix::
MatMatReleaseAugmentedData(Node *mNode) {

  INTEGER i = 0;
  Node *child = NULL;

  if (mNode != Root) {
    mNode->C.ReleaseAllMemory();
    mNode->D.ReleaseAllMemory();
  }

  child = mNode->LeftChild;
  for (i = 0; i < mNode->NumChild; i++, child = child->RightSibling) {
    MatMatReleaseAugmentedData(child);
  }

}


//--------------------------------------------------------------------------
void CMatrix::
MatMatUpward(const DMatrix &B, DMatrix &Y, Node *mNode, MatrixMode ModeA) {

  INTEGER j = 0, k = 0;
  Node *child = NULL, *sibling = NULL;
  Node *parent = mNode->Parent;
  INTEGER NCOL = B.GetN();

  if (mNode->NumChild == 0) {

    // C_i = U_i' * B_i
    DMatrix mB;
    B.GetBlock(mNode->start, mNode->n, 0, NCOL, mB);
    if (IsSym == true) {
      mNode->U.MatMat(mB, mNode->C, TRANSPOSE, NORMAL);
    }
    else { // If matrix unsymmetric, should use V_i instead of U_i
      mNode->V.MatMat(mB, mNode->C, TRANSPOSE, NORMAL);
    }

    // Y_i = A_{ii} * B_i
    DMatrix mY;
    mNode->A.MatMat(mB, mY, ModeA, NORMAL);
    Y.SetBlock(mNode->start, mNode->n, 0, NCOL, mY);

  }
  else {

    DMatrix mC;
    // C_i = W_i' * (sum_{j \in Ch(i)} C_j), part 0: init
    if (mNode != Root) {
      mC.Init(r, NCOL);
    }

    child = mNode->LeftChild;
    for (j = 0; j < mNode->NumChild; j++, child = child->RightSibling) {

      // Recurse on children
      MatMatUpward(B, Y, child, ModeA);
      if (mNode == Root) { continue; }

      // C_i = W_i' * (sum_{j \in Ch(i)} C_j), part 1: summing C_j
      mC.Add(child->C);

    }

    // C_i = W_i' * (sum_{j \in Ch(i)} C_j), part 2: left-mult W_i'
    if (mNode != Root) {
      if (IsSym == true) {
        mNode->W.MatMat(mC, mNode->C, TRANSPOSE, NORMAL);
      }
      else { // If matrix unsymmetric, should use Z instead of W
        mNode->Z.MatMat(mC, mNode->C, TRANSPOSE, NORMAL);
      }
    }

  }

  if (mNode == Root) { return; }

  // D_k = D_k + Sigma_p * C_i
  DMatrix Dk;
  parent->Sigma.MatMat(mNode->C, Dk, ModeA, NORMAL);

  sibling = parent->LeftChild;
  for (k = 0; k < parent->NumChild; k++, sibling = sibling->RightSibling) {
    if (sibling == mNode) { continue; }
    sibling->D.Add(Dk);
  }

}


//--------------------------------------------------------------------------
void CMatrix::
MatMatDownward(const DMatrix &B, DMatrix &Y, Node *mNode, MatrixMode ModeA) {

  INTEGER j = 0;
  Node *child = NULL;
  INTEGER NCOL = B.GetN();

  if (mNode->NumChild == 0) {

    // Y_i = Y_i + U_i * D_i
    DMatrix mY;
    Y.GetBlock(mNode->start, mNode->n, 0, NCOL, mY);
    mNode->U.DGEMM(mNode->D, mY, 1.0, 1.0, NORMAL, NORMAL);
    Y.SetBlock(mNode->start, mNode->n, 0, NCOL, mY);
    return;

  }

  // D_j = D_j + W_i * D_i
  DMatrix Dj;
  if (mNode != Root) {
    mNode->W.MatMat(mNode->D, Dj, NORMAL, NORMAL);
  }

  child = mNode->LeftChild;
  for (j = 0; j < mNode->NumChild; j++, child = child->RightSibling) {

    // D_j = D_j + W_i * D_i
    if (mNode != Root) {
      child->D.Add(Dj);
    }

    // Recurse on children
    MatMatDownward(B, Y, child, ModeA);

  }

}


//--------------------------------------------------------------------------
void CMatrix::
Invert(CMatrix &tA) {

  if (Root == NULL) {
    printf("CMatrix::Invert. Error: Matrix is empty. Function call takes no effect.\n");
    return;
  }

  tA = *this; // Init (keeps the IsSym tag unchanged)
  Node *RootTA = tA.GetRoot();
  InvertInitAugmentedData(RootTA);
  InvertUpward(Root, RootTA);
  InvertDownward(RootTA);
  InvertReleaseAugmentedData(RootTA);

}


//--------------------------------------------------------------------------
void CMatrix::
InvertInitAugmentedData(Node *mNodeTA) {

  INTEGER i = 0;
  Node *childTA = NULL;

  if (mNodeTA->NumChild != 0) {
    mNodeTA->E.Init(r);
    mNodeTA->Xi.Init(r);
  }
  if (mNodeTA != Root) {
    mNodeTA->Theta.Init(r);
  }

  childTA = mNodeTA->LeftChild;
  for (i = 0; i < mNodeTA->NumChild; i++, childTA = childTA->RightSibling) {
    InvertInitAugmentedData(childTA);
  }

}


//--------------------------------------------------------------------------
void CMatrix::
InvertReleaseAugmentedData(Node *mNodeTA) {

  INTEGER i = 0;
  Node *childTA = NULL;

  if (mNodeTA->NumChild != 0) {
    mNodeTA->E.ReleaseAllMemory();
    mNodeTA->Xi.ReleaseAllMemory();
  }
  if (mNodeTA != Root) {
    mNodeTA->Theta.ReleaseAllMemory();
  }

  childTA = mNodeTA->LeftChild;
  for (i = 0; i < mNodeTA->NumChild; i++, childTA = childTA->RightSibling) {
    InvertReleaseAugmentedData(childTA);
  }

}


//--------------------------------------------------------------------------
void CMatrix::
InvertUpward(Node *mNodeA, Node *mNodeTA) {

  INTEGER j = 0;
  Node *childA = NULL, *childTA = NULL;
  Node *parentA = mNodeA->Parent;

  if (mNodeA->NumChild == 0) {

    // B = A_{ii} - U_i * Sigma_p * U_i'
    DMatrix mUS, mUSU, mB;
    if (mNodeA != Root) {
      mNodeA->U.MatMat(parentA->Sigma, mUS, NORMAL, NORMAL);
      if (IsSym == true) { // If A symmetric, symmetrize result
        mUS.MatMat(mNodeA->U, mUSU, NORMAL, TRANSPOSE);
        mUSU.Symmetrize();
      }
      else { // If A unsymmetric, use V in place of U
        mUS.MatMat(mNodeA->V, mUSU, NORMAL, TRANSPOSE);
      }
      mNodeA->A.Subtract(mUSU, mB);
    }
    else {
      mB = mNodeA->A;
    }

    if (IsSym == true) { // If A symmetric, use Cholesky factorization

      // tA_{ii} = inv(B)
      //         = inv(G)' * inv(G) if B = GG'
      mB.DPOTRF(LOWER);
      DMatrix &mG = mB;
      DMatrix Eye(mNodeA->n);
      Eye.SetIdentity();
      DMatrix invG;
      mG.DTRSM(Eye, invG, NORMAL, LOWER);
      invG.MatMat(invG, mNodeTA->A, TRANSPOSE, NORMAL);

      // tU_i = inv(B) * U_i
      //      = inv(G)' * inv(G) * U_i = inv(G)' * invGU
      DMatrix invGU;
      mG.DTRSM(mNodeA->U, invGU, NORMAL, LOWER);
      mG.DTRSM(invGU, mNodeTA->U, TRANSPOSE, LOWER);

      // tTheta_i = U_i' * tU_i
      //          = invGU' * invGU
      invGU.MatMat(invGU, mNodeTA->Theta, TRANSPOSE, NORMAL);

      // det(B)
      mNodeA->mLogDet = mB.Det(SPD, CHOL_FACT);

    }
    else { // If A unsymmetric, use LU factorization

      // tA_{ii} = inv(B)
      mB.DGETRF();
      DMatrix Eye(mNodeA->n);
      Eye.SetIdentity();
      mB.DGETRS(Eye, mNodeTA->A, NORMAL);

      // tU_i = tA_{ii} * U_i
      // tV_i = tA_{ii}' * V_i
      mB.DGETRS(mNodeA->U, mNodeTA->U, NORMAL);
      mB.DGETRS(mNodeA->V, mNodeTA->V, TRANSPOSE);

      // tTheta_i = V_i' * tU_i
      mNodeA->V.MatMat(mNodeTA->U, mNodeTA->Theta, TRANSPOSE, NORMAL);

      // det(B)
      mNodeA->mLogDet = mB.Det(GENERAL, LU_FACT);

    }

    return;

  }

  childA = mNodeA->LeftChild;
  childTA = mNodeTA->LeftChild;
  for (j = 0; j < mNodeA->NumChild; j++, childA = childA->RightSibling,
         childTA = childTA->RightSibling) {

    // Recurse on children
    InvertUpward(childA, childTA);

    if (childA->NumChild == 0) { continue; }

    // tW_j = ( I + tSigma_j * tXi_j ) * W_j
    DMatrix mST;
    childTA->Sigma.MatMat(childTA->Xi, mST, NORMAL, NORMAL);
    mST.AddDiagonal(1.0);
    mST.MatMat(childA->W, childTA->W, NORMAL, NORMAL);

    if (IsSym == false) { // If unsymmetric, also need tZ_j
      childTA->Sigma.MatMat(childTA->Xi, mST, TRANSPOSE, TRANSPOSE);
      mST.AddDiagonal(1.0);
      mST.MatMat(childA->Z, childTA->Z, NORMAL, NORMAL);
    }

    // tTheta_j = W_j' * tXi_j * tW_j
    DMatrix mWX;
    if (IsSym == true) {
      childA->W.MatMat(childTA->Xi, mWX, TRANSPOSE, NORMAL);
    }
    else { // If A unsymmetric, replace W_j by Z_j
      childA->Z.MatMat(childTA->Xi, mWX, TRANSPOSE, NORMAL);
    }
    mWX.MatMat(childTA->W, childTA->Theta, NORMAL, NORMAL);
    if (IsSym == true) { // If A symmetric, tTheta_j is symmetric
      childTA->Theta.Symmetrize();
    }

  }

  // tXi_i = sum_j tTheta_j
  mNodeTA->Xi.Init(r);
  childTA = mNodeTA->LeftChild;
  for (j = 0; j < mNodeTA->NumChild; j++, childTA = childTA->RightSibling) {
    mNodeTA->Xi.Add(childTA->Theta);
  }

  // Lambda = Sigma_i - W_i * Sigma_p * W_i'
  DMatrix mLambda, mWS, mWSW;
  if (mNodeA != Root) {
    mNodeA->W.MatMat(parentA->Sigma, mWS, NORMAL, NORMAL);
    if (IsSym == true) { // If A symmetric, W_i * Sigma_p * W_i' is symmetric
      mWS.MatMat(mNodeA->W, mWSW, NORMAL, TRANSPOSE);
      mWSW.Symmetrize();
    }
    else { // If A unsymmetric, replace second W_i by Z_i
      mWS.MatMat(mNodeA->Z, mWSW, NORMAL, TRANSPOSE);
    }
    mNodeA->Sigma.Subtract(mWSW, mLambda);
  }
  else {
    mLambda = mNodeA->Sigma;
  }

  // H = I + Lambda * tXi_i
  DMatrix mH;
  mLambda.MatMat(mNodeTA->Xi, mH, NORMAL, NORMAL);
  mH.AddDiagonal(1.0);

  // det(H)
  mH.DGETRF();
  mNodeA->mLogDet = mH.Det(GENERAL, LU_FACT);

  // tSigma_i = -inv(H) * Lambda
  mH.DGETRS(mLambda, mNodeTA->Sigma, NORMAL);
  if (IsSym == true) { // If A symmetric, tSigma_i is symmetric
    mNodeTA->Sigma.Symmetrize();
  }
  mNodeTA->Sigma.Negate();

  // tE_j = tW_j * tSigma_i * tW_j'
  childTA = mNodeTA->LeftChild;
  for (j = 0; j < mNodeTA->NumChild; j++, childTA = childTA->RightSibling) {
    if (childTA->NumChild == 0) { continue; }
    childTA->W.MatMat(mNodeTA->Sigma, mWS, NORMAL, NORMAL);
    if (IsSym == true) { // If A symmetric, tE_j is symmetric
      mWS.MatMat(childTA->W, childTA->E, NORMAL, TRANSPOSE);
      childTA->E.Symmetrize();
    }
    else { // If A unsymmetric, replace the second tW_j by tZ_j
      mWS.MatMat(childTA->Z, childTA->E, NORMAL, TRANSPOSE);
    }
  }

  // No need to set tE for root because it is zero by initialization

}


//--------------------------------------------------------------------------
void CMatrix::
InvertDownward(Node *mNodeTA) {

  INTEGER j = 0;
  Node *childTA = NULL;
  Node *parentTA = mNodeTA->Parent;

  if (mNodeTA->NumChild == 0) {

    // tA_{ii} = tA_{ii} + tU_i * tSigma_p * tU_i'
    DMatrix mUS, mUSU;
    mNodeTA->U.MatMat(parentTA->Sigma, mUS, NORMAL, NORMAL);
    if (IsSym == true) { // If A symmetric, tU_i * tSigma_p * tU_i' is symmetric
      mUS.MatMat(mNodeTA->U, mUSU, NORMAL, TRANSPOSE);
      mUSU.Symmetrize();
    }
    else { // If A unsymmetric, replace second tU_i by tV_i
      mUS.MatMat(mNodeTA->V, mUSU, NORMAL, TRANSPOSE);
    }
    mNodeTA->A.Add(mUSU);
    return;

  }

  // tE_i = tE_i + tW_i * tE_p * tW_i'
  // tSigma_i = tSigma_i + tE_i
  DMatrix mWE, mWEW;
  if (parentTA != NULL) {
    mNodeTA->W.MatMat(parentTA->E, mWE, NORMAL, NORMAL);
    if (IsSym == true) { // If A symmetric, tW_i * tE_p * tW_i' is symmetric
      mWE.MatMat(mNodeTA->W, mWEW, NORMAL, TRANSPOSE);
      mWEW.Symmetrize();
    }
    else { // If A unsymmetric, replace second tW_i by tZ_i
      mWE.MatMat(mNodeTA->Z, mWEW, NORMAL, TRANSPOSE);
    }
    mNodeTA->E.Add(mWEW);
  }
  mNodeTA->Sigma.Add(mNodeTA->E);

  // Recurse on children
  childTA = mNodeTA->LeftChild;
  for (j = 0; j < mNodeTA->NumChild; j++, childTA = childTA->RightSibling) {
    InvertDownward(childTA);
  }

}


//--------------------------------------------------------------------------
LogDet CMatrix::
Det(void) const {

  LogDet mLogDet = { 0.0, 1 };

  if (Root == NULL) {
    printf("CMatrix::Det. Error: Matrix is empty. Function call takes no effect.\n");
    return mLogDet;
  }

  mLogDet = DetUpward(Root);

  return mLogDet;

}


//--------------------------------------------------------------------------
LogDet CMatrix::
DetUpward(Node *mNode) const {

  INTEGER j = 0;
  Node *child = NULL;
  LogDet mLogDet = { mNode->mLogDet.LogAbsDet, mNode->mLogDet.Sign };

  if (mNode->NumChild != 0) {
    child = mNode->LeftChild;
    for (j = 0; j < mNode->NumChild; j++, child = child->RightSibling) {
      LogDet mLogDet2 = DetUpward(child);
      mLogDet.LogAbsDet += mLogDet2.LogAbsDet;
      mLogDet.Sign *= mLogDet2.Sign;
    }
  }

  return mLogDet;

}


//--------------------------------------------------------------------------
void CMatrix::
Sqrt(CMatrix &G) const {

  if (Root == NULL) {
    printf("CMatrix::Sqrt. Error: Matrix is empty. Function call takes no effect.\n");
    return;
  }

  G = *this; // Init
  G.SetSymmetric(false);
  Node *RootG = G.GetRoot();
  SqrtInitAugmentedData(RootG);
  SqrtUpward(Root, RootG);
  SqrtDownward(RootG);
  SqrtReleaseAugmentedData(RootG);

}


//--------------------------------------------------------------------------
void CMatrix::
SqrtInitAugmentedData(Node *mNodeG) const {

  INTEGER i = 0;
  Node *childG = NULL;

  if (mNodeG->NumChild != 0) {
    mNodeG->E.Init(r);
    mNodeG->Xi.Init(r);
  }
  if (mNodeG != Root) {
    mNodeG->Theta.Init(r);
  }

  childG = mNodeG->LeftChild;
  for (i = 0; i < mNodeG->NumChild; i++, childG = childG->RightSibling) {
    SqrtInitAugmentedData(childG);
  }

}


//--------------------------------------------------------------------------
void CMatrix::
SqrtReleaseAugmentedData(Node *mNodeG) const {

  INTEGER i = 0;
  Node *childG = NULL;

  if (mNodeG->NumChild != 0) {
    mNodeG->E.ReleaseAllMemory();
    mNodeG->Xi.ReleaseAllMemory();
  }
  if (mNodeG != Root) {
    mNodeG->Theta.ReleaseAllMemory();
  }

  childG = mNodeG->LeftChild;
  for (i = 0; i < mNodeG->NumChild; i++, childG = childG->RightSibling) {
    SqrtReleaseAugmentedData(childG);
  }

}


//--------------------------------------------------------------------------
void CMatrix::
SqrtUpward(Node *mNodeA, Node *mNodeG) const {

  INTEGER j = 0;
  Node *childA = NULL, *childG = NULL;
  Node *parentA = mNodeA->Parent;

  if (mNodeA->NumChild == 0) {

    // B = A_{ii} - U_i * Sigma_p * U_i'
    DMatrix mUS, mUSU, mB;
    if (mNodeA != Root) {
      mNodeA->U.MatMat(parentA->Sigma, mUS, NORMAL, NORMAL);
      mUS.MatMat(mNodeA->U, mUSU, NORMAL, TRANSPOSE);
      mUSU.Symmetrize();
      mNodeA->A.Subtract(mUSU, mB);
    }
    else {
      mB = mNodeA->A;
    }

    // B = G_{ii} * G_{ii}'
    mB.Chol(mNodeG->A, LOWER);

    // V_i = G_{ii} \ U_i
    mNodeG->A.DTRSM(mNodeG->U, mNodeG->V, NORMAL, LOWER);
    // Note: The function DMatrix::Chol is implemented by first
    // calling DPOTRF, followed by zeroing out the irrelavant
    // triangular part. Hence, calling DTRSM here should be safe, as
    // long as the irrelavant triangular part is not referenced. An
    // advantage of DTRSM over Mldivide is that it is cheaper.

    // Theta_i = V_i' * V_i
    mNodeG->V.MatMat(mNodeG->V, mNodeG->Theta, TRANSPOSE, NORMAL);

    return;

  }

  childA = mNodeA->LeftChild;
  childG = mNodeG->LeftChild;
  for (j = 0; j < mNodeA->NumChild; j++, childA = childA->RightSibling,
         childG = childG->RightSibling) {

    // Recurse on children
    SqrtUpward(childA, childG);

    if (childA->NumChild == 0) { continue; }

    // Z_j = ( I + Omega_j * Xi_j ) \ W_j
    DMatrix mOX;
    childG->Sigma.MatMat(childG->Xi, mOX, NORMAL, NORMAL);
    mOX.AddDiagonal(1.0);
    mOX.Mldivide(childG->W, childG->Z, NORMAL, GENERAL);

    // Theta_j = Z_j' * Xi_j * Z_j
    DMatrix mZX;
    childG->Z.MatMat(childG->Xi, mZX, TRANSPOSE, NORMAL);
    mZX.MatMat(childG->Z, childG->Theta, NORMAL, NORMAL);
    childG->Theta.Symmetrize();

  }

  // Xi_i = sum_j Theta_j
  mNodeG->Xi.Init(r);
  childG = mNodeG->LeftChild;
  for (j = 0; j < mNodeG->NumChild; j++, childG = childG->RightSibling) {
    mNodeG->Xi.Add(childG->Theta);
  }

  // Lambda_i = Sigma_i - W_i * Sigma_p * W_i'
  DMatrix mLambda, mWS, mWSW;
  if (mNodeA != Root) {
    mNodeA->W.MatMat(parentA->Sigma, mWS, NORMAL, NORMAL);
    mWS.MatMat(mNodeA->W, mWSW, NORMAL, TRANSPOSE);
    mWSW.Symmetrize();
    mNodeA->Sigma.Subtract(mWSW, mLambda);
  }
  else {
    mLambda = mNodeA->Sigma;
  }

  // Solve Lambda_i = Omega_i' + Omega_i + Omega_i * Xi_i * Omega_i' for Omega_i
  SqrtSolveRiccati(mLambda, mNodeG->Xi, mNodeG->Sigma);

  // E_j = W_j * Omega_i * Z_j'
  childG = mNodeG->LeftChild;
  for (j = 0; j < mNodeG->NumChild; j++, childG = childG->RightSibling) {
    if (childG->NumChild == 0) { continue; }
    childG->W.MatMat(mNodeG->Sigma, mWS, NORMAL, NORMAL);
    mWS.MatMat(childG->Z, childG->E, NORMAL, TRANSPOSE);
  }

  // No need to set E for root because it is zero by initialization

}


//--------------------------------------------------------------------------
void CMatrix::
SqrtDownward(Node *mNodeG) const {

  INTEGER j = 0;
  Node *childG = NULL;
  Node *parentG = mNodeG->Parent;

  if (mNodeG->NumChild == 0) {

    // G_{ii} = G_{ii} + U_i * Omega_p * V_i'
    DMatrix mUO, mUOV;
    mNodeG->U.MatMat(parentG->Sigma, mUO, NORMAL, NORMAL);
    mUO.MatMat(mNodeG->V, mUOV, NORMAL, TRANSPOSE);
    mNodeG->A.Add(mUOV);
    return;

  }

  // E_i = E_i + W_i * E_p * Z_i'
  // Omega_i = Omega_i + E_i
  DMatrix mWE, mWEZ;
  if (parentG != NULL) {
    mNodeG->W.MatMat(parentG->E, mWE, NORMAL, NORMAL);
    mWE.MatMat(mNodeG->Z, mWEZ, NORMAL, TRANSPOSE);
    mNodeG->E.Add(mWEZ);
  }
  mNodeG->Sigma.Add(mNodeG->E);

  // Recurse on children
  childG = mNodeG->LeftChild;
  for (j = 0; j < mNodeG->NumChild; j++, childG = childG->RightSibling) {
    SqrtDownward(childG);
  }

}


//--------------------------------------------------------------------------
void CMatrix::
SqrtSolveRiccati(const DMatrix &Lambda, const DMatrix &Xi, DMatrix &D) const {

  DMatrix Eye(r), H(r*2), S, Q, Q21, Q11;
  Eye.SetIdentity();

  // H = [I Xi; Lambda -I]
  H.SetBlock(0, r, 0, r, Eye);
  H.SetBlock(0, r, r, r, Xi);
  H.SetBlock(r, r, 0, r, Lambda);
  Eye.Negate();
  H.SetBlock(r, r, r, r, Eye);

  // HQ = QS
  H.RealSchur(S, Q, RHP);
  Q.GetBlock(0, r, 0, r, Q11);
  Q.GetBlock(r, r, 0, r, Q21);

  // D = Q11' \ Q21'
  Q21.Transpose();
  Q11.Mldivide(Q21, D, TRANSPOSE, GENERAL);

}


//--------------------------------------------------------------------------
void CMatrix::
SamplingPivotsRegularGrid(const double *BBox, INTEGER d, DPointArray &P) {

  // Step: Compute grid dim
  Elm3 *GridInfo = NULL;
  New_1D_Array<Elm3, INTEGER>(&GridInfo, d);
  for (INTEGER i = 0; i < d; i++) {
    GridInfo[i].idx = i;
    GridInfo[i].len = BBox[i+d] - BBox[i];
    GridInfo[i].low = BBox[i];
  }
  INTEGER GridDimTotal = 1;
  for (INTEGER j = 0; j < d; j++) {
    double base = 1.0;
    for (INTEGER i = j; i < d; i++) {
      base *= GridInfo[i].len;
    }
    double alpha = pow((double)r/GridDimTotal/base, 1.0/(d-j));
    GridInfo[j].dim = (INTEGER)ceil(GridInfo[j].len * alpha);
    GridInfo[j].inc = GridInfo[j].len / GridInfo[j].dim;
    GridDimTotal *= GridInfo[j].dim;
  }

  // Step: Sort grid dim (in the increasing order)
  qsort(GridInfo, d, sizeof(Elm3), CompareElm3Dim);

  // Step: Compute GridDimCum
  INTEGER *GridDimCum = NULL;
  New_1D_Array<INTEGER, INTEGER>(&GridDimCum, d+1);
  GridDimCum[0] = 1;
  for (INTEGER j = 0; j < d; j++) {
    GridDimCum[j+1] = GridDimCum[j] * GridInfo[j].dim;
  }

  // Step: Setup a point array PP, whose columns correspond to the
  // sorted GridInfo.dim. The first dimension varies the fastest.
  DPointArray PP(GridDimTotal, d); // Init 0
  double *mPP = PP.GetPointer();

  // Step: Fill PP with coordinates. The exception is that the
  // first column is filled with integers, because they will
  // undergo additional processing in the next step. Grid points
  // are marked one by one, with the highest priority going to the
  // longest dimension (that is, the last column).
  INTEGER *CurValInt = NULL;
  New_1D_Array<INTEGER, INTEGER>(&CurValInt, d);
  double *CurValDbl = NULL;
  New_1D_Array<double, INTEGER>(&CurValDbl, d);
  for (INTEGER i = 0; i < d; i++) {
    CurValInt[i] = 1;
    CurValDbl[i] = 0.5 * GridInfo[i].inc + GridInfo[i].low;
  }
  for (INTEGER i = 0; i < r; i++) {
    INTEGER loc = 0;
    for (INTEGER j = 0; j < d; j++) {
      loc += GridDimCum[j] * (CurValInt[j]-1);
    }
    mPP[loc] = (double)CurValInt[0]; // integer for the 1st column
    for (INTEGER j = 1; j < d; j++) {
      mPP[loc + j*GridDimTotal] = (double)CurValDbl[j];
    }
    for (INTEGER j = d-1; j >= 0; j--) {
      if (CurValInt[j] != GridInfo[j].dim) {
        CurValInt[j] ++;
        CurValDbl[j] += GridInfo[j].inc;
        break;
      }
      else {
        CurValInt[j] = 1;
        CurValDbl[j] = 0.5 * GridInfo[j].inc + GridInfo[j].low;
      }
    }
  }

  // Step: Copy the information in PP to P, resolving the right
  // order of dimensions. Convert from temporary integer
  // coordinates to floating points for the first column of
  // PP. Along this dimension, not every grid point is marked; so
  // we need special calculation for the coordinates. Some rows of
  // PP are empty, because the corresponding grid points are not
  // marked.
  P.Init(r, d);
  double *mP = P.GetPointer();
  INTEGER lP = 0, lPP = 0;
  for (INTEGER i = 0; i < GridDimTotal/GridInfo[0].dim; i++) {
    INTEGER ActualTotal = GridInfo[0].dim;
    for (INTEGER j = GridInfo[0].dim-1; j >= 0; j--) {
      if (mPP[lPP+j] == 0.0) {
        ActualTotal--;
      }
      else {
        break;
      }
    }
    for (INTEGER j = 0; j < ActualTotal; j++) {
      GridInfo[0].inc = GridInfo[0].len / ActualTotal;
      mPP[lPP+j] = ( -0.5 + mPP[lPP+j] ) * GridInfo[0].inc + GridInfo[0].low;
    }
    for (INTEGER j = 0; j < ActualTotal; j++) {
      for (INTEGER k = 0; k < d; k++) {
        mP[lP+j+GridInfo[k].idx*r] = mPP[lPP+j+k*GridDimTotal];
      }
    }
    lP += ActualTotal;
    lPP += GridInfo[0].dim;
  }

  // Step: Clean up
  Delete_1D_Array<Elm3>(&GridInfo);
  Delete_1D_Array<INTEGER>(&GridDimCum);
  Delete_1D_Array<INTEGER>(&CurValInt);
  Delete_1D_Array<double>(&CurValDbl);

}
