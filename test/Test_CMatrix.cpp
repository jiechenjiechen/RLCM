// Unit test of the class CMatrix.

#include "LibCMatrix.hpp"

void PrintPivots(const Node *mNode);
void ApproximationError(const Node *mNode, const DMatrix &DA, const DMatrix &B);

int main(int argc, char **argv) {

  INTEGER NumThreads = String2Integer(argv[1]);
#ifdef USE_OPENBLAS
  openblas_set_num_threads(NumThreads);
#elif defined USE_OPENMP
  omp_set_num_threads(NumThreads);
#else
  NumThreads = 1; // To avoid compiler warining of unused variable
#endif
  printf("+---------------------------------------------------------------+\n");
  printf("|  Test_CMatrix.cpp                                             |\n");
  printf("+---------------------------------------------------------------+\n");

  char *mParString = argv[2];
  PartMethod mPar = RAND;
  if (0 == strcmp(mParString, "RAND")) {
    mPar = RAND;
  }
  else if (0 == strcmp(mParString, "PCA")) {
    mPar = PCA;
  }
  else if (0 == strcmp(mParString, "BBOX")) {
    mPar = BBOX;
  }
  else {
    printf("Unsupported partitioning method!\n");
    exit(1);
  }

  // Part 1: Build matrix from point set --------------------------------------
  printf("\nPart 1: Build matrix from point set\n\n");

  // Read in data
  DPointArray X;
  DVector y;
  INTEGER d = 2;
  LibSVM_IO::ReadData("Toy_bin.train", X, y, d);

  // Save the original data
  DPointArray X_original = X;

  // Build the (approximate) matrix
  CMatrix A;
  INTEGER N = X.GetN();
  INTEGER *Perm = NULL, *iPerm = NULL;
  New_1D_Array<INTEGER, INTEGER>(&Perm, N);
  New_1D_Array<INTEGER, INTEGER>(&iPerm, N);
  INTEGER r = N/4;
  INTEGER N0;
  switch (mPar) {
  case RAND:
  case PCA:
    N0 = N/4; // Leaf size
    break;
  case BBOX:
    N0 = 2;   // Number of levels
    break;
  }
  double mDiagCorrect = 1e-8;
  unsigned Seed = 1;
  A.BuildTree<DPoint, DPointArray>(X, Perm, iPerm, r, N0,
                                   mDiagCorrect, Seed, mPar);

  double s = 1.0;
  double sigma = 0.3;
  IsotropicGaussian mKernel(s, sigma);
  double lambda = 1e-4;
  A.BuildKernelMatrix<IsotropicGaussian, DPoint, DPointArray>
    (mKernel, X, lambda);

  // Build the nonapproximate matrix (must be done after building A,
  // because X is permuted)
  DMatrix K;
  K.BuildKernelMatrix<IsotropicGaussian, DPoint, DPointArray>
    (mKernel, X, lambda);

  // Inspect tree structure
  printf("Tree structure ----------\n");
  A.PrintTree();

  // Inspect pivots
  if (mPar == BBOX) {
    printf("Pivots ----------\n");
    PrintPivots(A.GetRoot());
  }

  // Check permutation
  printf("Check permutation ----------\n");
  double *X_perm_back = NULL;
  New_1D_Array<double, INTEGER>(&X_perm_back, N*d);
  double *X_ptr = X.GetPointer();
  for (INTEGER i = 0; i < N; i++) {
    for (INTEGER j = 0; j < d; j++) {
      X_perm_back[Perm[i]+j*N] = X_ptr[i+j*N];
    }
  }
  double derr = Diff1(X_perm_back, X_original.GetPointer(), N*d);
  printf("Discrepancy in permuted data points: %g\n", derr);

  X_ptr = X_original.GetPointer();
  for (INTEGER i = 0; i < N; i++) {
    for (INTEGER j = 0; j < d; j++) {
      X_perm_back[iPerm[i]+j*N] = X_ptr[i+j*N];
    }
  }
  derr = Diff1(X_perm_back, X.GetPointer(), N*d);
  printf("Discrepancy in permuted data points %g\n", derr);

  // Convert to dense matrix
  DMatrix DA;
  A.ConvertToDMatrix(DA);

  // Approximation error
  printf("Matrix approximation ----------\n");
  ApproximationError(A.GetRoot(), DA, K);

  // Matrix-vector multiplication
  printf("Matrix-vector multiplication ----------\n");
  printf("[K: Kernel matrix; A: RLCM; DA: dense form of A]\n");
  DVector b(N);
  b.SetStandardNormal();
  DVector y1, y2, y3, y4, ydelta;
  A.MatVec(b, y1, NORMAL);
  A.MatVec(b, y2, TRANSPOSE);
  DA.MatVec(b, y3, NORMAL);
  K.MatVec(b, y4, NORMAL);

  y1.Subtract(y2, ydelta);
  derr = ydelta.Norm2()/y1.Norm2();
  printf("Discprepancy between A*b and A'*b: %g\n", derr);
  y1.Subtract(y3, ydelta);
  derr = ydelta.Norm2()/y1.Norm2();
  printf("Discprepancy between A*b and DA*b: %g\n", derr);
  y1.Subtract(y4, ydelta);
  derr = ydelta.Norm2()/y1.Norm2();
  printf("Relative difference between A*b and K*b: %g\n", derr);

  // Matrix-matrix multiplication
  printf("Matrix-matrix multiplication ----------\n");
  printf("[K: Kernel matrix; A: RLCM; DA: dense form of A]\n");
  INTEGER M = 2;
  DMatrix B(N, M);
  B.SetStandardNormal();
  DMatrix Y1, Y2, Y3, Y4, Ydelta;
  A.MatMat(B, Y1, NORMAL);
  A.MatMat(B, Y2, TRANSPOSE);
  DA.MatMat(B, Y3, NORMAL, NORMAL);
  K.MatMat(B, Y4, NORMAL, NORMAL);

  Y1.Subtract(Y2, Ydelta);
  derr = Ydelta.NormF()/Y1.NormF();
  printf("Discprepancy between A*B and A'*B: %g\n", derr);
  Y1.Subtract(Y3, Ydelta);
  derr = Ydelta.NormF()/Y1.NormF();
  printf("Discprepancy between A*B and DA*B: %g\n", derr);
  Y1.Subtract(Y4, Ydelta);
  derr = Ydelta.NormF()/Y1.NormF();
  printf("Relative difference between A*B and K*B: %g\n", derr);

  // Inversion
  printf("Inversion ----------\n");
  CMatrix invA;
  A.Invert(invA);
  DMatrix DinvA, DD, Eye(N), Delta;
  invA.ConvertToDMatrix(DinvA);
  DA.MatMat(DinvA, DD, NORMAL, NORMAL);
  Eye.SetIdentity();
  DD.Subtract(Eye, Delta);
  derr = Delta.NormF()/sqrt((double)N);
  printf("Relative difference between A*invA and I: %g\n", derr);

  // Determinant
  printf("Determinant ----------\n");
  LogDet mLogDet = A.Det();
  LogDet mLogDet2 = DA.Det(GENERAL, UNFACT);
  if (mLogDet.Sign != mLogDet2.Sign) {
    printf("Error: Determinant has the wrong sign\n");
  }
  derr = fabs(mLogDet.LogAbsDet-mLogDet2.LogAbsDet)/fabs(mLogDet2.LogAbsDet);
  printf("Discrepancy in the log(|det|): %g\n", derr);

  // MatVecImplicit
  printf("MatVecImplicit ----------\n");
  DVector w(N), z, z2, z3;
  w.SetStandardNormal();
  A.MatVecImplicit<IsotropicGaussian, DPoint, DPointArray>
    (X, X, mKernel, w, z);

  DMatrix DA2;
  DA.SubtractDiagonal(lambda, DA2);
  DA2.MatVec(w, z2, TRANSPOSE);
  z.Subtract(z2, z3);
  derr = z3.Norm2()/z2.Norm2();
  printf("Discprepancy in return vector: %g\n", derr);

  // MatMatImplicit
  printf("MatMatImplicit ----------\n");
  DMatrix W(N, M), Z, Z2, Z3;
  W.SetStandardNormal();
  A.MatMatImplicit<IsotropicGaussian, DPoint, DPointArray>
    (X, X, mKernel, W, Z);

  DA2.MatMat(W, Z2, TRANSPOSE, NORMAL);
  Z.Subtract(Z2, Z3);
  derr = Z3.NormF()/Z2.NormF();
  printf("Discprepancy in return matrix: %g\n", derr);

  // BilinearImplicit
  printf("BilinearImplicit ----------\n");
  A.BilinearImplicit<IsotropicGaussian, DPoint, DPointArray>
    (X, X, mKernel, invA, z);

  z2.Init(N);
  DinvA.MatMat(DA2, W, NORMAL, NORMAL);
  for (INTEGER i = 0; i < N; i++) {
    W.GetColumn(i, w);
    DA2.GetColumn(i, y1);
    double z0 = y1.InProd(w);
    z2.SetEntry(i, z0);
  }

  z.Subtract(z2, z3);
  derr = z3.Norm2()/z2.Norm2();
  printf("Discprepancy in return vector: %g\n", derr);

  // Clean up
  Delete_1D_Array<INTEGER>(&Perm);
  Delete_1D_Array<INTEGER>(&iPerm);
  Delete_1D_Array<double>(&X_perm_back);


  // Part 1.5: Unsymmetric matrix ---------------------------------------------
  printf("\nPart 1.5: Unsymmetric matrix\n\n");

  // Sqrt
  printf("Sqrt ----------\n");
  CMatrix G;
  A.Sqrt(G);
  DMatrix DG, GG;
  G.ConvertToDMatrix(DG);
  DG.MatMat(DG, GG, NORMAL, TRANSPOSE);
  GG.Subtract(DA, Delta);
  derr = Delta.NormF()/DA.NormF();
  printf("Relative difference between GG' and A: %g\n", derr);

  // Matrix-vector multiplication
  printf("Matrix-vector multiplication ----------\n");
  printf("[G: RLCM; DG: dense form of G]\n");
  G.MatVec(b, y1, NORMAL);
  DG.MatVec(b, y3, NORMAL);

  y1.Subtract(y3, ydelta);
  derr = ydelta.Norm2()/y1.Norm2();
  printf("Discprepancy between G*b and DG*b: %g\n", derr);

  // Matrix-matrix multiplication
  printf("Matrix-matrix multiplication ----------\n");
  printf("[G: RLCM; DG: dense form of G]\n");
  G.MatMat(B, Y1, NORMAL);
  DG.MatMat(B, Y3, NORMAL, NORMAL);

  Y1.Subtract(Y3, Ydelta);
  derr = Ydelta.NormF()/Y1.NormF();
  printf("Discprepancy between G*B and DG*B: %g\n", derr);

  // Inversion
  printf("Inversion ----------\n");
  CMatrix invG;
  G.Invert(invG);
  DMatrix DinvG;
  invG.ConvertToDMatrix(DinvG);
  DG.MatMat(DinvG, DD, NORMAL, NORMAL);
  DD.Subtract(Eye, Delta);
  derr = Delta.NormF()/sqrt((double)N);
  printf("Relative difference between G*invG and I: %g\n", derr);

  // Determinant
  printf("Determinant ----------\n");
  mLogDet = G.Det();
  mLogDet2 = DG.Det(GENERAL, UNFACT);
  if (mLogDet.Sign != mLogDet2.Sign) {
    printf("Error: Determinant has the wrong sign\n");
  }
  derr = fabs(mLogDet.LogAbsDet-mLogDet2.LogAbsDet)/fabs(mLogDet2.LogAbsDet);
  printf("Discrepancy in the log(|det|): %g\n", derr);

  return 0;

}


//-----------------------------------------------------------------------------
void PrintPivots(const Node *mNode) {

  if (mNode->NumChild == 0) {
    return;
  }

  printf("Bounding box: ");
  for (INTEGER i = 0; i < mNode->Dim; i++) {
    printf("[%g, %g]", mNode->BBox[i], mNode->BBox[i + mNode->Dim]);
    if (i != mNode->Dim-1) {
      printf(" x ");
    }
    else {
      printf("\n");
    }
  }

  mNode->P.PrintPointArray("Pivots");

  PrintPivots(mNode->LeftChild);
  PrintPivots(mNode->LeftChild->RightSibling);

}


//-----------------------------------------------------------------------------
void ApproximationError(const Node *mNode, const DMatrix &DA, const DMatrix &B) {

  DMatrix DA_sub, B_sub, Delta;
  double derr;

  // Leaf node
  if (mNode->NumChild == 0) {

    DA.GetBlock(mNode->start, mNode->n, mNode->start, mNode->n, DA_sub);
    B.GetBlock(mNode->start, mNode->n, mNode->start, mNode->n, B_sub);
    DA_sub.Subtract(B_sub, Delta);
    derr = Delta.NormF()/B_sub.NormF();
    printf("Relative difference in Frobenius norm A[%ld:%ld, %ld:%ld]: %g\n",
           (long)mNode->start, (long)mNode->start + mNode->n - 1,
           (long)mNode->start, (long)mNode->start + mNode->n - 1, derr);

    return;
  }

  // Internal node
  const Node *left = mNode->LeftChild;
  const Node *right = mNode->LeftChild->RightSibling;
  ApproximationError(left, DA, B);
  ApproximationError(right, DA, B);

  DA.GetBlock(left->start, left->n, right->start, right->n, DA_sub);
  B.GetBlock(left->start, left->n, right->start, right->n, B_sub);
  DA_sub.Subtract(B_sub, Delta);
  derr = Delta.NormF()/B_sub.NormF();
  printf("Relative difference in Frobenius norm A[%ld:%ld, %ld:%ld]: %g\n",
         (long)left->start, (long)left->start + left->n - 1,
         (long)right->start, (long)right->start + right->n - 1, derr);

}
