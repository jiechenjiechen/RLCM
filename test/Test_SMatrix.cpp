// Unit test of the class SMatrix.

#include "LibCMatrix.hpp"

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
  printf("|  Test_SMatrix.cpp                                             |\n");
  printf("+---------------------------------------------------------------+\n");

  // Test GetBlock() ----------------------------------------------------------
  //
  //     [  1.0  0.0  0.3  0.4 ]
  // A = [  0.0  0.0  1.0  0.2 ]
  //     [ -0.1  1.0  0.0  0.0 ]
  //
  // A(1:2,1:2) = [ 0.0  1.0 ]
  //              [ 1.0  0.0 ]
  INTEGER mA = 3;
  INTEGER nA = 4;
  INTEGER nnzA = 7;
  INTEGER startA[] = { 0, 3, 5, 7 };
  INTEGER idxA[] = { 0, 2, 3, 2, 3, 0, 1 };
  double SA[] = { 1.0, 0.3, 0.4, 1.0, 0.2, -0.1, 1.0 };
  SMatrix A(mA, nA, nnzA);
  memcpy(A.GetPointerStart(), startA, (mA+1)*sizeof(INTEGER));
  memcpy(A.GetPointerIdx(), idxA, nnzA*sizeof(INTEGER));
  memcpy(A.GetPointerA(), SA, nnzA*sizeof(double));
  INTEGER RowStart = 1;
  INTEGER nRow = 2;
  INTEGER ColStart = 1;
  INTEGER nCol = 2;
  DMatrix B;
  A.GetBlock(RowStart, nRow, ColStart, nCol, B);
  double B_truth[] = { 0.0, 1.0, 1.0, 0.0 };
  double derr = Diff1(B.GetPointer(), B_truth, 4);
  printf("Test GetBlock() : Discrepancy in return matrix %g\n", derr);

  // A([0 2],[1 3]) = [ 0.0  0.4 ]
  //                  [ 1.0  0.0 ]
  INTEGER IdxRow[] = {0, 2};
  INTEGER IdxCol[] = {1, 3};
  A.GetBlock(IdxRow, nRow, IdxCol, nCol, B);
  double B_truth2[] = { 0.0, 1.0, 0.4, 0.0 };
  derr = Diff1(B.GetPointer(), B_truth2, 4);
  printf("Test GetBlock() : Discrepancy in return matrix %g\n", derr);

  // A([0 2],1:2) = [ 0.0  0.3 ]
  //                [ 1.0  0.0 ]
  A.GetBlock(IdxRow, nRow, ColStart, nCol, B);
  double B_truth3[] = { 0.0, 1.0, 0.3, 0.0 };
  derr = Diff1(B.GetPointer(), B_truth3, 4);
  printf("Test GetBlock() : Discrepancy in return matrix %g\n", derr);

  // Test MatVec() ------------------------------------------------------------
  //
  // [  0.29 ]   [  1.0  0.0  0.3  0.4 ]   [  0.3 ]
  // [  0.08 ] = [  0.0  0.0  1.0  0.2 ] * [  0.0 ]
  // [ -0.03 ]   [ -0.1  1.0  0.0  0.0 ]   [  0.1 ]
  //                                       [ -0.1 ]
  // ModeA = NORMAL
  double mb[] = { 0.3, 0.0, 0.1, -0.1 };
  DVector b(nA);
  memcpy(b.GetPointer(), mb, nA*sizeof(double));
  DVector y;
  A.MatVec(b, y, NORMAL);
  double my_truth[] = { 0.29, 0.08, -0.03 };
  derr = Diff1(y.GetPointer(), my_truth, mA);
  printf("Test MatVec() : Discrepancy in return vector %g\n", derr);

  // ModeA = TRANSPOSE
  INTEGER startAT[] = { 0, 2, 3, 5, 7 };
  INTEGER idxAT[] = { 0, 2, 2, 0, 1, 0, 1 };
  double SAT[] = { 1.0, -0.1, 1.0, 0.3, 1.0, 0.4, 0.2 };
  SMatrix AT(nA, mA, nnzA);
  memcpy(AT.GetPointerStart(), startAT, (nA+1)*sizeof(INTEGER));
  memcpy(AT.GetPointerIdx(), idxAT, nnzA*sizeof(INTEGER));
  memcpy(AT.GetPointerA(), SAT, nnzA*sizeof(double));
  AT.MatVec(b, y, TRANSPOSE);
  derr = Diff1(y.GetPointer(), my_truth, mA);
  printf("Test MatVec() : Discrepancy in return vector %g\n", derr);

  // Test MatMat() ------------------------------------------------------------
  //
  // [  0.03   0.60 ]   [  1.0  0.0  0.3  0.4 ]   [  0.0   1.0 ]
  // [  0.10  -0.20 ] = [  0.0  0.0  1.0  0.2 ] * [ -0.4   0.3 ]
  // [ -0.40   0.20 ]   [ -0.1  1.0  0.0  0.0 ]   [  0.1   0.0 ]
  //                                              [  0.0  -1.0 ]
  INTEGER MB2 = 4, NB2 = 2;
  double mB2[] = { 0.0, -0.4, 0.1, 0.0, 1.0, 0.3, 0.0, -1.0 };
  DMatrix B2(MB2, NB2);
  memcpy(B2.GetPointer(), mB2, MB2*NB2*sizeof(double));
  DMatrix C;
  A.MatMat(B2, C, NORMAL, NORMAL);
  double matmat_truth[] = { 0.03, 0.10, -0.40, 0.60, -0.20, 0.20 };
  derr = Diff1(C.GetPointer(), matmat_truth, mA*NB2);
  printf("Test MatMat() : Discrepancy in return matrix %g\n", derr);

  double mB2T[] = { 0.0, 1.0, -0.4, 0.3, 0.1, 0.0, 0.0, -1.0 };
  DMatrix B2T(NB2, MB2);
  memcpy(B2T.GetPointer(), mB2T, MB2*NB2*sizeof(double));
  A.MatMat(B2T, C, NORMAL, TRANSPOSE);
  derr = Diff1(C.GetPointer(), matmat_truth, mA*NB2);
  printf("Test MatMat() : Discrepancy in return matrix %g\n", derr);

  AT.MatMat(B2, C, TRANSPOSE, NORMAL);
  derr = Diff1(C.GetPointer(), matmat_truth, mA*NB2);
  printf("Test MatMat() : Discrepancy in return matrix %g\n", derr);

  AT.MatMat(B2T, C, TRANSPOSE, TRANSPOSE);
  derr = Diff1(C.GetPointer(), matmat_truth, mA*NB2);
  printf("Test MatMat() : Discrepancy in return matrix %g\n", derr);

  // Test Diag() --------------------------------------------------------------
  //
  //     [  1.0  0.0  0.3 ]
  // A = [  0.0  0.0  1.0 ]
  //     [ -0.1  1.0  0.0 ]
  INTEGER nA2 = 3;
  INTEGER nnzA2 = 5;
  INTEGER startA2[] = { 0, 2, 3, 5 };
  INTEGER idxA2[] = { 0, 2, 2, 0, 1 };
  double SA2[] = { 1.0, 0.3, 1.0, -0.1, 1.0 };
  SMatrix A2(nA2, nnzA2);
  memcpy(A2.GetPointerStart(), startA2, (nA2+1)*sizeof(INTEGER));
  memcpy(A2.GetPointerIdx(), idxA2, nnzA2*sizeof(INTEGER));
  memcpy(A2.GetPointerA(), SA2, nnzA2*sizeof(double));
  double b2[3];
  A2.Diag(b2);
  double b2_truth[] = { 1.0, 0.0, 0.0 };
  derr = Diff1(b2, b2_truth, nA2);
  printf("Test Diag() : Discrepancy in return vector %g\n", derr);

  A2.Diag(b);
  derr = Diff1(b.GetPointer(), b2_truth, nA2);
  printf("Test Diag() : Discrepancy in return vector %g\n", derr);

  // Test Divide() ------------------------------------------------------------
  //
  // Reuse above A2. b = 2
  memcpy(A2.GetPointerA(), SA2, nnzA2*sizeof(double));
  double b3 = 2.0;
  A2.Divide(b3);
  double divide_truth[] = { 0.5, 0.15, 0.5, -0.05, 0.5 };
  derr = Diff1(A2.GetPointerA(), divide_truth, nnzA2);
  printf("Test Divide() : Discrepancy in return matrix %g\n", derr);

  memcpy(A2.GetPointerA(), SA2, nnzA2*sizeof(double));
  SMatrix A3;
  A2.Divide(b3, A3);
  derr = Diff1(A3.GetPointerA(), divide_truth, nnzA2);
  printf("Test Divide() : Discrepancy in return matrix %g\n", derr);

  // Test SymmetricDivide() ---------------------------------------------------
  //
  // Reuse above A2
  //
  // [  1.0  0.0  0.3 ]          [ 2.0 ]
  // [  0.0  0.0  1.0 ] , symdiv [ 3.0 ].
  // [ -0.1  1.0  0.0 ]          [ 4.0 ]
  //
  //          [  0.25    0.0                     0.0375 ]
  // result = [  0.0     0.0                     0.08333333333333333333 ]
  //          [ -0.0125  0.08333333333333333333  0.0 ]
  memcpy(A2.GetPointerA(), SA2, nnzA2*sizeof(double));
  double mb4[] = { 2.0, 3.0, 4.0 };
  DVector b4(nA2);
  memcpy(b4.GetPointer(), mb4, nA2*sizeof(double));
  A2.SymmetricDivide(b4);
  double sym_divide_truth[] = { 0.25, 0.0375, 0.08333333333333333333, -0.0125, 0.08333333333333333333 };
  derr = Diff1(A2.GetPointerA(), sym_divide_truth, nnzA2);
  printf("Test SymmetricDivide() : Discrepancy in return matrix %g\n", derr);

  memcpy(A2.GetPointerA(), SA2, nnzA2*sizeof(double));
  A2.SymmetricDivide(b4, A3);
  derr = Diff1(A3.GetPointerA(), sym_divide_truth, nnzA2);
  printf("Test SymmetricDivide() : Discrepancy in return matrix %g\n", derr);

  return 0;

}
