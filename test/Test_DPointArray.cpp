// Unit test of the class DPointArray.

#include "LibCMatrix.hpp"

void Test_Bipartition(const DPointArray &X, const double *mX,
                      const char *FuncName, INTEGER start, INTEGER n,
                      const INTEGER *perm, const DPoint &normal, double offset,
                      INTEGER m1);

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
  printf("|  Test_DPointArray.cpp                                         |\n");
  printf("+---------------------------------------------------------------+\n");

  // Test GetPoint() ----------------------------------------------------------
  //
  //     [  1.0  0.0  0.3  0.4 ]
  // X = [  0.0  0.0  1.0  0.2 ]
  //     [ -0.1  1.0  0.0  0.0 ]
  //
  // X(1,:) = [  0.0  0.0  1.0  0.2 ]
  INTEGER N = 3;
  INTEGER d = 4;
  double mX[] = { 1.0, 0.0, -0.1, 0.0, 0.0, 1.0, 0.3, 1.0, 0.0, 0.4, 0.2, 0.0 };
  DPointArray X(N, d);
  memcpy(X.GetPointer(), mX, N*d*sizeof(double));
  DPoint x;
  X.GetPoint(1, x);
  double x_truth[] = { 0.0, 0.0, 1.0, 0.2 };
  double derr = Diff1(x.GetPointer(), x_truth, d);
  printf("Test GetPoint() : Discrepancy in return vector %g\n", derr);

  // Test GetSubset() ---------------------------------------------------------
  //
  // Reuse X
  //
  // X(1:2,:) = [  0.0  0.0  1.0  0.2 ]
  //            [ -0.1  1.0  0.0  0.0 ]
  INTEGER start = 1;
  INTEGER n = 2;
  DPointArray Y;
  X.GetSubset(start, n, Y);
  double Y_truth[] = { 0.0, -0.1, 0.0, 1.0, 1.0, 0.0, 0.2, 0.0 };
  derr = Diff1(Y.GetPointer(), Y_truth, n*d);
  printf("Test GetSubset() : Discrepancy in return matrix %g\n", derr);

  // Test GetSubset() ---------------------------------------------------------
  //
  // Reuse X
  //
  // X([0 2],:) = [  1.0  0.0  0.3  0.4 ]
  //              [ -0.1  1.0  0.0  0.0 ]
  INTEGER idx[] = { 0, 2 };
  X.GetSubset(idx, n, Y);
  double Y_truth2[] = { 1.0, -0.1, 0.0, 1.0, 0.3, 0.0, 0.4, 0.0 };
  derr = Diff1(Y.GetPointer(), Y_truth2, n*d);
  printf("Test GetSubset() : Discrepancy in return matrix %g\n", derr);

  // Test AsDMatrix() ---------------------------------------------------------
  //
  //           [  1.0  0.0  0.3  0.4 ]
  // Reuse X = [  0.0  0.0  1.0  0.2 ]
  //           [ -0.1  1.0  0.0  0.0 ]
  DMatrix A;
  X.AsDMatrix(A);
  double A_truth[] = { 1.0, 0.0, -0.1, 0.0, 0.0, 1.0, 0.3, 1.0, 0.0, 0.4, 0.2, 0.0 };
  derr = Diff1(X.GetPointer(), A_truth, N*d);
  printf("Test AsDMatrix() : Discrepancy in return matrix %g\n", derr);

  // Test Center() ------------------------------------------------------------
  //
  // Reuse X
  //
  // c = [  0.3  0.3333333333333333  0.4333333333333333  0.2 ]
  DPoint c;
  X.Center(c);
  double c_truth[] = { 0.3, 0.3333333333333333, 0.4333333333333333, 0.2 };
  derr = Diff1(c.GetPointer(), c_truth, d);
  printf("Test Center() : Discrepancy in return vector %g\n", derr);

  // Test StdvDist() ----------------------------------------------------------
  //
  // Reuse X
  //
  // stdv = 0.819213715162967
  double stdv = X.StdvDist();
  double stdv_truth = 0.819213715162967;
  derr = fabs(stdv-stdv_truth);
  printf("Test StdvDist() : Discrepancy in return value %g\n", derr);

  // Test MaxPointNorm2() -----------------------------------------------------
  //
  // Reuse X
  //
  // maxp = 1.25
  double maxp = X.MaxPointNorm2();
  double maxp_truth = 1.25;
  derr = fabs(maxp-maxp_truth);
  printf("Test MaxPointNorm2() : Discrepancy in return value %g\n", derr);

  // Test MatVec() ------------------------------------------------------------
  //
  // [  0.29 ]   [  1.0  0.0  0.3  0.4 ]   [  0.3 ]
  // [  0.08 ] = [  0.0  0.0  1.0  0.2 ] * [  0.0 ]
  // [ -0.03 ]   [ -0.1  1.0  0.0  0.0 ]   [  0.1 ]
  //                                       [ -0.1 ]
  // ModeX = NORMAL
  double mb[] = { 0.3, 0.0, 0.1, -0.1 };
  DVector Vb(d);
  memcpy(Vb.GetPointer(), mb, d*sizeof(double));
  DVector y;
  X.MatVec(Vb, y, NORMAL);
  double y_truth[] = { 0.29, 0.08, -0.03 };
  derr = Diff1(y.GetPointer(), y_truth, N);
  printf("Test MatVec() : Discrepancy in return vector %g\n", derr);

  DPoint Pb(d);
  memcpy(Pb.GetPointer(), mb, d*sizeof(double));
  X.MatVec(Pb, y, NORMAL);
  derr = Diff1(y.GetPointer(), y_truth, N);
  printf("Test MatVec() : Discrepancy in return vector %g\n", derr);

  // ModeX = TRANSPOSE
  double mXT[] = { 1.0, 0.0, 0.3, 0.4, 0.0, 0.0, 1.0, 0.2, -0.1, 1.0, 0.0, 0.0 };
  DPointArray XT(d, N);
  memcpy(XT.GetPointer(), mXT, N*d*sizeof(double));
  XT.MatVec(Vb, y, TRANSPOSE);
  derr = Diff1(y.GetPointer(), y_truth, N);
  printf("Test MatVec() : Discrepancy in return vector %g\n", derr);

  XT.MatVec(Pb, y, TRANSPOSE);
  derr = Diff1(y.GetPointer(), y_truth, N);
  printf("Test MatVec() : Discrepancy in return vector %g\n", derr);

  // Test MatMat() ------------------------------------------------------------
  //
  // [  0.03   0.60 ]   [  1.0  0.0  0.3  0.4 ]   [  0.0   1.0 ]
  // [  0.10  -0.20 ] = [  0.0  0.0  1.0  0.2 ] * [ -0.4   0.3 ]
  // [ -0.40   0.20 ]   [ -0.1  1.0  0.0  0.0 ]   [  0.1   0.0 ]
  //                                              [  0.0  -1.0 ]
  // ModeX = NORMAL, ModeB = NORMAL
  INTEGER MB = 4, NB = 2;
  double mB[] = { 0.0, -0.4, 0.1, 0.0, 1.0, 0.3, 0.0, -1.0 };
  DMatrix MatB(MB, NB);
  memcpy(MatB.GetPointer(), mB, MB*NB*sizeof(double));
  DMatrix MatY;
  X.MatMat(MatB, MatY, NORMAL, NORMAL);
  double MatY_truth[] = { 0.03, 0.10, -0.40, 0.60, -0.20, 0.20 };
  derr = Diff1(MatY.GetPointer(), MatY_truth, N*NB);
  printf("Test MatMat() : Discrepancy in return matrix %g\n", derr);

  DPointArray PtArryB(MB, NB);
  memcpy(PtArryB.GetPointer(), mB, MB*NB*sizeof(double));
  X.MatMat(PtArryB, MatY, NORMAL, NORMAL);
  derr = Diff1(MatY.GetPointer(), MatY_truth, N*NB);
  printf("Test MatMat() : Discrepancy in return matrix %g\n", derr);

  // ModeX = NORMAL, ModeB = TRANSPOSE
  double mBT[] = { 0.0, 1.0, -0.4, 0.3, 0.1, 0.0, 0.0, -1.0 };
  DMatrix MatBT(NB, MB);
  memcpy(MatBT.GetPointer(), mBT, MB*NB*sizeof(double));
  X.MatMat(MatBT, MatY, NORMAL, TRANSPOSE);
  derr = Diff1(MatY.GetPointer(), MatY_truth, N*NB);
  printf("Test MatMat() : Discrepancy in return matrix %g\n", derr);

  DPointArray PtArryBT(NB, MB);
  memcpy(PtArryBT.GetPointer(), mBT, MB*NB*sizeof(double));
  X.MatMat(PtArryBT, MatY, NORMAL, TRANSPOSE);
  derr = Diff1(MatY.GetPointer(), MatY_truth, N*NB);
  printf("Test MatMat() : Discrepancy in return matrix %g\n", derr);

  // ModeX = TRANSPOSE, ModeB = NORMAL
  XT.MatMat(MatB, MatY, TRANSPOSE, NORMAL);
  derr = Diff1(MatY.GetPointer(), MatY_truth, N*NB);
  printf("Test MatMat() : Discrepancy in return matrix %g\n", derr);

  XT.MatMat(PtArryB, MatY, TRANSPOSE, NORMAL);
  derr = Diff1(MatY.GetPointer(), MatY_truth, N*NB);
  printf("Test MatMat() : Discrepancy in return matrix %g\n", derr);

  // ModeX = TRANSPOSE, ModeB = TRANSPOSE
  XT.MatMat(MatBT, MatY, TRANSPOSE, TRANSPOSE);
  derr = Diff1(MatY.GetPointer(), MatY_truth, N*NB);
  printf("Test MatMat() : Discrepancy in return matrix %g\n", derr);

  XT.MatMat(PtArryBT, MatY, TRANSPOSE, TRANSPOSE);
  derr = Diff1(MatY.GetPointer(), MatY_truth, N*NB);
  printf("Test MatMat() : Discrepancy in return matrix %g\n", derr);

  // Test RandomBipartition() -------------------------------------------------
  //
  //     [ 1.0   0.0  -0.1 ]
  //     [ 0.0   0.0   1.0 ] ---+
  //     [ 0.3   1.0   0.0 ]    |
  // X = [ 0.4   0.2   0.0 ]    +--- partition this part
  //     [ 1.2   0.9  -0.1 ]    |
  //     [ 0.4  -1.0   0.0 ] ---+
  //     [ 0.0   0.0   0.7 ]
  N = 7;
  d = 3;
  double mX2[] = {  1.0, 0.0, 0.3, 0.4,  1.2,  0.4, 0.0,
                    0.0, 0.0, 1.0, 0.2,  0.9, -1.0, 0.0,
                   -0.1, 1.0, 0.0, 0.0, -0.1,  0.0, 0.7 };
  X.Init(N, d);
  memcpy(X.GetPointer(), mX2, N*d*sizeof(double));
  start = 1;
  n = 5;
  INTEGER N0 = 2;
  INTEGER perm[] = { 0, 1, 2, 3, 4, 5, 6 };
  DPoint normal;
  double offset;
  INTEGER m1 = X.RandomBipartition(start, n, N0, perm, normal, offset);

  const char *FuncName = "RandomBipartition";
  Test_Bipartition(X, mX2, FuncName, start, n, perm, normal, offset, m1);

  // Test PCABipartition() ----------------------------------------------------
  //
  // Reuse above X
  memcpy(X.GetPointer(), mX2, N*d*sizeof(double));
  for (INTEGER i = 0; i < N; i++) {
    perm[i] = i;
  }
  m1 = X.PCABipartition(start, n, N0, perm, normal, offset);

  // Check that the normal is correct
  // normal = [ 0.304486005685167; 0.928148426475952; -0.214076553531924 ],
  // up to sign
  double normal_truth[] = { 0.304486005685167, 0.928148426475952, -0.214076553531924 };
  double *normal_ptr = normal.GetPointer();
  if (normal_ptr[0]*normal_truth[0] < 0) {
    for (INTEGER i = 0; i < d; i++) {
      normal_truth[i] *= -1.0;
    }
  }
  derr = Diff1(normal_ptr, normal_truth, d);
  printf("Test PCABipartition() : Part 0, Discrepancy in normal %g\n", derr);

  // Other checks
  const char *FuncName2 = "PCABipartition";
  Test_Bipartition(X, mX2, FuncName2, start, n, perm, normal, offset, m1);

  // Test RandomBipartition() -------------------------------------------------
  //
  // Another X
  //
  //     [ 1.0   0.0  -0.1   1.0 ]
  //     [ 0.0   0.0   1.0  -0.2 ] ---+
  // X = [ 0.4   0.2   0.0   0.7 ]    +--- partition this part
  //     [ 0.4  -1.0   0.0   0.0 ] ---+
  //     [ 0.0   0.0   0.7   0.0 ]
  N = 5;
  d = 4;
  double mX3[] = { 1.0,  0.0, 0.4,  0.4, 0.0,
                   0.0,  0.0, 0.2, -1.0, 0.0,
                  -0.1,  1.0, 0.0,  0.0, 0.7,
                   1.0, -0.2, 0.7,  0.0, 0.0};
  X.Init(N, d);
  memcpy(X.GetPointer(), mX3, N*d*sizeof(double));
  start = 1;
  n = 3;
  N0 = 1;
  for (INTEGER i = 0; i < N; i++) {
    perm[i] = i;
  }
  m1 = X.RandomBipartition(start, n, N0, perm, normal, offset);

  const char *FuncName3 = "RandomBipartition";
  Test_Bipartition(X, mX3, FuncName3, start, n, perm, normal, offset, m1);

  // Test PCABipartition() ----------------------------------------------------
  //
  // Reuse above X
  memcpy(X.GetPointer(), mX3, N*d*sizeof(double));
  for (INTEGER i = 0; i < N; i++) {
    perm[i] = i;
  }
  m1 = X.PCABipartition(start, n, N0, perm, normal, offset);

  // Check that the normal is correct
  // normal = [ -0.289048928337944;
  //             0.586654813940892;
  //             0.722622320844860;
  //            -0.223838843196287 ],
  // up to sign
  double normal_truth2[] = { -0.289048928337944,
                              0.586654813940892,
                              0.722622320844860,
                             -0.223838843196287 };
  normal_ptr = normal.GetPointer();
  if (normal_ptr[0]*normal_truth2[0] < 0) {
    for (INTEGER i = 0; i < d; i++) {
      normal_truth2[i] *= -1.0;
    }
  }
  derr = Diff1(normal_ptr, normal_truth2, d);
  printf("Test PCABipartition() : Part 0, Discrepancy in normal %g\n", derr);

  // Other checks
  const char *FuncName4 = "PCABipartition";
  Test_Bipartition(X, mX3, FuncName4, start, n, perm, normal, offset, m1);

  // Test ComputeBBox() and BBoxBipartition() ---------------------------------
  //
  //     [ 1.0   0.0  -0.1 ]
  //     [ 0.0  -0.9   1.0 ] ---+
  //     [ 0.3   1.0   0.0 ]    |
  // X = [ 0.4   0.2   0.0 ]    +--- partition this part
  //     [ 1.2   0.9  -0.1 ]    |
  //     [ 0.4  -1.0   0.0 ] ---+
  //     [ 0.0   0.0   0.7 ]
  N = 7;
  d = 3;
  double mX4[] = {  1.0,  0.0, 0.3, 0.4,  1.2,  0.4, 0.0,
                    0.0, -0.9, 1.0, 0.2,  0.9, -1.0, 0.0,
                   -0.1,  1.0, 0.0, 0.0, -0.1,  0.0, 0.7 };
  X.Init(N, d);
  memcpy(X.GetPointer(), mX4, N*d*sizeof(double));
  start = 1;
  n = 5;
  DPointArray Xsub;
  X.GetSubset(start, n, Xsub);

  // Check that the bounding box (of Xsub) is correct
  // bbox = [0.0; -1.0; -0.1; 1.2; 1.0; 1.0]
  double bbox_truth[] = {0.0, -1.0, -0.1, 1.2, 1.0, 1.0};
  double *bbox = NULL;
  INTEGER dim;
  Xsub.ComputeBBox(&bbox, dim);
  int ierr = abs((int)dim - (int)d);
  printf("Test ComputeBBox() : Discrepancy in return value %d\n", ierr);
  derr = Diff1(bbox, bbox_truth, d*2);
  printf("Test ComputeBBox() : Discrepancy in return vector %g\n", derr);

  // Run BBoxBipartition()
  N0 = 2;
  for (INTEGER i = 0; i < N; i++) {
    perm[i] = i;
  }
  INTEGER which_dim;
  m1 = X.BBoxBipartition(start, n, N0, perm, normal, offset, bbox, which_dim);
  Delete_1D_Array<double>(&bbox);

  // Check results of BBoxBipartition()
  const char *FuncName5 = "BBoxBipartition";
  Test_Bipartition(X, mX4, FuncName5, start, n, perm, normal, offset, m1);

  // Additionally, check normal ([0.0; 1.0; 0.0])
  normal_truth[0] = 0.0;
  normal_truth[1] = 1.0;
  normal_truth[2] = 0.0;
  derr = Diff1(normal.GetPointer(), normal_truth, d);
  printf("Test BBoxBipartition() : Part 4, Discrepancy in normal %g\n", derr);

  // Additionally, check offset (0.0)
  double offset_truth = 0.0;
  derr = fabs(offset-offset_truth);
  printf("Test BBoxBipartition() : Part 5, Discrepancy in offset %g\n", derr);

  // Additionally, check which_dim (1)
  INTEGER which_dim_truth = 1;
  ierr = abs((int)which_dim - (int)which_dim_truth);
  printf("Test BBoxBipartition() : Part 6, Error in which_dim %d\n", ierr);

  // Test SetRegularGrid() ----------------------------------------------------
  //
  // X = [ 1 3 6;
  //       2 3 6;
  //       1 4 6;
  //       2 4 6;
  //       1 5 6;
  //       2 5 6;
  //       1 3 7;
  //       2 3 7;
  //       1 4 7;
  //       2 4 7;
  //       1 5 7;
  //       2 5 7 ]
  INTEGER Dim[] = {2, 3, 2};
  double Lower[] = {1, 3, 6};
  double Upper[] = {2, 5, 7};
  X.SetRegularGrid(3, Dim, Lower, Upper);
  double mX5[] = { 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2,
                   3, 3, 4, 4, 5, 5, 3, 3, 4, 4, 5, 5,
                   6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7 };
  derr = Diff1(X.GetPointer(), mX5, 2*3*2*3);
  printf("Test SetRegularGrid() : Discrepancy in X %g\n", derr);

  return 0;

}


//-----------------------------------------------------------------------------
void Test_Bipartition(const DPointArray &X, const double *mX,
                      const char *FuncName, INTEGER start, INTEGER n,
                      const INTEGER *perm, const DPoint &normal, double offset,
                      INTEGER m1) {

  // Check that the permuted array contains the same content as the
  // original one, by using perm
  INTEGER N = X.GetN();
  INTEGER d = X.GetD();
  double *X_perm_back = NULL;
  New_1D_Array<double, INTEGER>(&X_perm_back, N*d);
  double *X_ptr = X.GetPointer();
  for (INTEGER i = 0; i < N; i++) {
    for (INTEGER j = 0; j < d; j++) {
      X_perm_back[perm[i]+j*N] = X_ptr[i+j*N];
    }
  }
  double derr = Diff1(X_perm_back, mX, N*d);
  printf("Test %s() : Part 1, Error in X %g\n", FuncName, derr);

  // Check that outside [start, start+n-1], perm[i] = i
  derr = 0.0;
  for (INTEGER i = 0; i < start; i++) {
    if (perm[i] != i) {
      derr += 1.0;
    }
  }
  for (INTEGER i = start+n; i < N; i++) {
    if (perm[i] != i) {
      derr += 1.0;
    }
  }
  printf("Test %s() : Part 2, Error in perm %g\n", FuncName, derr);

  // Check that the signed distances indeed have the correct sign
  DVector y;
  X.MatVec(normal, y, NORMAL);
  y.Subtract(offset);
  double *y_ptr = y.GetPointer();
  derr = 0.0;
  for (INTEGER i = start; i < start+m1; i++) {
    if (y_ptr[i] > 0) {
      derr += 1.0;
    }
  }
  for (INTEGER i = start+m1; i < start+n; i++) {
    if (y_ptr[i] < 0) {
      derr += 1.0;
    }
  }
  printf("Test %s() : Part 3, Error in sign %g\n", FuncName, derr);

  // Clean up
  Delete_1D_Array<double>(&X_perm_back);

}
