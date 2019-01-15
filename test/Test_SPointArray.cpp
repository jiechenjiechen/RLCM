// Unit test of the class SPointArray.

#include "LibCMatrix.hpp"

void Test_Bipartition(const SPointArray &X, const INTEGER *start,
                      const INTEGER *idx, const double *val,
                      const char *FuncName, INTEGER istart, INTEGER n,
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
  printf("|  Test_SPointArray.cpp                                         |\n");
  printf("+---------------------------------------------------------------+\n");

  // Test assignment operator (dense to sparse) -------------------------------
  //
  //     [  1.0  0.0  0.3  0.4 ]
  // X = [  0.0  0.0  1.0  0.2 ]
  //     [ -0.1  1.0  0.0  0.0 ]
  INTEGER N = 3;
  INTEGER d = 4;
  INTEGER nnz = 7;
  INTEGER start[] = { 0, 3, 5, 7 };
  INTEGER idx[] = { 0, 2, 3, 2, 3, 0, 1 };
  double val[] = { 1.0, 0.3, 0.4, 1.0, 0.2, -0.1, 1.0 };
  SPointArray X(N, d, nnz);
  memcpy(X.GetPointerStart(), start, (N+1)*sizeof(INTEGER));
  memcpy(X.GetPointerIdx(), idx, nnz*sizeof(INTEGER));
  memcpy(X.GetPointerX(), val, nnz*sizeof(double));

  double mDX[] = {1.0, 0.0, -0.1, 0.0, 0.0, 1.0, 0.3, 1.0, 0.0, 0.4, 0.2, 0.0};
  DPointArray DX(N, d);
  memcpy(DX.GetPointer(), mDX, N*d*sizeof(double));
  SPointArray SX;
  SX = DX;
  double derr = Diff1(SX.GetPointerStart(), X.GetPointerStart(), X.GetN()+1);
  printf("Test operator= (dense to sparse) : Discrepancy in start array %g\n", derr);
  derr = Diff1(SX.GetPointerIdx(), X.GetPointerIdx(), X.GetNNZ());
  printf("Test operator= (dense to sparse) : Discrepancy in idx array %g\n", derr);
  derr = Diff1(SX.GetPointerX(), X.GetPointerX(), X.GetNNZ());
  printf("Test operator= (dense to sparse) : Discrepancy in val array %g\n", derr);

  // Test GetPoint() ----------------------------------------------------------
  //
  // Reuse X
  //
  // X(1,:) = [  0.0  0.0  1.0  0.2 ]
  SPoint Sx;
  X.GetPoint(1, Sx);
  INTEGER Sx_idx_truth[] = { 2, 3 };
  double Sx_val_truth[] = { 1.0, 0.2 };
  derr = Diff1(Sx.GetPointerIdx(), Sx_idx_truth, Sx.GetNNZ());
  printf("Test GetPoint() : Discrepancy in idx array %g\n", derr);
  derr = Diff1(Sx.GetPointerX(), Sx_val_truth, Sx.GetNNZ());
  printf("Test GetPoint() : Discrepancy in val array %g\n", derr);

  DPoint Dx;
  X.GetPoint(1, Dx);
  double Dx_truth[] = { 0.0, 0.0, 1.0, 0.2 };
  derr = Diff1(Dx.GetPointer(), Dx_truth, d);
  printf("Test GetPoint() : Discrepancy in return vector %g\n", derr);

  // Test GetSubset() ---------------------------------------------------------
  //
  // Reuse X
  //
  // X(1:2,:) = [  0.0  0.0  1.0  0.2 ]
  //            [ -0.1  1.0  0.0  0.0 ]
  INTEGER istart = 1;
  INTEGER n = 2;
  SPointArray SY;
  X.GetSubset(istart, n, SY);
  INTEGER SY_start_truth[] = { 0, 2, 4 };
  INTEGER SY_idx_truth[] = { 2, 3, 0, 1 };
  double SY_val_truth[] = { 1.0, 0.2, -0.1, 1.0 };
  derr = Diff1(SY.GetPointerStart(), SY_start_truth, n+1);
  printf("Test GetSubset() : Discrepancy in start array %g\n", derr);
  derr = Diff1(SY.GetPointerIdx(), SY_idx_truth, SY.GetNNZ());
  printf("Test GetSubset() : Discrepancy in idx array %g\n", derr);
  derr = Diff1(SY.GetPointerX(), SY_val_truth, SY.GetNNZ());
  printf("Test GetSubset() : Discrepancy in val array %g\n", derr);

  DPointArray DY;
  X.GetSubset(istart, n, DY);
  double DY_truth[] = { 0.0, -0.1, 0.0, 1.0, 1.0, 0.0, 0.2, 0.0 };
  derr = Diff1(DY.GetPointer(), DY_truth, n*d);
  printf("Test GetSubset() : Discrepancy in return matrix %g\n", derr);

  // Test GetSubset() ---------------------------------------------------------
  //
  // Reuse X
  //
  // X([0 2],:) = [  1.0  0.0  0.3  0.4 ]
  //              [ -0.1  1.0  0.0  0.0 ]
  INTEGER iidx[] = { 0, 2 };
  X.GetSubset(iidx, n, SY);
  INTEGER SY_start_truth2[] = { 0, 3, 5 };
  INTEGER SY_idx_truth2[] = { 0, 2, 3, 0, 1 };
  double SY_val_truth2[] = { 1.0, 0.3, 0.4, -0.1, 1.0 };
  derr = Diff1(SY.GetPointerStart(), SY_start_truth2, n+1);
  printf("Test GetSubset() : Discrepancy in start array %g\n", derr);
  derr = Diff1(SY.GetPointerIdx(), SY_idx_truth2, SY.GetNNZ());
  printf("Test GetSubset() : Discrepancy in idx array %g\n", derr);
  derr = Diff1(SY.GetPointerX(), SY_val_truth2, SY.GetNNZ());
  printf("Test GetSubset() : Discrepancy in val array %g\n", derr);

  X.GetSubset(iidx, n, DY);
  double DY_truth2[] = { 1.0, -0.1, 0.0, 1.0, 0.3, 0.0, 0.4, 0.0 };
  derr = Diff1(DY.GetPointer(), DY_truth2, n*d);
  printf("Test GetSubset() : Discrepancy in return matrix %g\n", derr);

  // Test AsDMatrix() ---------------------------------------------------------
  //
  //           [  1.0  0.0  0.3  0.4 ]
  // Reuse X = [  0.0  0.0  1.0  0.2 ]
  //           [ -0.1  1.0  0.0  0.0 ]
  DMatrix A;
  X.AsDMatrix(A);
  double A_truth[] = { 1.0, 0.0, -0.1, 0.0, 0.0, 1.0, 0.3, 1.0, 0.0, 0.4, 0.2, 0.0 };
  derr = Diff1(mDX, A_truth, N*d);
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
  INTEGER startT[] = { 0, 2, 3, 5, 7 };
  INTEGER idxT[] = { 0, 2, 2, 0, 1, 0, 1 };
  double valT[] = { 1.0, -0.1, 1.0, 0.3, 1.0, 0.4, 0.2 };
  SPointArray XT(d, N, nnz);
  memcpy(XT.GetPointerStart(), startT, (d+1)*sizeof(INTEGER));
  memcpy(XT.GetPointerIdx(), idxT, nnz*sizeof(INTEGER));
  memcpy(XT.GetPointerX(), valT, nnz*sizeof(double));
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

  INTEGER PtArryB_nnz = 5;
  INTEGER PtArryB_start[] = { 0, 1, 3, 4, 5 };
  INTEGER PtArryB_idx[] = { 1, 0, 1, 0, 1 };
  double PtArryB_val[] = { 1.0, -0.4, 0.3, 0.1, -1.0 };
  SPointArray PtArryB(MB, NB, PtArryB_nnz);
  memcpy(PtArryB.GetPointerStart(), PtArryB_start, (MB+1)*sizeof(INTEGER));
  memcpy(PtArryB.GetPointerIdx(), PtArryB_idx, PtArryB_nnz*sizeof(INTEGER));
  memcpy(PtArryB.GetPointerX(), PtArryB_val, PtArryB_nnz*sizeof(double));
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

  INTEGER PtArryBT_nnz = 5;
  INTEGER PtArryBT_start[] = { 0, 2, 5 };
  INTEGER PtArryBT_idx[] = { 1, 2, 0, 1, 3 };
  double PtArryBT_val[] = { -0.4, 0.1, 1.0, 0.3, -1.0 };
  SPointArray PtArryBT(NB, MB, PtArryBT_nnz);
  memcpy(PtArryBT.GetPointerStart(), PtArryBT_start, (NB+1)*sizeof(INTEGER));
  memcpy(PtArryBT.GetPointerIdx(), PtArryBT_idx, PtArryBT_nnz*sizeof(INTEGER));
  memcpy(PtArryBT.GetPointerX(), PtArryBT_val, PtArryBT_nnz*sizeof(double));
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
  nnz = 13;
  INTEGER start2[] = { 0, 2, 3, 5, 7, 10, 12, 13 };
  INTEGER idx2[] = { 0, 2, 2, 0, 1, 0, 1, 0, 1, 2, 0, 1, 2 };
  double val2[] = { 1.0, -0.1, 1.0, 0.3, 1.0, 0.4, 0.2, 1.2, 0.9, -0.1, 0.4, -1.0, 0.7 };
  X.Init(N, d, nnz);
  memcpy(X.GetPointerStart(), start2, (N+1)*sizeof(INTEGER));
  memcpy(X.GetPointerIdx(), idx2, nnz*sizeof(INTEGER));
  memcpy(X.GetPointerX(), val2, nnz*sizeof(double));
  istart = 1;
  n = 5;
  INTEGER N0 = 2;
  INTEGER perm[] = { 0, 1, 2, 3, 4, 5, 6 };
  DPoint normal;
  double offset;
  INTEGER m1 = X.RandomBipartition(istart, n, N0, perm, normal, offset);

  const char *FuncName = "RandomBipartition";
  Test_Bipartition(X, start2, idx2, val2, FuncName, istart, n, perm, normal, offset, m1);

  // Test PCABipartition() ----------------------------------------------------
  //
  // Reuse above X
  memcpy(X.GetPointerStart(), start2, (N+1)*sizeof(INTEGER));
  memcpy(X.GetPointerIdx(), idx2, nnz*sizeof(INTEGER));
  memcpy(X.GetPointerX(), val2, nnz*sizeof(double));
  for (INTEGER i = 0; i < N; i++) {
    perm[i] = i;
  }
  m1 = X.PCABipartition(istart, n, N0, perm, normal, offset);

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
  Test_Bipartition(X, start2, idx2, val2, FuncName2, istart, n, perm, normal, offset, m1);

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
  nnz = 11;
  INTEGER start3[] = { 0, 3, 5, 8, 10, 11 };
  INTEGER idx3[] = { 0, 2, 3, 2, 3, 0, 1, 3, 0, 1, 2 };
  double val3[] = { 1.0, -0.1, 1.0, 1.0, -0.2, 0.4, 0.2, 0.7, 0.4, -1.0, 0.7 };
  X.Init(N, d, nnz);
  memcpy(X.GetPointerStart(), start3, (N+1)*sizeof(INTEGER));
  memcpy(X.GetPointerIdx(), idx3, nnz*sizeof(INTEGER));
  memcpy(X.GetPointerX(), val3, nnz*sizeof(double));
  istart = 1;
  n = 3;
  N0 = 1;
  for (INTEGER i = 0; i < N; i++) {
    perm[i] = i;
  }
  m1 = X.RandomBipartition(istart, n, N0, perm, normal, offset);

  const char *FuncName3 = "RandomBipartition";
  Test_Bipartition(X, start3, idx3, val3, FuncName3, istart, n, perm, normal, offset, m1);

  // Test PCABipartition() ----------------------------------------------------
  //
  // Reuse above X
  memcpy(X.GetPointerStart(), start3, (N+1)*sizeof(INTEGER));
  memcpy(X.GetPointerIdx(), idx3, nnz*sizeof(INTEGER));
  memcpy(X.GetPointerX(), val3, nnz*sizeof(double));
  for (INTEGER i = 0; i < N; i++) {
    perm[i] = i;
  }
  m1 = X.PCABipartition(istart, n, N0, perm, normal, offset);

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
  Test_Bipartition(X, start3, idx3, val3, FuncName4, istart, n, perm, normal, offset, m1);

  return 0;

}


//-----------------------------------------------------------------------------
void Test_Bipartition(const SPointArray &X, const INTEGER *start,
                      const INTEGER *idx, const double *val,
                      const char *FuncName, INTEGER istart, INTEGER n,
                      const INTEGER *perm, const DPoint &normal, double offset,
                      INTEGER m1) {

  // Check that the permuted array contains the same content as the
  // original one, by using perm
  INTEGER N = X.GetN();
  INTEGER nnz = X.GetNNZ();
  INTEGER *X_perm_back_start = NULL;
  INTEGER *X_perm_back_idx = NULL;
  double *X_perm_back_val = NULL;
  New_1D_Array<INTEGER, INTEGER>(&X_perm_back_start, N+1);
  New_1D_Array<INTEGER, INTEGER>(&X_perm_back_idx, nnz);
  New_1D_Array<double, INTEGER>(&X_perm_back_val, nnz);
  INTEGER *X_ptr_start = X.GetPointerStart();
  INTEGER *X_ptr_idx = X.GetPointerIdx();
  double *X_ptr_val = X.GetPointerX();

  X_ptr_start[0] = 0;
  INTEGER current = 0;
  for (INTEGER i = 0; i < N; i++) {
    INTEGER j = perm[i];
    INTEGER start_j = X_ptr_start[j];
    for (INTEGER k = 0; k < X_ptr_start[j+1] - X_ptr_start[j]; k++) {
      X_perm_back_idx[current] = X_ptr_idx[start_j+k];
      X_perm_back_val[current] = X_ptr_val[start_j+k];
      current++;
    }
    X_perm_back_start[i+1] = X_perm_back_start[i] + X_ptr_start[j+1] - X_ptr_start[j];
  }

  double derr = Diff1(X_perm_back_start, start, N+1);
  printf("Test %s() : Part 1, Error in start %g\n", FuncName, derr);
  derr = Diff1(X_perm_back_idx, idx, nnz);
  printf("Test %s() : Part 1, Error in idx %g\n", FuncName, derr);
  derr = Diff1(X_perm_back_val, val, nnz);
  printf("Test %s() : Part 1, Error in val %g\n", FuncName, derr);

  // Check that outside [istart, istart+n-1], perm[i] = i
  derr = 0.0;
  for (INTEGER i = 0; i < istart; i++) {
    if (perm[i] != i) {
      derr += 1.0;
    }
  }
  for (INTEGER i = istart+n; i < N; i++) {
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
  for (INTEGER i = istart; i < istart+m1; i++) {
    if (y_ptr[i] > 0) {
      derr += 1.0;
    }
  }
  for (INTEGER i = istart+m1; i < istart+n; i++) {
    if (y_ptr[i] < 0) {
      derr += 1.0;
    }
  }
  printf("Test %s() : Part 3, Error in sign %g\n", FuncName, derr);

  // Clean up
  Delete_1D_Array<INTEGER>(&X_perm_back_start);
  Delete_1D_Array<INTEGER>(&X_perm_back_idx);
  Delete_1D_Array<double>(&X_perm_back_val);

}
