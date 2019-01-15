// Unit test of the sparse linear algebra routines.

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
  printf("|  Test_spblas.cpp                                              |\n");
  printf("+---------------------------------------------------------------+\n");

  // Test sp_n_ds2d() ---------------------------------------------------------
  //
  // sx = [ 0.0; 0.4; 0.3; 0.0 ]
  // dy_truth = [ 0.0; 0.4; 0.3; 0.0 ]
  INTEGER nnzx = 2;
  INTEGER idxx[4] = { 1, 2 };
  double sx[4] = { 0.4, 0.3 };
  INTEGER n = 4;
  double dy[4] = { 0.0 };
  double dy_truth[] = { 0.0, 0.4, 0.3, 0.0 };
  sp_n_ds2d(n, nnzx, idxx, sx, dy);
  double derr = Diff1(dy, dy_truth, n);
  printf("Test sp_n_ds2d() : Discrepancy in return vector %g\n", derr);

  // Test sp_d_ddot() ---------------------------------------------------------
  //
  // sx: reuse above
  // dy = [ 0.1; 0.0; -0.2; 1.0 ]
  // z = -0.06
  dy[0] = 0.1; dy[1] = 0.0; dy[2] = -0.2; dy[3] = 1.0;
  double z = sp_d_ddot(n, nnzx, idxx, sx, dy);
  double z_truth = -0.06;
  derr = fabs(z - z_truth);
  printf("Test sp_d_ddot() : Discrepancy in return value %g\n", derr);

  // Test sp_s_ddot() ---------------------------------------------------------
  //
  // sx: reuse above
  // sy = [ 0.1; 0.0; -0.2; 1.0 ]
  // z = -0.06
  INTEGER nnzy = 3;
  INTEGER idxy[4] = { 0, 2, 3 };
  double sy[4] = { 0.1, -0.2, 1.0 };
  z = sp_s_ddot(n, nnzx, idxx, sx, nnzy, idxy, sy);
  z_truth = -0.06;
  derr = fabs(z - z_truth);
  printf("Test sp_s_ddot() : Discrepancy in return value %g\n", derr);

  // Test sp_d_ddist() --------------------------------------------------------
  //
  // sx: reuse above
  // dy: reuse above
  // func = fabs
  // z = 2.0
  z = sp_d_ddist(n, nnzx, idxx, sx, dy, fabs);
  z_truth = 2.0;
  derr = fabs(z - z_truth);
  printf("Test sp_d_ddist() : Discrepancy in return value %g\n", derr);

  // Test sp_s_ddist() --------------------------------------------------------
  //
  // sx: reuse above
  // sy: reuse above
  // func = Square
  // z = 1.42
  z = sp_s_ddist(n, nnzx, idxx, sx, nnzy, idxy, sy, Square);
  z_truth = 1.42;
  derr = fabs(z - z_truth);
  printf("Test sp_s_ddist() : Discrepancy in return value %g\n", derr);

  // Test sp_d_ddists() --------------------------------------------------------
  //
  // sx: reuse above
  // dy: reuse above
  // func = fabs
  // sigma = [ 2.0; 1.0; 0.5; 2.0 ]
  // z = 1.95
  double sigma[] = { 2.0, 1.0, 0.5, 2.0 };
  z = sp_d_ddists(n, nnzx, idxx, sx, dy, sigma, fabs);
  z_truth = 1.95;
  derr = fabs(z - z_truth);
  printf("Test sp_d_ddists() : Discrepancy in return value %g\n", derr);

  // Test sp_s_ddists() --------------------------------------------------------
  //
  // sx: reuse above
  // sy: reuse above
  // sigma: reuse above
  // func = Square
  // z = 1.4125
  z = sp_s_ddists(n, nnzx, idxx, sx, nnzy, idxy, sy, sigma, Square);
  z_truth = 1.4125;
  derr = fabs(z - z_truth);
  printf("Test sp_s_ddists() : Discrepancy in return value %g\n", derr);

  // Test sp_d_dchi2() --------------------------------------------------------
  //
  // sx = [ 0.0; 0.4; 0.3; 0.0 ] (reuse above)
  // dy = [ 0.1; 0.0; 0.2; 0.0 ]
  // z = 0.24
  dy[0] = 0.1; dy[1] = 0.0; dy[2] = 0.2; dy[3] = 0.0;
  z = sp_d_dchi2(n, nnzx, idxx, sx, dy);
  z_truth = 0.24;
  derr = fabs(z - z_truth);
  printf("Test sp_d_dchi2() : Discrepancy in return value %g\n", derr);

  // Test sp_s_dchi2() --------------------------------------------------------
  //
  // sx = [ 0.0; 0.4; 0.3; 0.0 ] (reuse above)
  // sy = [ 0.1; 0.0; 0.2; 0.0 ]
  // z = 0.24
  nnzy = 2;
  idxy[0] = 0; idxy[1] = 2; 
  sy[0] = 0.1; sy[1] = 0.2;
  z = sp_s_dchi2(n, nnzx, idxx, sx, nnzy, idxy, sy);
  z_truth = 0.24;
  derr = fabs(z - z_truth);
  printf("Test sp_s_dchi2() : Discrepancy in return value %g\n", derr);

  // Test sp_d_dgemv() --------------------------------------------------------
  //
  // [  0.29 ]   [  1.0  0.0  0.3  0.4 ]   [  0.3 ]
  // [  0.08 ] = [  0.0  0.0  1.0  0.2 ] * [  0.0 ]
  // [ -0.03 ]   [ -0.1  1.0  0.0  0.0 ]   [  0.1 ]
  //                                       [ -0.1 ]
  // TRANS = 'N'
  INTEGER mA1 = 3;
  INTEGER nA1 = 4;
  INTEGER nnzA1 = 7;
  INTEGER startA1[] = { 0, 3, 5, 7 };
  INTEGER idxA1[] = { 0, 2, 3, 2, 3, 0, 1 };
  double SA1[] = { 1.0, 0.3, 0.4, 1.0, 0.2, -0.1, 1.0 };
  double dx[] = { 0.3, 0.0, 0.1, -0.1 };
  dy_truth[0] = 0.29; dy_truth[1] = 0.08; dy_truth[2] = -0.03;
  sp_d_dgemv('N', mA1, nA1, nnzA1, startA1, idxA1, SA1, dx, dy);
  derr = Diff1(dy, dy_truth, 3);
  printf("Test sp_d_dgemv() : Discrepancy in return vector %g\n", derr);

  // TRANS = 'T'
  INTEGER mA2 = 4;
  INTEGER nA2 = 3;
  INTEGER nnzA2 = 7;
  INTEGER startA2[] = { 0, 2, 3, 5, 7 };
  INTEGER idxA2[] = { 0, 2, 2, 0, 1, 0, 1 };
  double SA2[] = { 1.0, -0.1, 1.0, 0.3, 1.0, 0.4, 0.2 };
  sp_d_dgemv('T', mA2, nA2, nnzA2, startA2, idxA2, SA2, dx, dy);
  derr = Diff1(dy, dy_truth, 3);
  printf("Test sp_d_dgemv() : Discrepancy in return vector %g\n", derr);

  // Test sp_s_dgemv() --------------------------------------------------------
  //
  // same as above
  // TRANS = 'N'
  nnzx = 3;
  idxx[0] = 0; idxx[1] = 2; idxx[2] = 3;
  sx[0] = 0.3; sx[1] = 0.1; sx[2] = -0.1;
  // Arrays idxx and dx declared above has sufficient length
  sp_s_dgemv('N', mA1, nA1, nnzA1, startA1, idxA1, SA1, nnzx, idxx, sx, dy);
  derr = Diff1(dy, dy_truth, 3);
  printf("Test sp_s_dgemv() : Discrepancy in return vector %g\n", derr);

  // TRANS = 'T'
  sp_s_dgemv('T', mA2, nA2, nnzA2, startA2, idxA2, SA2, nnzx, idxx, sx, dy);
  derr = Diff1(dy, dy_truth, 3);
  printf("Test sp_s_dgemv() : Discrepancy in return vector %g\n", derr);

  // Test sp_d_dgemm() --------------------------------------------------------
  //
  // [  0.03   0.60 ]   [  1.0  0.0  0.3  0.4 ]   [  0.0   1.0 ]
  // [  0.10  -0.20 ] = [  0.0  0.0  1.0  0.2 ] * [ -0.4   0.3 ]
  // [ -0.40   0.20 ]   [ -0.1  1.0  0.0  0.0 ]   [  0.1   0.0 ]
  //                                              [  0.0  -1.0 ]
  // TRANSA = 'N', TRANSB = 'N'
  INTEGER m = 3;
  n = 2;
  INTEGER k = 4;
  double DB1[] = { 0.0, -0.4, 0.1, 0.0, 1.0, 0.3, 0.0, -1.0 };
  double DC[6];
  double DC_truth[] = { 0.03, 0.10, -0.40, 0.60, -0.20, 0.20 };
  sp_d_dgemm('N', 'N', m, n, k, nnzA1, startA1, idxA1, SA1, DB1, DC);
  derr = Diff1(DC, DC_truth, 6);
  printf("Test sp_d_dgemm() : Discrepancy in return matrix %g\n", derr);

  // TRANSA = 'N', TRANSB = 'T'
  double DB2[] = { 0.0, 1.0, -0.4, 0.3, 0.1, 0.0, 0.0, -1.0 };
  sp_d_dgemm('N', 'T', m, n, k, nnzA1, startA1, idxA1, SA1, DB2, DC);
  derr = Diff1(DC, DC_truth, 6);
  printf("Test sp_d_dgemm() : Discrepancy in return matrix %g\n", derr);

  // TRANSA = 'T', TRANSB = 'N'
  sp_d_dgemm('T', 'N', m, n, k, nnzA2, startA2, idxA2, SA2, DB1, DC);
  derr = Diff1(DC, DC_truth, 6);
  printf("Test sp_d_dgemm() : Discrepancy in return matrix %g\n", derr);

  // TRANSA = 'T', TRANSB = 'T'
  sp_d_dgemm('T', 'T', m, n, k, nnzA2, startA2, idxA2, SA2, DB2, DC);
  derr = Diff1(DC, DC_truth, 6);
  printf("Test sp_d_dgemm() : Discrepancy in return matrix %g\n", derr);

  // Test sp_s_dgemm() --------------------------------------------------------
  //
  // same as above
  // TRANSA = 'N', TRANSB = 'N'
  INTEGER nnzB1 = 5;
  INTEGER startB1[] = { 0, 1, 3, 4, 5 };
  INTEGER idxB1[] = { 1, 0, 1, 0, 1 };
  double SB1[] = { 1.0, -0.4, 0.3, 0.1, -1.0 };
  sp_s_dgemm('N', 'N', m, n, k, nnzA1, startA1, idxA1, SA1, nnzB1, startB1, idxB1, SB1, DC);
  derr = Diff1(DC, DC_truth, 6);
  printf("Test sp_s_dgemm() : Discrepancy in return matrix %g\n", derr);

  // TRANSA = 'N', TRANSB = 'T'
  INTEGER nnzB2 = 5;
  INTEGER startB2[] = { 0, 2, 5 };
  INTEGER idxB2[] = { 1, 2, 0, 1, 3 };
  double SB2[] = { -0.4, 0.1, 1.0, 0.3, -1.0 };
  sp_s_dgemm('N', 'T', m, n, k, nnzA1, startA1, idxA1, SA1, nnzB2, startB2, idxB2, SB2, DC);
  derr = Diff1(DC, DC_truth, 6);
  printf("Test sp_s_dgemm() : Discrepancy in return matrix %g\n", derr);

  // TRANSA = 'T', TRANSB = 'N'
  sp_s_dgemm('T', 'N', m, n, k, nnzA2, startA2, idxA2, SA2, nnzB1, startB1, idxB1, SB1, DC);
  derr = Diff1(DC, DC_truth, 6);
  printf("Test sp_s_dgemm() : Discrepancy in return matrix %g\n", derr);

  // TRANSA = 'T', TRANSB = 'T'
  sp_s_dgemm('T', 'T', m, n, k, nnzA2, startA2, idxA2, SA2, nnzB2, startB2, idxB2, SB2, DC);
  derr = Diff1(DC, DC_truth, 6);
  printf("Test sp_s_dgemm() : Discrepancy in return matrix %g\n", derr);

  return 0;

}
