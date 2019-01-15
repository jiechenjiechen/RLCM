// Unit test of the utility routines implemented in Common.hpp.

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
  printf("|  Test_Common.cpp                                              |\n");
  printf("+---------------------------------------------------------------+\n");

  // Test Square() ------------------------------------------------------------
  //
  // (2.1)^2 = 4.41
  double dx = 2.1;
  double dy = Square(dx);
  double dy_truth = 4.41;
  double derr = fabs(dy - dy_truth);
  printf("Test Square() : Discrepancy in return value %g\n", derr);

  // Test Parity() ------------------------------------------------------------
  //
  // a = [2 5 2 4 4 3 1 6 8 7], n = 10, parity = 1
  INTEGER n = 10;
  INTEGER a[10] = {2, 5, 2, 4, 4, 3, 1, 6, 8, 7};
  INTEGER p = Parity(a, n);
  int p_truth = 1;
  int ierr = abs((int)p - p_truth);
  printf("Test Parity() : Discrepancy in return value %ld\n", (long)ierr);

  // Test UniformRandom01() ---------------------------------------------------
  //
  n = 10000;
  double b[10000] = {0};
  UniformRandom01(b, n);
  const char *out_file1 = "Test_Common_UniformRandom01.txt";
  FILE *fid = fopen(out_file1, "w");
  for (INTEGER i = 0; i < n; i++) {
    fprintf(fid, "%g\n", b[i]);
  }
  fclose(fid);
  printf("Test UniformRandom01() : Run Test_Common_RandomNumbers.m to verify correctness.\n");

  // Test StandardNormal() ----------------------------------------------------
  //
  StandardNormal(b, n);
  const char *out_file2 = "Test_Common_StandardNormal.txt";
  fid = fopen(out_file2, "w");
  for (INTEGER i = 0; i < n; i++) {
    fprintf(fid, "%g\n", b[i]);
  }
  fclose(fid);
  printf("Test StandardNormal() : Run Test_Common_RandomNumbers.m to verify correctness.\n");

  // Test StudentT1() ---------------------------------------------------------
  //
  StudentT1(b, n);
  const char *out_file3 = "Test_Common_StudentT1.txt";
  fid = fopen(out_file3, "w");
  for (INTEGER i = 0; i < n; i++) {
    fprintf(fid, "%g\n", b[i]);
  }
  fclose(fid);
  printf("Test StudentT1() : Run Test_Common_RandomNumbers.m to verify correctness.\n");

  // Test MultivariateStudentT1() ---------------------------------------------
  //
  const char *out_file4 = "Test_Common_MultivariateStudentT1.txt";
  fid = fopen(out_file4, "w");
  for (INTEGER i = 0; i < n; i++) {
    MultivariateStudentT1(b, 2);
    fprintf(fid, "%g %g\n", b[0], b[1]);
  }
  fclose(fid);
  printf("Test MultivariateStudentT1() : Run Test_Common_RandomNumbers.m to verify correctness.\n");

  // Test RandomSech() --------------------------------------------------------
  //
  RandomSech(b, n);
  const char *out_file5 = "Test_Common_RandomSech.txt";
  fid = fopen(out_file5, "w");
  for (INTEGER i = 0; i < n; i++) {
    fprintf(fid, "%g\n", b[i]);
  }
  fclose(fid);
  printf("Test RandomSech() : Run Test_Common_RandomNumbers.m to verify correctness.\n");

  // Test RandPerm() ----------------------------------------------------------
  //
  // n = 10, k = 5. There must show 5 non-repeating integers betwenn 0 and 9
  n = 10;
  INTEGER k = 5;
  RandPerm(n, k, a); // Array a declared above has sufficient length
  printf("Test RandPerm() :");
  for (INTEGER i = 0; i < k; i++) {
    printf(" %ld", (long)a[i]);
  }
  printf(". (There must show 5 non-repeating integers betwenn 0 and 9.)\n");

  // Test CompareNaturalOrderLess() -------------------------------------------
  //
  // x = 1.0, y = 2.0, return -1
  // x = 3.0, y = 2.5, return  1
  // x = 1.5, y = 1.5, return  0
  int iret, iret_truth;
  dx = 1.0; dy = 2.0; iret_truth = -1;
  iret = CompareNaturalOrderLess(&dx, &dy);
  ierr = abs(iret - iret_truth);
  printf("Test CompareNaturalOrderLess() : Discrepancy in return value %ld\n", (long)ierr);

  dx = 3.0; dy = 2.5; iret_truth =  1;
  iret = CompareNaturalOrderLess(&dx, &dy);
  ierr = abs(iret - iret_truth);
  printf("Test CompareNaturalOrderLess() : Discrepancy in return value %ld\n", (long)ierr);

  dx = 1.5; dy = 1.5; iret_truth =  0;
  iret = CompareNaturalOrderLess(&dx, &dy);
  ierr = abs(iret - iret_truth);
  printf("Test CompareNaturalOrderLess() : Discrepancy in return value %ld\n", (long)ierr);

  // Test CompareNaturalOrderGreater() ----------------------------------------
  //
  // x = 1.0, y = 2.0, return  1
  // x = 3.0, y = 2.5, return -1
  // x = 1.5, y = 1.5, return  0
  dx = 1.0; dy = 2.0; iret_truth =  1;
  iret = CompareNaturalOrderGreater(&dx, &dy);
  ierr = abs(iret - iret_truth);
  printf("Test CompareNaturalOrderGreater() : Discrepancy in return value %ld\n", (long)ierr);

  dx = 3.0; dy = 2.5; iret_truth = -1;
  iret = CompareNaturalOrderGreater(&dx, &dy);
  ierr = abs(iret - iret_truth);
  printf("Test CompareNaturalOrderGreater() : Discrepancy in return value %ld\n", (long)ierr);

  dx = 1.5; dy = 1.5; iret_truth =  0;
  iret = CompareNaturalOrderGreater(&dx, &dy);
  ierr = abs(iret - iret_truth);
  printf("Test CompareNaturalOrderGreater() : Discrepancy in return value %ld\n", (long)ierr);

  // Test CompareIntegerNaturalOrderLess() ------------------------------------
  //
  // x = 1, y = 2, return -1
  // x = 3, y = 2, return  1
  // x = 1, y = 1, return  0
  INTEGER ix, iy;
  ix = 1; iy = 2; iret_truth = -1;
  iret = CompareIntegerNaturalOrderLess(&ix, &iy);
  ierr = abs(iret - iret_truth);
  printf("Test CompareIntegerNaturalOrderLess() : Discrepancy in return value %ld\n", (long)ierr);

  ix = 3; iy = 2; iret_truth =  1;
  iret = CompareIntegerNaturalOrderLess(&ix, &iy);
  ierr = abs(iret - iret_truth);
  printf("Test CompareIntegerNaturalOrderLess() : Discrepancy in return value %ld\n", (long)ierr);

  ix = 1; iy = 1; iret_truth =  0;
  iret = CompareIntegerNaturalOrderLess(&ix, &iy);
  ierr = abs(iret - iret_truth);
  printf("Test CompareIntegerNaturalOrderLess() : Discrepancy in return value %ld\n", (long)ierr);

  // Test CompareIntegerNaturalOrderGreater() ---------------------------------
  //
  // x = 1, y = 2, return  1
  // x = 3, y = 2, return -1
  // x = 1, y = 1, return  0
  ix = 1; iy = 2; iret_truth =  1;
  iret = CompareIntegerNaturalOrderGreater(&ix, &iy);
  ierr = abs(iret - iret_truth);
  printf("Test CompareIntegerNaturalOrderGreater() : Discrepancy in return value %ld\n", (long)ierr);

  ix = 3; iy = 2; iret_truth = -1;
  iret = CompareIntegerNaturalOrderGreater(&ix, &iy);
  ierr = abs(iret - iret_truth);
  printf("Test CompareIntegerNaturalOrderGreater() : Discrepancy in return value %ld\n", (long)ierr);

  ix = 1; iy = 1; iret_truth =  0;
  iret = CompareIntegerNaturalOrderGreater(&ix, &iy);
  ierr = abs(iret - iret_truth);
  printf("Test CompareIntegerNaturalOrderGreater() : Discrepancy in return value %ld\n", (long)ierr);

  // Test CompareByMagnitudeLess() --------------------------------------------
  //
  // x =  1.0, y = -2.0, return -1
  // x = -3.0, y =  2.0, return  1
  // x = -1.0, y =  1.0, return  0
  dx =  1.0; dy = -2.0; iret_truth = -1;
  iret = CompareByMagnitudeLess(&dx, &dy);
  ierr = abs(iret - iret_truth);
  printf("Test CompareByMagnitudeLess() : Discrepancy in return value %ld\n", (long)ierr);

  dx = -3.0; dy =  2.0; iret_truth =  1;
  iret = CompareByMagnitudeLess(&dx, &dy);
  ierr = abs(iret - iret_truth);
  printf("Test CompareByMagnitudeLess() : Discrepancy in return value %ld\n", (long)ierr);

  dx = -1.0; dy =  1.0; iret_truth =  0;
  iret = CompareByMagnitudeLess(&dx, &dy);
  ierr = abs(iret - iret_truth);
  printf("Test CompareByMagnitudeLess() : Discrepancy in return value %ld\n", (long)ierr);

  // Test CompareByMagnitudeGreater() -----------------------------------------
  //
  // x =  1.0, y = -2.0, return  1
  // x = -3.0, y =  2.0, return -1
  // x = -1.0, y =  1.0, return  0
  dx =  1.0; dy = -2.0; iret_truth =  1;
  iret = CompareByMagnitudeGreater(&dx, &dy);
  ierr = abs(iret - iret_truth);
  printf("Test CompareByMagnitudeGreater() : Discrepancy in return value %ld\n", (long)ierr);

  dx = -3.0; dy =  2.0; iret_truth = -1;
  iret = CompareByMagnitudeGreater(&dx, &dy);
  ierr = abs(iret - iret_truth);
  printf("Test CompareByMagnitudeGreater() : Discrepancy in return value %ld\n", (long)ierr);

  dx = -1.0; dy =  1.0; iret_truth =  0;
  iret = CompareByMagnitudeGreater(&dx, &dy);
  ierr = abs(iret - iret_truth);
  printf("Test CompareByMagnitudeGreater() : Discrepancy in return value %ld\n", (long)ierr);

  // Test SelectLeftHalfPlane() -----------------------------------------------
  //
  // WR =  1.0, WI = -1.0, return 0
  // WR = -1.0, WI =  2.0, return 1
  double WR, WI;
  WR =  1.0; WI = -1.0; iret_truth = 0;
  iret = SelectLeftHalfPlane(&WR, &WI);
  ierr = abs(iret - iret_truth);
  printf("Test SelectLeftHalfPlane() : Discrepancy in return value %ld\n", (long)ierr);

  WR = -1.0; WI =  2.0; iret_truth = 1;
  iret = SelectLeftHalfPlane(&WR, &WI);
  ierr = abs(iret - iret_truth);
  printf("Test SelectLeftHalfPlane() : Discrepancy in return value %ld\n", (long)ierr);

  // Test SelectRightHalfPlane() ----------------------------------------------
  //
  // WR = -1.0, WI =  2.0, return 0
  // WR =  1.0, WI = -1.0, return 1
  WR = -1.0; WI =  2.0; iret_truth = 0;
  iret = SelectRightHalfPlane(&WR, &WI);
  ierr = abs(iret - iret_truth);
  printf("Test SelectRightHalfPlane() : Discrepancy in return value %ld\n", (long)ierr);

  WR =  1.0; WI = -1.0; iret_truth = 1;
  iret = SelectRightHalfPlane(&WR, &WI);
  ierr = abs(iret - iret_truth);
  printf("Test SelectRightHalfPlane() : Discrepancy in return value %ld\n", (long)ierr);

  // Test CompareElm1Idx() ----------------------------------------------------
  //
  // x.idx = 6, x.val = 1.0, y.idx = 7, y.val = 3.0, return -1
  // x.idx = 7, x.val = 2.0, y.idx = 4, y.val = 2.0, return  1
  // x.idx = 5, x.val = 3.0, y.idx = 5, y.val = 1.0, return  0
  Elm1 e1x, e1y;
  e1x = { 6, 1.0 }; e1y = { 7, 3.0 }; iret_truth = -1;
  iret = CompareElm1Idx(&e1x, &e1y);
  ierr = abs(iret - iret_truth);
  printf("Test CompareElm1Idx() : Discrepancy in return value %ld\n", (long)ierr);

  e1x = { 7, 2.0 }; e1y = { 4, 2.0 }; iret_truth =  1;
  iret = CompareElm1Idx(&e1x, &e1y);
  ierr = abs(iret - iret_truth);
  printf("Test CompareElm1Idx() : Discrepancy in return value %ld\n", (long)ierr);

  e1x = { 5, 3.0 }; e1y = { 5, 1.0 }; iret_truth =  0;
  iret = CompareElm1Idx(&e1x, &e1y);
  ierr = abs(iret - iret_truth);
  printf("Test CompareElm1Idx() : Discrepancy in return value %ld\n", (long)ierr);

  // Test CompareElm2Idx() ----------------------------------------------------
  //
  // x.idx = 6, x.above = 1, x.below = 3, y.idx = 7, y.above = 3, y.below = 1, 
  //                                                                  return -1
  // x.idx = 7, x.above = 2, x.below = 2, y.idx = 4, y.above = 2, x.below = 2, 
  //                                                                  return  1
  // x.idx = 5, x.above = 3, x.below = 1, y.idx = 5, y.above = 1, x.below = 3, 
  //                                                                  return  0
  Elm2 e2x, e2y;
  e2x = { 6, 1, 3 }; e2y = { 7, 3, 1 }; iret_truth = -1;
  iret = CompareElm2Idx(&e2x, &e2y);
  ierr = abs(iret - iret_truth);
  printf("Test CompareElm2Idx() : Discrepancy in return value %ld\n", (long)ierr);

  e2x = { 7, 2, 2 }; e2y = { 4, 2, 2 }; iret_truth =  1;
  iret = CompareElm2Idx(&e2x, &e2y);
  ierr = abs(iret - iret_truth);
  printf("Test CompareElm2Idx() : Discrepancy in return value %ld\n", (long)ierr);

  e2x = { 5, 3, 1 }; e2y = { 5, 1, 3 }; iret_truth =  0;
  iret = CompareElm2Idx(&e2x, &e2y);
  ierr = abs(iret - iret_truth);
  printf("Test CompareElm2Idx() : Discrepancy in return value %ld\n", (long)ierr);

  // Test CompareElm2Above() --------------------------------------------------
  //
  // x.idx = 6, x.above = 1, x.below = 3, y.idx = 7, y.above = 3, y.below = 1, 
  //                                                                  return -1
  // x.idx = 7, x.above = 2, x.below = 2, y.idx = 4, y.above = 2, x.below = 2, 
  //                                                                  return  0
  // x.idx = 5, x.above = 3, x.below = 1, y.idx = 5, y.above = 1, x.below = 3, 
  //                                                                  return  1
  e2x = { 6, 1, 3 }; e2y = { 7, 3, 1 }; iret_truth = -1;
  iret = CompareElm2Above(&e2x, &e2y);
  ierr = abs(iret - iret_truth);
  printf("Test CompareElm2Above() : Discrepancy in return value %ld\n", (long)ierr);

  e2x = { 7, 2, 2 }; e2y = { 4, 2, 2 }; iret_truth =  0;
  iret = CompareElm2Above(&e2x, &e2y);
  ierr = abs(iret - iret_truth);
  printf("Test CompareElm2Above() : Discrepancy in return value %ld\n", (long)ierr);

  e2x = { 5, 3, 1 }; e2y = { 5, 1, 3 }; iret_truth =  1;
  iret = CompareElm2Above(&e2x, &e2y);
  ierr = abs(iret - iret_truth);
  printf("Test CompareElm2Above() : Discrepancy in return value %ld\n", (long)ierr);

  // Test CompareElm2Below() --------------------------------------------------
  //
  // x.idx = 6, x.above = 1, x.below = 3, y.idx = 7, y.above = 3, y.below = 1, 
  //                                                                  return  1
  // x.idx = 7, x.above = 2, x.below = 2, y.idx = 4, y.above = 2, x.below = 2, 
  //                                                                  return  0
  // x.idx = 5, x.above = 3, x.below = 1, y.idx = 5, y.above = 1, x.below = 3, 
  //                                                                  return -1
  e2x = { 6, 1, 3 }; e2y = { 7, 3, 1 }; iret_truth =  1;
  iret = CompareElm2Below(&e2x, &e2y);
  ierr = abs(iret - iret_truth);
  printf("Test CompareElm2Below() : Discrepancy in return value %ld\n", (long)ierr);

  e2x = { 7, 2, 2 }; e2y = { 4, 2, 2 }; iret_truth =  0;
  iret = CompareElm2Below(&e2x, &e2y);
  ierr = abs(iret - iret_truth);
  printf("Test CompareElm2Below() : Discrepancy in return value %ld\n", (long)ierr);

  e2x = { 5, 3, 1 }; e2y = { 5, 1, 3 }; iret_truth = -1;
  iret = CompareElm2Below(&e2x, &e2y);
  ierr = abs(iret - iret_truth);
  printf("Test CompareElm2Below() : Discrepancy in return value %ld\n", (long)ierr);

  // Test CompareElm3Dim() ----------------------------------------------------
  //
  // x.idx = 6, x.dim = 1, x.len = 0.1, x.inc = 0.1, x.low = 0.1
  // y.idx = 7, y.dim = 3, y.len = 0.1, y.inc = 0.1, y.low = 0.1, return -1
  // x.idx = 7, x.dim = 2, x.len = 0.1, x.inc = 0.1, x.low = 0.1
  // y.idx = 4, y.dim = 2, y.len = 0.1, y.inc = 0.1, y.low = 0.1, return  0
  // x.idx = 5, x.dim = 3, x.len = 0.1, x.inc = 0.1, x.low = 0.1
  // y.idx = 5, y.dim = 1, y.len = 0.1, y.inc = 0.1, y.low = 0.1, return  1
  Elm3 e3x, e3y;
  e3x = { 6, 1, 0.1, 0.1, 0.1 }; e3y = { 7, 3, 0.1, 0.1, 0.1 }; iret_truth = -1;
  iret = CompareElm3Dim(&e3x, &e3y);
  ierr = abs(iret - iret_truth);
  printf("Test CompareElm3Dim() : Discrepancy in return value %ld\n", (long)ierr);

  e3x = { 7, 2, 0.1, 0.1, 0.1 }; e3y = { 4, 2, 0.1, 0.1, 0.1 }; iret_truth =  0;
  iret = CompareElm3Dim(&e3x, &e3y);
  ierr = abs(iret - iret_truth);
  printf("Test CompareElm3Dim() : Discrepancy in return value %ld\n", (long)ierr);

  e3x = { 5, 3, 0.1, 0.1, 0.1 }; e3y = { 5, 1, 0.1, 0.1, 0.1 }; iret_truth =  1;
  iret = CompareElm3Dim(&e3x, &e3y);
  ierr = abs(iret - iret_truth);
  printf("Test CompareElm3Dim() : Discrepancy in return value %ld\n", (long)ierr);

  // Test Swap() --------------------------------------------------------------
  //
  // x = 1.0, y = 0.2
  dx = 1.0; dy = 0.2;
  double dx_after = dy, dy_after = dx;
  Swap<double>(dx, dy);
  double derr_dx = fabs(dx - dx_after);
  double derr_dy = fabs(dy - dy_after);
  printf("Test Swap() : Discrepancy in results %g, %g\n", derr_dx, derr_dy);

  // Test Max() ---------------------------------------------------------------
  //
  // Reuse x = 1.0, y = 0.2
  double max = Max<double>(dx, dy);
  double max_truth = 1.0;
  derr = fabs(max - max_truth);
  printf("Test Max() : Discrepancy in results %g\n", derr);

  // Test Min() ---------------------------------------------------------------
  //
  // Reuse x = 1.0, y = 0.2
  double min = Min<double>(dx, dy);
  double min_truth = 0.2;
  derr = fabs(min - min_truth);
  printf("Test Min() : Discrepancy in results %g\n", derr);

  // Test zbesk_() ------------------------------------------------------------
  //
  // nu = 3.1, x = 1.2
#ifdef USE_MATERN
#define LEN 56
  double fnu = 3.1;
  double zr = 1.2;
  double zi = 0;
  int kode = 1;
  int nn = 1;
  double cyr[1];
  double cyi[1];
  int nz;
  zbesk_(&zr, &zi, &fnu, &kode, &nn, cyr, cyi, &nz, &ierr);
  char output[LEN+1] = "cyr[0] = 4.55423672200606, cyi[0] = 0, nz = 0, ierr = 0\n";
  char output2[LEN+1];
  snprintf(output2, LEN, "cyr[0] = %.14f, cyi[0] = %g, nz = %d, ierr = %d\n",
           cyr[0], cyi[0], nz, ierr);
  output2[LEN-1] = '\n';
  output2[LEN] = '\0';
  if (strncmp(output, output2, LEN) != 0) {
    printf("Test zbesk_() : Failed.\n");
    printf("Correct output: %s", output);
    printf("Computed output: %s", output2);
  }
  else {
    printf("Test zbesk_() : Success.\n");
  }
#endif

  return 0;

}
