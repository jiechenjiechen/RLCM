// Unit test of the class DPoint.

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
  printf("|  Test_DPoint.cpp                                              |\n");
  printf("+---------------------------------------------------------------+\n");

  // Test Normalize() ---------------------------------------------------------
  //
  // x = [ 0.0; 0.4; 0.3; 0.0 ]
  // normalized = [ 0.0; 0.8; 0.6; 0.0 ]
  INTEGER d = 4;
  double mx[] = { 0.0, 0.4, 0.3, 0.0 };
  DPoint x;
  x.SetPoint(mx, d);
  x.Normalize();
  double normalize_truth[] = { 0.0, 0.8, 0.6, 0.0 };
  double derr = Diff1(x.GetPointer(), normalize_truth, d);
  printf("Test Normalize() : Discrepancy in return vector %g\n", derr);

  // Test InProd() ------------------------------------------------------------
  //
  // x = [ 0.0; 0.4; 0.3; 0.0 ]
  // y = [ 0.1; 0.0; -0.2; 1.0 ]
  // inprod = -0.06
  x.SetPoint(mx, d);
  double my[] = { 0.1, 0.0, -0.2, 1.0 };
  DPoint y;
  y.SetPoint(my, d);
  double inprod = x.InProd(y);
  double inprod_truth = -0.06;
  derr = fabs(inprod-inprod_truth);
  printf("Test InProd() : Discrepancy in return value %g\n", derr);

  // Test Dist1() -------------------------------------------------------------
  //
  // Use the above x and y
  // dist1 = 2.0
  double dist1 = x.Dist1(y);
  double dist1_truth = 2.0;
  derr = fabs(dist1-dist1_truth);
  printf("Test Dist1() : Discrepancy in return value %g\n", derr);

  // Test Dist1() -------------------------------------------------------------
  //
  // Use the above x and y
  // sigma = [ 2.0; 1.0; 0.5; 2.0 ]
  // dist1 = 1.95
  double sigma[] = { 2.0, 1.0, 0.5, 2.0 };
  dist1 = x.Dist1(y, sigma);
  dist1_truth = 1.95;
  derr = fabs(dist1-dist1_truth);
  printf("Test Dist1() : Discrepancy in return value %g\n", derr);

  // Test Dist2() -------------------------------------------------------------
  //
  // Use the above x and y
  // dist1 = 1.42
  double dist2 = x.Dist2(y);
  double dist2_truth = 1.42;
  derr = fabs(dist2-dist2_truth);
  printf("Test Dist2() : Discrepancy in return value %g\n", derr);

  // Test Dist2() -------------------------------------------------------------
  //
  // Use the above x, y and sigma
  // dist1 = 1.4125
  dist2 = x.Dist2(y, sigma);
  dist2_truth = 1.4125;
  derr = fabs(dist2-dist2_truth);
  printf("Test Dist2() : Discrepancy in return value %g\n", derr);

  // Test Subtract() ----------------------------------------------------------
  //
  // Use the above x and y
  // z = [ -0.1; 0.4; 0.5; -1.0 ]
  x.Subtract(y);
  double subtract_truth[] = { -0.1, 0.4, 0.5, -1.0 };
  derr = Diff1(x.GetPointer(), subtract_truth, d);
  printf("Test Subtract() : Discrepancy in return vector %g\n", derr);

  x.SetPoint(mx, d);
  DPoint z;
  x.Subtract(y, z);
  derr = Diff1(z.GetPointer(), subtract_truth, d);
  printf("Test Subtract() : Discrepancy in return vector %g\n", derr);

  // Test AverageWith() -------------------------------------------------------
  //
  // Use the above x and y
  // z = [ 0.05; 0.2; 0.05; 0.5 ]
  x.SetPoint(mx, d);
  x.AverageWith(y);
  double ave_truth[] = { 0.05, 0.2, 0.05, 0.5 };
  derr = Diff1(x.GetPointer(), ave_truth, d);
  printf("Test AverageWith() : Discrepancy in return vector %g\n", derr);

  x.SetPoint(mx, d);
  x.AverageWith(y, z);
  derr = Diff1(z.GetPointer(), ave_truth, d);
  printf("Test AverageWith() : Discrepancy in return vector %g\n", derr);

  // Test KernelFuncChi2() ----------------------------------------------------
  //
  // x = [ 0.0; 0.0; 0.3; 0.0 ]
  // y = [ 0.1; 0.0; 1.2; 0.1 ]
  // kval = 0.48
  double mx2[] = { 0.0, 0.0, 0.3, 0.0 };
  x.SetPoint(mx2, d);
  double my2[] = { 0.1, 0.0, 1.2, 0.1 };
  y.SetPoint(my2, d);
  double kval = x.KernelFuncChi2(y);
  double kval_truth = 0.48;
  derr = fabs(kval-kval_truth);
  printf("Test KernelFuncChi2() : Discrepancy in return value %g\n", derr);

  return 0;

}
