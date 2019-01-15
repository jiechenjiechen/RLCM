// Unit test of the class IsotropicGaussain.

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
  printf("|  Test_IsotropicGaussian.cpp                                   |\n");
  printf("+---------------------------------------------------------------+\n");

  // Test Eval() --------------------------------------------------------------
  //
  // x = [ 0.0; 0.4; -0.3 ]
  // y = [ 0.1; 0.0; 1.0 ]
  // s = 2.0
  // sigma = 2.0
  // lambda = 0
  // z = 1.585099498693556
  double s = 2.0;
  double sigma = 2.0;
  IsotropicGaussian mKernel(s, sigma);
  INTEGER d = 3;
  DPoint x(d), y(d);
  double mx[] = { 0.0, 0.4, -0.3 };
  double my[] = { 0.1, 0.0, 1.0 };
  x.SetPoint(mx, d);
  y.SetPoint(my, d);
  double z = mKernel.Eval<DPoint>(x, y);
  double z_truth = 1.585099498693556;
  double derr = fabs(z - z_truth);
  printf("Test Eval() : Discrepancy in return value %g\n", derr);

  // Test Eval() --------------------------------------------------------------
  //
  // lambda = 0.1
  // Because x != y, lambda will not take effect
  double lambda = 0.1;
  z = mKernel.Eval<DPoint>(x, y, lambda);
  derr = fabs(z - z_truth);
  printf("Test Eval() : Discrepancy in return value %g\n", derr);

  return 0;

}
