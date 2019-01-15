// Unit test of the class Polynomial.

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
  printf("|  Test_Polynomial.cpp                                          |\n");
  printf("+---------------------------------------------------------------+\n");

  // Test Eval() --------------------------------------------------------------
  //
  // x = [ 0.0; 0.0; 0.4; 0.3 ]
  // y = [ 0.0; 0.1; 0.0; 1.0 ]
  // s = 2.0
  // a = 0.5
  // c = 0.1
  // deg = 3
  // lambda = 0
  // z = 0.03125
  double s = 2.0;
  double a = 0.5;
  double c = 0.1;
  double deg = 3;
  Polynomial mKernel(s, a, c, deg);
  INTEGER d = 4;
  DPoint x(d), y(d);
  double mx[] = { 0.0, 0.0, 0.4, 0.3 };
  double my[] = { 0.0, 0.1, 0.0, 1.0 };
  x.SetPoint(mx, d);
  y.SetPoint(my, d);
  double z = mKernel.Eval<DPoint>(x, y);
  double z_truth = 0.03125;
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
