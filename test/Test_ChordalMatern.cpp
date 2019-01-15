// Unit test of the class ChordalMatern.

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
  printf("|  Test_ChordalMatern.cpp                                       |\n");
  printf("+---------------------------------------------------------------+\n");

  // Test Eval() --------------------------------------------------------------
  //
  // x = [ pi/7; 3*pi/5 ]
  // y = [ -pi/2; pi/2 ]
  // s = 2.0
  // nu = 1.3
  // ell = 2.0
  // lambda = 0
  // z = 1.107058457456521
  double s = 2.0;
  double nu = 1.3;
  double ell = 2.0;
  ChordalMatern mKernel(s, nu, ell);
  INTEGER d = 2;
  DPoint x(d), y(d);
  double mx[] = { M_PI/7.0, 3.0*M_PI/5.0 };
  double my[] = { -M_PI/2.0, M_PI/2.0 };
  x.SetPoint(mx, d);
  y.SetPoint(my, d);
  double z = mKernel.Eval<DPoint>(x, y);
  double z_truth = 1.107058457456521;
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
