// Unit test of the class NoName.

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
  printf("|  Test_NoName.cpp                                              |\n");
  printf("+---------------------------------------------------------------+\n");

  // Test Eval() --------------------------------------------------------------
  //
  // x = [ 0.4; -0.3 ]
  NoName mTestFunction;
  INTEGER d = 2;
  DPoint x(d);
  double mx[] = { 0.4, -0.3 };
  x.SetPoint(mx, d);
  double y = mTestFunction.Eval(x);
  double y_truth = -0.305158110371360;
  double derr = fabs(y - y_truth);
  printf("Test Eval() : Discrepancy in return value %g\n", derr);

  return 0;

}
