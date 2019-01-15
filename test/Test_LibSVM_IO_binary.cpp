// Unit test of the namespace LibSVM_IO_binary.

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
  printf("|  Test_LibSVM_IO_binary.cpp                                    |\n");
  printf("+---------------------------------------------------------------+\n");

  const char *file = "Test_LibSVM_IO_binary.dat";
  INTEGER d = 5;
  DVector y;
  bool ret;

  // Test ReadData() ----------------------------------------------------------
  //
  // DX = [ 0.0       0.0        0.00872489 0.0023827 0.0
  //        0.0102457 0.00150432 0.0        0.0016495 0.0
  //        0.0       0.0        0.0102457  0.0       0.0 ]
  // y = [ 2.0; 20.0; 1.0 ]
  DPointArray DX;
  ret = LibSVM_IO_binary::ReadData(file, DX, y, d);
  if (ret == false) {
    printf("Test ReadData() : Mysterious error returned\n");
  }

  INTEGER N = 3;
  double DX_truth[] = { 0.0, 0.0102457, 0.0,
                        0.0, 0.00150432, 0.0,
                        0.00872489, 0.0, 0.0102457,
                        0.0023827, 0.0016495, 0.0,
                        0.0, 0.0, 0.0 };
  double derr = Diff1(DX.GetPointer(), DX_truth, N*d);
  printf("Test ReadData() : Discrepancy in X %g\n", derr);

  double y_truth[] = { 2.0, 20.0, 1.0 };
  derr = Diff1(y.GetPointer(), y_truth, N);
  printf("Test ReadData() : Discrepancy in y %g\n", derr);

  return 0;

}
