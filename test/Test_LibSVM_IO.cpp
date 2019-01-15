// Unit test of the namespace LibSVM_IO.

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
  printf("|  Test_LibSVM_IO.cpp                                           |\n");
  printf("+---------------------------------------------------------------+\n");

  const char *file = "Test_LibSVM_IO.dat";
  INTEGER d = 5;
  DVector y;
  bool ret;

  // Test ReadData() for dense format -----------------------------------------
  //
  // DX = [ 0.0       0.0        0.00872489 0.0023827 0.0
  //        0.0102457 0.00150432 0.0        0.0016495 0.0
  //        0.0       0.0        0.0102457  0.0       0.0 ]
  // y = [ 2.0; 20.0; 1.0 ]
  DPointArray DX;
  ret = LibSVM_IO::ReadData(file, DX, y, d);
  if (ret == false) {
    printf("Test ReadData() for dense format : Mysterious error returned\n");
  }

  INTEGER N = 3;
  double DX_truth[] = { 0.0, 0.0102457, 0.0,
                        0.0, 0.00150432, 0.0,
                        0.00872489, 0.0, 0.0102457,
                        0.0023827, 0.0016495, 0.0,
                        0.0, 0.0, 0.0 };
  double derr = Diff1(DX.GetPointer(), DX_truth, N*d);
  printf("Test ReadData() for dense format : Discrepancy in X %g\n", derr);

  double y_truth[] = { 2.0, 20.0, 1.0 };
  derr = Diff1(y.GetPointer(), y_truth, N);
  printf("Test ReadData() for dense format : Discrepancy in y %g\n", derr);

  // Test ReadData() for sparse format ----------------------------------------
  SPointArray SX;
  ret = LibSVM_IO::ReadData(file, SX, y, d);
  if (ret == false) {
    printf("Test ReadData() for sparse format : Mysterious error returned\n");
  }

  INTEGER nnz = 6;
  INTEGER SX_truth_start[] = { 0, 2, 5, 6 };
  derr = Diff1(SX.GetPointerStart(), SX_truth_start, N+1);
  printf("Test ReadData() for sparse format : Discrepancy in start %g\n", derr);

  INTEGER SX_truth_idx[] = { 2, 3, 0, 1, 3, 2 };
  derr = Diff1(SX.GetPointerIdx(), SX_truth_idx, nnz);
  printf("Test ReadData() for sparse format : Discrepancy in idx %g\n", derr);

  double SX_truth_X[] = { 0.00872489, 0.0023827, 0.0102457,
                          0.00150432, 0.0016495, 0.0102457 };
  derr = Diff1(SX.GetPointerX(), SX_truth_X, nnz);
  printf("Test ReadData() for sparse format : Discrepancy in X %g\n", derr);

  derr = Diff1(y.GetPointer(), y_truth, N);
  printf("Test ReadData() for sparse format : Discrepancy in y %g\n", derr);

  return 0;

}
