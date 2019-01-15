// Unit test of the namespace Raw_binary.

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
  printf("|  Test_Raw_binary.cpp                                          |\n");
  printf("+---------------------------------------------------------------+\n");

  // Test ReadData() ----------------------------------------------------------
  //
  // y = [ 0.537667139546100;
  //       1.833885014595086;
  //      -2.258846861003648;
  //       0.862173320368121;
  //       0.318765239858981 ]
  const char *file = "Test_Raw_binary.dat";
  DVector y;
  bool ret = Raw_binary::ReadData(file, y);
  if (ret == false) {
    printf("Test ReadData() : Mysterious error returned\n");
  }

  INTEGER N = 5;
  double y_truth[] = { 0.537667139546100,
                       1.833885014595086,
                      -2.258846861003648,
                       0.862173320368121,
                       0.318765239858981 };
  double derr = Diff1(y.GetPointer(), y_truth, N);
  printf("Test ReadData() : Discrepancy in y %g\n", derr);

  return 0;

}
