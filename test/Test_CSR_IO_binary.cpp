// Unit test of the namespace CSR_IO_binary.

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
  printf("|  Test_CSR_IO_binary.cpp                                       |\n");
  printf("+---------------------------------------------------------------+\n");

  // The graph is borrowed from Figure 2(b) of the METIS manual (page
  // 11). Self loops with a constant weight 1.0 are added.
  const char *file = "Test_CSR_IO_binary.dat";
  bool ret;
  double derr;

  // Test ReadData() ----------------------------------------------------------
  //
  // A = [ 1 1 2 0 1 0 0
  //       1 1 2 1 0 0 0
  //       2 2 1 2 3 0 0
  //       0 1 2 1 0 2 5
  //       1 0 3 0 1 2 0
  //       0 0 0 2 2 1 6
  //       0 0 0 5 0 6 1 ]
  SMatrix A;
  ret = CSR_IO_binary::ReadData(file, A);
  if (ret == false) {
    printf("Test ReadData() : Mysterious error returned\n");
  }

  INTEGER N = 7, nnz = 29;
  INTEGER A_truth_start[] = { 0, 4, 8, 13, 18, 22, 26, 29 };
  derr = Diff1(A.GetPointerStart(), A_truth_start, N+1);
  printf("Test ReadData() : Discrepancy in start %g\n", derr);

  INTEGER A_truth_idx[] = { 0, 1, 2, 4,
                            0, 1, 2, 3,
                            0, 1, 2, 3, 4,
                            1, 2, 3, 5, 6,
                            0, 2, 4, 5,
                            3, 4, 5, 6,
                            3, 5, 6 };
  derr = Diff1(A.GetPointerIdx(), A_truth_idx, nnz);
  printf("Test ReadData() : Discrepancy in idx %g\n", derr);

  double A_truth_A[] = { 1, 1, 2, 1,
                         1, 1, 2, 1,
                         2, 2, 1, 2, 3,
                         1, 2, 1, 2, 5,
                         1, 3, 1, 2,
                         2, 2, 1, 6,
                         5, 6, 1 };
  derr = Diff1(A.GetPointerA(), A_truth_A, nnz);
  printf("Test ReadData() : Discrepancy in A %g\n", derr);

  return 0;

}
