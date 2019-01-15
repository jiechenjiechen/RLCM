// Unit test of the class DVector.

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
  printf("|  Test_DVector.cpp                                             |\n");
  printf("+---------------------------------------------------------------+\n");

  // a = [ 0.1; 0.4; 0.1; -0.0; -1.0; 1.0; -0.1; 0.0; 0.0; 0.1 ]
  double ma[] = { 0.1, 0.4, 0.1, -0.0, -1.0, 1.0, -0.1, 0.0, 0.0, 0.1 };
  INTEGER N = 10;
  DVector a(N);
  memcpy(a.GetPointer(), ma, N*sizeof(double));

  // Test GetEntry() ----------------------------------------------------------
  //
  // a[3] = 0.0
  INTEGER i = 3;
  double ai = a.GetEntry(i);
  double ai_truth = 0.0;
  double derr = fabs(ai - ai_truth);
  printf("Test GetEntry() : Discrepancy in return value %g\n", derr);

  // Test GetBlock() ----------------------------------------------------------
  //
  // a[3:7] = [ -0.0; -1.0; 1.0; -0.1; 0.0 ]
  INTEGER RowStart = 3;
  INTEGER nRow = 5;
  double mb[10];
  a.GetBlock(RowStart, nRow, mb);
  double mb_truth[10] = { 0.0, -1.0, 1.0, -0.1, 0.0 };
  derr = Diff1(mb, mb_truth, nRow);
  printf("Test GetBlock() : Discrepancy in return vector %g\n", derr);

  // Test GetBlock() ----------------------------------------------------------
  //
  // a[3,7,1,4] = [ -0.0; 0.0; 0.4; -1.0 ]
  INTEGER idx[10] = { 3, 7, 1, 4 };
  INTEGER n = 4;
  DVector b;
  a.GetBlock(idx, n, b);
  mb_truth[0] = 0.0; mb_truth[1] = 0.0; mb_truth[2] = 0.4; mb_truth[3] = -1.0;
  derr = Diff1(b.GetPointer(), mb_truth, n);
  printf("Test GetBlock() : Discrepancy in return vector %g\n", derr);

  // Test SetBlock() ----------------------------------------------------------
  //
  // Set a[3:7] = [ -1.0; -1.4; 0.0; 0.1; 0.8 ]
  RowStart = 3;
  nRow = 5;
  mb[0] = -1.0; mb[1] = -1.4; mb[2] = 0.0; mb[3] = 0.1; mb[4] = 0.8;
  a.SetBlock(RowStart, nRow, mb);
  double ma_truth[] = { 0.1, 0.4, 0.1, -1.0, -1.4, 0.0, 0.1, 0.8, 0.0, 0.1 };
  derr = Diff1(a.GetPointer(), ma_truth, N);
  printf("Test GetBlock() : Discrepancy in return vector %g\n", derr);

  // Test Permute() -----------------------------------------------------------
  //
  // a = [ 0.1; 0.4; 0.1; -0.0; -1.0; 1.0; -0.1; 0.0; 0.0; 0.1 ]
  // Perm = [ 2 5 3 8 4 9 7 0 1 6 ]
  // after permute, a = [ 0.1; 1.0; 0.0; 0.0; -1.0; 0.1; 0.0; 0.1; 0.4; -0.1 ]
  INTEGER Perm[] = { 2, 5, 3, 8, 4, 9, 7, 0, 1, 6 };
  memcpy(a.GetPointer(), ma, N*sizeof(double));
  a.Permute(Perm, N);
  double a_perm_truth[] = { 0.1, 1.0, 0.0, 0.0, -1.0, 0.1, 0.0, 0.1, 0.4, -0.1 };
  derr = Diff1(a.GetPointer(), a_perm_truth, N);
  printf("Test Permute() : Discrepancy in return vector %g\n", derr);

  // Test iPermute() ----------------------------------------------------------
  //
  // a = [ 0.1; 0.4; 0.1; -0.0; -1.0; 1.0; -0.1; 0.0; 0.0; 0.1 ]
  // iPerm = [ 2; 5; 3; 8; 4; 9; 7; 0; 1; 6 ]
  // after ipermute, a = [ 0.0; 0.0; 0.1; 0.1; -1.0; 0.4; 0.1; -0.1; 0.0; 1.0 ]
  INTEGER iPerm[] = { 2, 5, 3, 8, 4, 9, 7, 0, 1, 6 };
  memcpy(a.GetPointer(), ma, N*sizeof(double));
  a.iPermute(iPerm, N);
  double a_iperm_truth[] = { 0.0, 0.0, 0.1, 0.1, -1.0, 0.4, 0.1, -0.1, 0.0, 1.0 };
  derr = Diff1(a.GetPointer(), a_iperm_truth, N);
  printf("Test iPermute() : Discrepancy in return vector %g\n", derr);

  // Test Sort() --------------------------------------------------------------
  //
  // a = [ 0.1; 0.4; 0.1; -0.0; -1.0; 1.0; -0.1; 0.0; 0.0; 0.1 ]
  // after sorting (in the ascending order)
  // a = [ -1.0; -0.1; 0.0; 0.0; 0.0; 0.1; 0.1; 0.1; 0.4; 1.0 ]
  memcpy(a.GetPointer(), ma, N*sizeof(double));
  a.Sort(ASCEND);
  double a_sort_truth[] = { -1.0, -0.1, 0.0, 0.0, 0.0, 0.1, 0.1, 0.1, 0.4, 1.0 };
  derr = Diff1(a.GetPointer(), a_sort_truth, N);
  printf("Test Sort() : Discrepancy in return vector %g\n", derr);

  // Test SortByMagnitude() ---------------------------------------------------
  //
  // a = [ 0.1; 0.4; 0.1; -0.0; -1.0; 1.0; -0.1; 0.0; 0.0; 0.1 ]
  // after sorting (in the descending order)
  // a = [ -1.0; -0.1; 0.0; 0.0; 0.0; 0.1; 0.1; 0.1; 0.4; 1.0 ]
  memcpy(a.GetPointer(), ma, N*sizeof(double));
  a.SortByMagnitude(DESCEND);
  double *ptr_a = a.GetPointer();
  bool correct = true;
  for (i = 0; i < N-1; i++) {
    if (fabs(ptr_a[i]) < fabs(ptr_a[i+1])) {
      correct = false;
      break;
    }
  }
  printf("Test SortByMagnitude() : Result is %s\n", correct?"correct":"wrong");

  // Test FindLargerThan() ----------------------------------------------------
  //
  // a = [ 0.1; 0.4; 0.1; -0.0; -1.0; 1.0; -0.1; 0.0; 0.0; 0.1 ]
  // tol = 0.1
  // return idx = [ 1; 5 ]
  memcpy(a.GetPointer(), ma, N*sizeof(double));
  double tol = 0.1;
  INTEGER num = a.FindLargerThan(tol, idx);
  INTEGER num_truth = 2;
  INTEGER idx_truth[10] = { 1, 5 };
  correct = (num==num_truth && idx[0]==idx_truth[0] && idx[1]==idx_truth[1]);
  printf("Test FindLargerThan() : Result is %s\n", correct?"correct":"wrong");

  // Test Negate() ------------------------------------------------------------
  //
  // a = [ 0.1; 0.4; 0.1; -0.0; -1.0; 1.0; -0.1; 0.0; 0.0; 0.1 ]
  // b = [ -0.1; -0.4; -0.1; 0.0; 1.0; -1.0; 0.1; -0.0; -0.0; -0.1 ]
  memcpy(a.GetPointer(), ma, N*sizeof(double));
  a.Negate();
  double a_negate_truth[] = { -0.1, -0.4, -0.1, 0.0, 1.0, -1.0, 0.1, -0.0, -0.0, -0.1 };
  derr = Diff1(a.GetPointer(), a_negate_truth, N);
  printf("Test Negate() : Discrepancy in return vector %g\n", derr);
  
  memcpy(a.GetPointer(), ma, N*sizeof(double));
  a.Negate(b);
  derr = Diff1(b.GetPointer(), a_negate_truth, N);
  printf("Test Negate() : Discrepancy in return vector %g\n", derr);
  
  // Test Add() ---------------------------------------------------------------
  //
  // a = [ 0.1; 0.4; 0.1; -0.0; -1.0; 1.0; -0.1; 0.0; 0.0; 0.1 ]
  // b = 1.0
  // c = [ 1.1; 1.4; 1.1; 1.0; 0.0; 2.0; 0.9; 1.0; 1.0; 1.1 ]
  memcpy(a.GetPointer(), ma, N*sizeof(double));
  double b2 = 1.0;
  a.Add(b2);
  double a_add_truth[] = { 1.1, 1.4, 1.1, 1.0, 0.0, 2.0, 0.9, 1.0, 1.0, 1.1 };
  derr = Diff1(a.GetPointer(), a_add_truth, N);
  printf("Test Add() : Discrepancy in return vector %g\n", derr);

  memcpy(a.GetPointer(), ma, N*sizeof(double));
  DVector c;
  a.Add(b2, c);
  derr = Diff1(c.GetPointer(), a_add_truth, N);
  printf("Test Add() : Discrepancy in return vector %g\n", derr);
  
  // Test Add() ---------------------------------------------------------------
  //
  // a = [ 0.1; 0.4; 0.1; -0.0; -1.0; 1.0; -0.1; 0.0; 0.0; 0.1 ]
  // b = [ 0.4; 0.1; -0.1; 0.1; 1.0; 1.9; -0.3; 1.2; 0.0; -0.4 ]
  // c = [ 0.5; 0.5; 0.0; 0.1; 0.0; 2.9; -0.4; 1.2; 0.0; -0.3 ]
  memcpy(a.GetPointer(), ma, N*sizeof(double));
  double mb3[] = { 0.4, 0.1, -0.1, 0.1, 1.0, 1.9, -0.3, 1.2, 0.0, -0.4 };
  memcpy(b.GetPointer(), mb3, N*sizeof(double));
  a.Add(b);
  double a_add_truth2[] = { 0.5, 0.5, 0.0, 0.1, 0.0, 2.9, -0.4, 1.2, 0.0, -0.3 };
  derr = Diff1(a.GetPointer(), a_add_truth2, N);
  printf("Test Add() : Discrepancy in return vector %g\n", derr);
  
  memcpy(a.GetPointer(), ma, N*sizeof(double));
  a.Add(b, c);
  derr = Diff1(c.GetPointer(), a_add_truth2, N);
  printf("Test Add() : Discrepancy in return vector %g\n", derr);
  
  // Test Subtract() ----------------------------------------------------------
  //
  // a = [ 0.1; 0.4; 0.1; -0.0; -1.0; 1.0; -0.1; 0.0; 0.0; 0.1 ]
  // b = 1.0
  // c = [ -0.9; -0.6; -0.9; -1.0; -2.0; 0.0; -1.1; -1.0; -1.0; -0.9 ]
  memcpy(a.GetPointer(), ma, N*sizeof(double));
  b2 = 1.0;
  a.Subtract(b2);
  double a_subtract_truth[] = { -0.9, -0.6, -0.9, -1.0, -2.0, 0.0, -1.1, -1.0, -1.0, -0.9 };
  derr = Diff1(a.GetPointer(), a_subtract_truth, N);
  printf("Test Subtract() : Discrepancy in return vector %g\n", derr);

  memcpy(a.GetPointer(), ma, N*sizeof(double));
  a.Subtract(b2, c);
  derr = Diff1(c.GetPointer(), a_subtract_truth, N);
  printf("Test Subtract() : Discrepancy in return vector %g\n", derr);
  
  // Test Subtract() ----------------------------------------------------------
  //
  // a = [ 0.1; 0.4; 0.1; -0.0; -1.0; 1.0; -0.1; 0.0; 0.0; 0.1 ]
  // b = [ 0.4; 0.1; -0.1; 0.1; 1.0; 1.9; -0.3; 1.2; 0.0; -0.4 ] (same as old)
  // c = [ -0.3; 0.3; 0.2; -0.1; -2.0; -0.9; 0.2; -1.2; 0.0; 0.5 ]
  memcpy(a.GetPointer(), ma, N*sizeof(double));
  memcpy(b.GetPointer(), mb3, N*sizeof(double));
  a.Subtract(b);
  double a_subtract_truth2[] = { -0.3, 0.3, 0.2, -0.1, -2.0, -0.9, 0.2, -1.2, 0.0, 0.5 };
  derr = Diff1(a.GetPointer(), a_subtract_truth2, N);
  printf("Test Subtract() : Discrepancy in return vector %g\n", derr);
  
  memcpy(a.GetPointer(), ma, N*sizeof(double));
  a.Subtract(b, c);
  derr = Diff1(c.GetPointer(), a_subtract_truth2, N);
  printf("Test Subtract() : Discrepancy in return vector %g\n", derr);
  
  // Test Multiply() ----------------------------------------------------------
  //
  // a = [ 0.1; 0.4; 0.1; -0.0; -1.0; 1.0; -0.1; 0.0; 0.0; 0.1 ]
  // b = 2.0
  // c = [ 0.2; 0.8; 0.2; 0.0; -2.0; 2.0; -0.2; 0.0; 0.0; 0.2 ]
  memcpy(a.GetPointer(), ma, N*sizeof(double));
  b2 = 2.0;
  a.Multiply(b2);
  double a_multiply_truth[] = { 0.2, 0.8, 0.2, 0.0, -2.0, 2.0, -0.2, 0.0, 0.0, 0.2 };
  derr = Diff1(a.GetPointer(), a_multiply_truth, N);
  printf("Test Multiply() : Discrepancy in return vector %g\n", derr);

  memcpy(a.GetPointer(), ma, N*sizeof(double));
  a.Multiply(b2, c);
  derr = Diff1(c.GetPointer(), a_multiply_truth, N);
  printf("Test Multiply() : Discrepancy in return vector %g\n", derr);
  
  // Test Multiply() ----------------------------------------------------------
  //
  // a = [ 0.1; 0.4; 0.1; -0.0; -1.0; 1.0; -0.1; 0.0; 0.0; 0.1 ]
  // b = [ 0.4; 0.1; -0.1; 0.1; 1.0; 1.9; -0.3; 1.2; 0.0; -0.4 ] (same as old)
  // c = [ 0.04; 0.04; -0.01; 0.0; -1.0; 1.9; 0.03; 0.0; 0.0; -0.04 ]
  memcpy(a.GetPointer(), ma, N*sizeof(double));
  memcpy(b.GetPointer(), mb3, N*sizeof(double));
  a.Multiply(b);
  double a_multiply_truth2[] = { 0.04, 0.04, -0.01, 0.0, -1.0, 1.9, 0.03, 0.0, 0.0, -0.04 };
  derr = Diff1(a.GetPointer(), a_multiply_truth2, N);
  printf("Test Multiply() : Discrepancy in return vector %g\n", derr);
  
  memcpy(a.GetPointer(), ma, N*sizeof(double));
  a.Multiply(b, c);
  derr = Diff1(c.GetPointer(), a_multiply_truth2, N);
  printf("Test Multiply() : Discrepancy in return vector %g\n", derr);
  
  // Test Divide() ------------------------------------------------------------
  //
  // a = [ 0.1; 0.4; 0.1; -0.0; -1.0; 1.0; -0.1; 0.0; 0.0; 0.1 ]
  // b = 2.0
  // c = [ 0.05; 0.2; 0.05; 0.0; -0.5; 0.5; -0.05; 0.0; 0.0; 0.05 ]
  memcpy(a.GetPointer(), ma, N*sizeof(double));
  b2 = 2.0;
  a.Divide(b2);
  double a_divide_truth[] = { 0.05, 0.2, 0.05, 0.0, -0.5, 0.5, -0.05, 0.0, 0.0, 0.05 };
  derr = Diff1(a.GetPointer(), a_divide_truth, N);
  printf("Test Divide() : Discrepancy in return vector %g\n", derr);

  memcpy(a.GetPointer(), ma, N*sizeof(double));
  a.Divide(b2, c);
  derr = Diff1(c.GetPointer(), a_divide_truth, N);
  printf("Test Divide() : Discrepancy in return vector %g\n", derr);
  
  // Test Divide() ------------------------------------------------------------
  //
  // a = [ 0.1; 0.4; 0.1; -0.0; -1.0; 1.0; -0.1; 0.0; 0.0; 0.1 ]
  // b = [ 0.4; 0.1; -0.1; 0.1; 1.0; 1.9; -0.3; 1.2; 0.8; -0.4 ] (b[8] changed)
  // c = [ 1.0/4; 4; -1; 0; -1; 1.0/1.9; 1.0/3; 0; 0; -1.0/4 ]
  memcpy(a.GetPointer(), ma, N*sizeof(double));
  memcpy(b.GetPointer(), mb3, N*sizeof(double));
  b.SetEntry(8, 0.8);
  a.Divide(b);
  double a_divide_truth2[] = { 1.0/4, 4, -1, 0, -1, 1.0/1.9, 1.0/3, 0, 0, -1.0/4 };
  derr = Diff1(a.GetPointer(), a_divide_truth2, N);
  printf("Test Divide() : Discrepancy in return vector %g\n", derr);
  
  memcpy(a.GetPointer(), ma, N*sizeof(double));
  a.Divide(b, c);
  derr = Diff1(c.GetPointer(), a_divide_truth2, N);
  printf("Test Divide() : Discrepancy in return vector %g\n", derr);
  
  // Test InProd() ------------------------------------------------------------
  //
  // a = [ 0.1; 0.4; 0.1; -0.0; -1.0; 1.0; -0.1; 0.0; 0.0; 0.1 ]
  // b = [ 0.4; 0.1; -0.1; 0.1; 1.0; 1.9; -0.3; 1.2; 0.0; -0.4 ] (same as old)
  // inprod = 0.96
  memcpy(a.GetPointer(), ma, N*sizeof(double));
  memcpy(b.GetPointer(), mb3, N*sizeof(double));
  double inp = a.InProd(b);
  double inp_truth = 0.96;
  derr = fabs(inp - inp_truth);
  printf("Test InProd() : Discrepancy in return value %g\n", derr);
  
  // Test Abs() ---------------------------------------------------------------
  //
  // a = [ 0.1; 0.4; 0.1; -0.0; -1.0; 1.0; -0.1; 0.0; 0.0; 0.1 ]
  // b = [ 0.1; 0.4; 0.1; 0.0; 1.0; 1.0; 0.1; 0.0; 0.0; 0.1 ]
  memcpy(a.GetPointer(), ma, N*sizeof(double));
  a.Abs();
  double a_abs_truth[] = { 0.1, 0.4, 0.1, 0.0, 1.0, 1.0, 0.1, 0.0, 0.0, 0.1 };
  derr = Diff1(a.GetPointer(), a_abs_truth, N);
  printf("Test Abs() : Discrepancy in return vector %g\n", derr);
  
  memcpy(a.GetPointer(), ma, N*sizeof(double));
  a.Abs(b);
  derr = Diff1(b.GetPointer(), a_abs_truth, N);
  printf("Test Abs() : Discrepancy in return vector %g\n", derr);
  
  // Test Sqrt() --------------------------------------------------------------
  //
  // a = [ 0.01; 0.16; 0.01; 0.0; 1.21; 1.21; 0.01; 0.0; 0.0; 0.01 ]
  // b = [ 0.1; 0.4; 0.1; 0.0; 1.1; 1.1; 0.1; 0.0; 0.0; 0.1 ]
  double ma2[] = { 0.01, 0.16, 0.01, 0.0, 1.21, 1.21, 0.01, 0.0, 0.0, 0.01 };
  memcpy(a.GetPointer(), ma2, N*sizeof(double));
  a.Sqrt();
  double a_sqrt_truth[] = { 0.1, 0.4, 0.1, 0.0, 1.1, 1.1, 0.1, 0.0, 0.0, 0.1 };
  derr = Diff1(a.GetPointer(), a_sqrt_truth, N);
  printf("Test Sqrt() : Discrepancy in return vector %g\n", derr);
  
  memcpy(a.GetPointer(), ma2, N*sizeof(double));
  a.Sqrt(b);
  derr = Diff1(b.GetPointer(), a_sqrt_truth, N);
  printf("Test Sqrt() : Discrepancy in return vector %g\n", derr);
  
  // Test Square() ------------------------------------------------------------
  //
  // Reverse the role of a and b in the above
  memcpy(a.GetPointer(), a_sqrt_truth, N*sizeof(double));
  a.Square();
  derr = Diff1(a.GetPointer(), ma2, N);
  printf("Test Square() : Discrepancy in return vector %g\n", derr);
  
  memcpy(a.GetPointer(), a_sqrt_truth, N*sizeof(double));
  a.Square(b);
  derr = Diff1(b.GetPointer(), ma2, N);
  printf("Test Square() : Discrepancy in return vector %g\n", derr);
  
  // Test Inv() ---------------------------------------------------------------
  //
  // a = [ 0.1; 0.4; 0.1; -2.0; -1.0; 1.0; -0.1; 0.9; 2.0; 0.1 ]
  // b = [ 10; 2.5; 10; -0.5; -1.0; 1.0; -10.0; 10.0/9; 0.5; 10.0 ]
  double ma3[] = { 0.1, 0.4, 0.1, -2.0, -1.0, 1.0, -0.1, 0.9, 2.0, 0.1 };
  memcpy(a.GetPointer(), ma3, N*sizeof(double));
  a.Inv();
  double a_inv_truth[] = { 10, 2.5, 10, -0.5, -1.0, 1.0, -10.0, 10.0/9, 0.5, 10.0 };
  derr = Diff1(a.GetPointer(), a_inv_truth, N);
  printf("Test Inv() : Discrepancy in return vector %g\n", derr);
  
  memcpy(a.GetPointer(), ma3, N*sizeof(double));
  a.Inv(b);
  derr = Diff1(b.GetPointer(), a_inv_truth, N);
  printf("Test Inv() : Discrepancy in return vector %g\n", derr);
  
  // Test Norm2() -------------------------------------------------------------
  //
  // a = [ 0.1; 0.4; 0.1; -2.0; -1.0; 1.0; -0.1; 0.9; 2.0; 0.1 ]
  // norm2 = 3.318132004607412
  memcpy(a.GetPointer(), ma3, N*sizeof(double));
  double norm2 = a.Norm2();
  double norm2_truth = 3.318132004607412;
  derr = fabs(norm2 - norm2_truth);
  printf("Test Norm2() : Discrepancy in return value %g\n", derr);
  
  // Test Norm1() -------------------------------------------------------------
  //
  // a = [ 0.1; 0.4; 0.1; -2.0; -1.0; 1.0; -0.1; 0.9; 2.0; 0.1 ]
  // norm1 = 7.7
  memcpy(a.GetPointer(), ma3, N*sizeof(double));
  double norm1 = a.Norm1();
  double norm1_truth = 7.7;
  derr = fabs(norm1 - norm1_truth);
  printf("Test Norm1() : Discrepancy in return value %g\n", derr);
  
  // Test Min() ---------------------------------------------------------------
  //
  // a = [ 0.1; 0.4; 0.1; -2.0; -1.0; 1.0; -0.1; 0.9; 2.0; 0.1 ]
  // min = -2.0, idx = 3
  memcpy(a.GetPointer(), ma3, N*sizeof(double));
  INTEGER idx2;
  double min = a.Min(idx2);
  double min_truth = -2.0;
  INTEGER idx2_truth = 3;
  derr = fabs(min - min_truth);
  correct = (min==min_truth && idx2==idx2_truth);
  printf("Test Min() : Result is %s\n", correct?"correct":"wrong");
  
  // Test Max() ---------------------------------------------------------------
  //
  // a = [ 0.1; 0.4; 0.1; -2.0; -1.0; 1.0; -0.1; 0.9; 2.0; 0.1 ]
  // max = 2.0, idx = 8
  memcpy(a.GetPointer(), ma3, N*sizeof(double));
  double max = a.Max(idx2);
  double max_truth = 2.0;
  idx2_truth = 8;
  derr = fabs(max - max_truth);
  correct = (max==max_truth && idx2==idx2_truth);
  printf("Test Max() : Result is %s\n", correct?"correct":"wrong");
  
  // Test Sum() ---------------------------------------------------------------
  //
  // a = [ 0.1; 0.4; 0.1; -2.0; -1.0; 1.0; -0.1; 0.9; 2.0; 0.1 ]
  // sum = 1.5
  memcpy(a.GetPointer(), ma3, N*sizeof(double));
  double sum = a.Sum();
  double sum_truth = 1.5;
  derr = fabs(sum - sum_truth);
  printf("Test Sum() : Discrepancy in return value %g\n", derr);
  
  // Test Mean() --------------------------------------------------------------
  //
  // a = [ 0.1; 0.4; 0.1; -2.0; -1.0; 1.0; -0.1; 0.9; 2.0; 0.1 ]
  // mean = 0.15
  memcpy(a.GetPointer(), ma3, N*sizeof(double));
  double mean = a.Mean();
  double mean_truth = 0.15;
  derr = fabs(mean - mean_truth);
  printf("Test Mean() : Discrepancy in return value %g\n", derr);
  
  // Test BuildResponseVector() -----------------------------------------------
  //
  //     [ -1.06  -0.2  ]        [  0.924169173333334 ]
  // X = [  0.7    1.2  ],   y = [  2.241168000000000 ]
  //     [  0.75   0.07 ]        [ -0.912325381783667 ]
  //
  // Each row of X is a point; but the matrix is stored in column major order.
  // Test function: Six humps.
  INTEGER NX = 3;
  INTEGER DIM = 2;
  double mPX[] = { -1.06, 0.7, 0.75, -0.2, 1.2, 0.07 };
  DPointArray PX(NX, DIM);
  memcpy(PX.GetPointer(), mPX, NX*DIM*sizeof(double));

  SixHumps mTestFunction;

  DVector y;
  y.BuildResponseVector<SixHumps, DPoint, DPointArray>(mTestFunction, PX);
  double y_truth[] = { 0.924169173333334,
                       2.241168000000000,
                      -0.912325381783667 };
  derr = Diff1(y.GetPointer(), y_truth, NX);
  printf("Test BuildResponseVector() : Discrepancy in return vector %g\n", derr);

  return 0;

}
