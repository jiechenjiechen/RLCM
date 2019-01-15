// Unit test of the class DMatrix.

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
  printf("|  Test_DMatrix.cpp                                             |\n");
  printf("+---------------------------------------------------------------+\n");

  // Test GetEntry() ----------------------------------------------------------
  //
  //     [  1.0  0.0  0.3  0.4 ]
  // A = [  0.0  0.0  1.0  0.2 ]
  //     [ -0.1  1.0  0.0  0.0 ]
  //
  // A(1,2) = 1.0
  double mA[] = { 1.0, 0.0, -0.1, 0.0, 0.0, 1.0, 0.3, 1.0, 0.0, 0.4, 0.2, 0.0 };
  INTEGER M = 3, N = 4;
  DMatrix A(M, N);
  memcpy(A.GetPointer(), mA, M*N*sizeof(double));
  INTEGER i = 1, j = 2;
  double aij = A.GetEntry(i, j);
  double aij_truth = 1.0;
  double derr = fabs(aij - aij_truth);
  printf("Test GetEntry() : Discrepancy in return value %g\n", derr);

  // Test GetColumn() ---------------------------------------------------------
  //
  // Reuse A
  //
  // A(:,2) = [ 0.3; 1.0; 0.0 ]
  DVector b;
  A.GetColumn(j, b);
  double A_col_truth[] = { 0.3, 1.0, 0.0 };
  derr = Diff1(b.GetPointer(), A_col_truth, 3);
  printf("Test GetColumn() : Discrepancy in return vector %g\n", derr);

  // Test GetColumns() --------------------------------------------------------
  //
  // Reuse A
  //
  //              [ 0.3  0.0 ]
  // A(:,[2,1]) = [ 1.0  0.0 ]
  //              [ 0.0  1.0 ]
  DMatrix B;
  INTEGER IdxCol[] = { 2, 1 };
  INTEGER nCol = 2;
  A.GetColumns(IdxCol, nCol, B);
  double A_cols_truth[] = { 0.3, 1.0, 0.0, 0.0, 0.0, 1.0 };
  derr = Diff1(B.GetPointer(), A_cols_truth, 6);
  printf("Test GetColumns() : Discrepancy in return matrix %g\n", derr);

  // Test GetRow() ------------------------------------------------------------
  //
  // Reuse A
  //
  // A(1,:) = [ 0.0  0.0  1.0  0.2 ]
  A.GetRow(i, b);
  double A_row_truth[] = { 0.0, 0.0, 1.0, 0.2 };
  derr = Diff1(b.GetPointer(), A_row_truth, 4);
  printf("Test GetRow() : Discrepancy in return vector %g\n", derr);

  // Test GetBlock() ----------------------------------------------------------
  //
  // Reuse A
  //
  // A(1:2,1:3) = [ 0.0  1.0  0.2 ]
  //              [ 1.0  0.0  0.0 ]
  INTEGER RowStart = 1;
  INTEGER nRow = 2;
  INTEGER ColStart = 1;
  nCol = 3;
  A.GetBlock(RowStart, nRow, ColStart, nCol, B);
  double A_block_truth[] = { 0.0, 1.0, 1.0, 0.0, 0.2, 0.0 };
  derr = Diff1(B.GetPointer(), A_block_truth, 6);
  printf("Test GetBlock() : Discrepancy in return matrix %g\n", derr);

  // Test SetEntry() ----------------------------------------------------------
  //
  // Set A(1,2) = -1.2
  memcpy(A.GetPointer(), mA, M*N*sizeof(double));
  A.SetEntry(i, j, -1.2);
  double A_set_entry_truth[] = { 1.0, 0.0, -0.1, 0.0, 0.0, 1.0, 0.3, -1.2, 0.0, 0.4, 0.2, 0.0 };
  derr = Diff1(A.GetPointer(), A_set_entry_truth, M*N);
  printf("Test SetEntry() : Discrepancy in return matrix %g\n", derr);

  // Test AddToEntry() --------------------------------------------------------
  //
  // A(1,2) + (-1.2) = -0.2
  memcpy(A.GetPointer(), mA, M*N*sizeof(double));
  A.AddToEntry(i, j, -1.2);
  double A_add_entry_truth[] = { 1.0, 0.0, -0.1, 0.0, 0.0, 1.0, 0.3, -0.2, 0.0, 0.4, 0.2, 0.0 };
  derr = Diff1(A.GetPointer(), A_add_entry_truth, M*N);
  printf("Test AddToEntry() : Discrepancy in return matrix %g\n", derr);

  // Test SetColumn() ---------------------------------------------------------
  //
  // Reuse A
  //
  // Set A(:,2) = [ -3.3; 1.2; -0.1 ]
  memcpy(A.GetPointer(), mA, M*N*sizeof(double));
  DVector new_col(M);
  double m_new_col[] = { -3.3, 1.2, -0.1 };
  memcpy(new_col.GetPointer(), m_new_col, M*sizeof(double));
  A.SetColumn(j, new_col);
  double A_set_col_truth[] = { 1.0, 0.0, -0.1, 0.0, 0.0, 1.0, -3.3, 1.2, -0.1, 0.4, 0.2, 0.0 };
  derr = Diff1(A.GetPointer(), A_set_col_truth, M*N);
  printf("Test SetColumn() : Discrepancy in return matrix %g\n", derr);

  // Test SetRow() ------------------------------------------------------------
  //
  // Reuse A
  //
  // Set A(1,:) = [ 1.0  0.4  -1.0  2.1 ]
  memcpy(A.GetPointer(), mA, M*N*sizeof(double));
  DVector new_row(N);
  double m_new_row[] = { 1.0, 0.4, -1.0, 2.1 };
  memcpy(new_row.GetPointer(), m_new_row, N*sizeof(double));
  A.SetRow(i, new_row);
  double A_set_row_truth[] = { 1.0, 1.0, -0.1, 0.0, 0.4, 1.0, 0.3, -1.0, 0.0, 0.4, 2.1, 0.0 };
  derr = Diff1(A.GetPointer(), A_set_row_truth, M*N);
  printf("Test SetRow() : Discrepancy in return matrix %g\n", derr);

  // Test SetBlock() ----------------------------------------------------------
  //
  // Reuse A
  //
  // Set A(1:2,1:3) = [ 0.6  -1.0  -0.2 ]
  //                  [ 1.9   0.7   1.0 ]
  memcpy(A.GetPointer(), mA, M*N*sizeof(double));
  DMatrix new_block(nRow, nCol);
  double m_new_block[] = { 0.6, 1.9, -1.0, 0.7, -0.2, 1.0 };
  memcpy(new_block.GetPointer(), m_new_block, nRow*nCol*sizeof(double));
  A.SetBlock(RowStart, nRow, ColStart, nCol, new_block);
  double A_set_block_truth[] = { 1.0, 0.0, -0.1, 0.0, 0.6, 1.9, 0.3, -1.0, 0.7, 0.4, -0.2, 1.0 };
  derr = Diff1(A.GetPointer(), A_set_block_truth, M*N);
  printf("Test SetBlock() : Discrepancy in return matrix %g\n", derr);

  // Test SetIdentity() -------------------------------------------------------
  //
  // A2 is 3*3
  INTEGER N2 = 3;
  DMatrix A2(N2);
  A2.SetIdentity();
  double id_truth[] = { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };
  derr = Diff1(A2.GetPointer(), id_truth, N2*N2);
  printf("Test SetIdentity() : Discrepancy in return matrix %g\n", derr);

  // Test SetMultIdentity() ---------------------------------------------------
  //
  // A2 = lambda * eye(3), lambda = 2.3
  double lambda = 2.3;
  A2.SetMultIdentity(lambda);
  double mult_id_truth[] = { 2.3, 0.0, 0.0, 0.0, 2.3, 0.0, 0.0, 0.0, 2.3 };
  derr = Diff1(A2.GetPointer(), mult_id_truth, N2*N2);
  printf("Test SetMultIdentity() : Discrepancy in return matrix %g\n", derr);

  // Test SetConstVal() -------------------------------------------------------
  //
  // A2 = lambda * ones(3), lambda = 2.3
  A2.SetConstVal(lambda);
  double const_val_truth[] = { 2.3, 2.3, 2.3, 2.3, 2.3, 2.3, 2.3, 2.3, 2.3 };
  derr = Diff1(A2.GetPointer(), const_val_truth, N2*N2);
  printf("Test SetConstVal() : Discrepancy in return matrix %g\n", derr);

  // Test MakeDiag() ----------------------------------------------------------
  //
  // A2 = diag([ 1.0  -3.1  0.2 ])
  double md[] = { 1.0, -3.1, 0.2 };
  DVector d(N2);
  memcpy(d.GetPointer(), md, N2*sizeof(double));
  A2.MakeDiag(d);
  double diag_truth[] = { 1.0, 0.0, 0.0, 0.0, -3.1, 0.0, 0.0, 0.0, 0.2 };
  derr = Diff1(A2.GetPointer(), diag_truth, N2*N2);
  printf("Test MakeDiag() : Discrepancy in return matrix %g\n", derr);

  // Test Transpose() ---------------------------------------------------------
  //
  //           [  1.0  0.0  0.3  0.4 ]
  // Reuse A = [  0.0  0.0  1.0  0.2 ]
  //           [ -0.1  1.0  0.0  0.0 ]
  memcpy(A.GetPointer(), mA, M*N*sizeof(double));
  A.Transpose();
  double trans_truth[] = { 1.0, 0.0, 0.3, 0.4, 0.0, 0.0, 1.0, 0.2, -0.1, 1.0, 0.0, 0.0 };
  derr = Diff1(A.GetPointer(), trans_truth, M*N);
  printf("Test Transpose() : Discrepancy in return matrix %g\n", derr);

  // Test Transpose() ---------------------------------------------------------
  //
  // Reuse A
  A.Init(M, N);
  memcpy(A.GetPointer(), mA, M*N*sizeof(double));
  DMatrix AT;
  A.Transpose(AT);
  derr = Diff1(AT.GetPointer(), trans_truth, M*N);
  printf("Test Transpose() : Discrepancy in return matrix %g\n", derr);

  // Test Symmetrize() --------------------------------------------------------
  //
  //                [ 0.3  -1.0   0.2 ]
  // Set a new A2 = [ 0.4   0.3   0.0 ]
  //                [ 0.0  -1.0  -2.1 ]
  double mA2[] = { 0.3, 0.4, 0.0, -1.0, 0.3, -1.0, 0.2, 0.0, -2.1 };
  memcpy(A2.GetPointer(), mA2, N2*N2*sizeof(double));
  A2.Symmetrize();
  double sym_truth[] = { 0.3, -0.3, 0.1, -0.3, 0.3, -0.5, 0.1, -0.5, -2.1 };
  derr = Diff1(A2.GetPointer(), sym_truth, N2*N2);
  printf("Test Symmetrize() : Discrepancy in return matrix %g\n", derr);

  // Test Symmetrize() --------------------------------------------------------
  //
  // Reuse above A2
  memcpy(A2.GetPointer(), mA2, N2*N2*sizeof(double));
  DMatrix A2S;
  A2.Symmetrize(A2S);
  derr = Diff1(A2S.GetPointer(), sym_truth, N2*N2);
  printf("Test Symmetrize() : Discrepancy in return matrix %g\n", derr);

  // Test Negate() ------------------------------------------------------------
  //
  //           [  1.0  0.0  0.3  0.4 ]
  // Reuse A = [  0.0  0.0  1.0  0.2 ]
  //           [ -0.1  1.0  0.0  0.0 ]
  memcpy(A.GetPointer(), mA, M*N*sizeof(double));
  A.Negate();
  double negate_truth[] = { -1.0, -0.0, 0.1, -0.0, -0.0, -1.0, -0.3, -1.0, -0.0, -0.4, -0.2, -0.0 };
  derr = Diff1(A.GetPointer(), negate_truth, M*N);
  printf("Test Negate() : Discrepancy in return matrix %g\n", derr);

  // Test Negate() ------------------------------------------------------------
  //
  // Reuse A
  memcpy(A.GetPointer(), mA, M*N*sizeof(double));
  A.Negate(B);
  derr = Diff1(B.GetPointer(), negate_truth, M*N);
  printf("Test Negate() : Discrepancy in return matrix %g\n", derr);

  // Test Add() ---------------------------------------------------------------
  //
  //           [  1.0  0.0  0.3  0.4 ]        [ -0.5   0.3  -2.0   1.1 ]
  // Reuse A = [  0.0  0.0  1.0  0.2 ],   B = [  0.4  -3.0   1.3  -0.3 ]
  //           [ -0.1  1.0  0.0  0.0 ]        [ -0.6   1.7   2.9   1.1 ]
  //
  //     [  0.5   0.3  -1.7   1.5 ]
  // C = [  0.4  -3.0   2.3  -0.1 ]
  //     [ -0.7   2.7   2.9   1.1 ]
  memcpy(A.GetPointer(), mA, M*N*sizeof(double));
  double mB[] = { -0.5, 0.4, -0.6, 0.3, -3.0, 1.7, -2.0, 1.3, 2.9, 1.1, -0.3, 1.1 };
  // B already has a compatible size because of Negate()
  memcpy(B.GetPointer(), mB, M*N*sizeof(double));
  A.Add(B);
  double add_truth[] = { 0.5, 0.4, -0.7, 0.3, -3.0, 2.7, -1.7, 2.3, 2.9, 1.5, -0.1, 1.1 };
  derr = Diff1(A.GetPointer(), add_truth, M*N);
  printf("Test Add() : Discrepancy in return matrix %g\n", derr);

  // Test Add() ---------------------------------------------------------------
  //
  // Reuse A and B
  memcpy(A.GetPointer(), mA, M*N*sizeof(double));
  DMatrix C;
  A.Add(B, C);
  derr = Diff1(C.GetPointer(), add_truth, M*N);
  printf("Test Add() : Discrepancy in return matrix %g\n", derr);

  // Test Subtract() ----------------------------------------------------------
  //
  //           [  1.0  0.0  0.3  0.4 ]         [ -0.5   0.3  -2.0   1.1 ]
  // Reuse A = [  0.0  0.0  1.0  0.2 ] and B = [  0.4  -3.0   1.3  -0.3 ]
  //           [ -0.1  1.0  0.0  0.0 ]         [ -0.6   1.7   2.9   1.1 ]
  //
  //     [  1.5  -0.3   2.3  -0.7 ]
  // C = [ -0.4   3.0  -0.3   0.5 ]
  //     [  0.5  -0.7  -2.9  -1.1 ]
  memcpy(A.GetPointer(), mA, M*N*sizeof(double));
  A.Subtract(B);
  double subtract_truth[] = { 1.5, -0.4, 0.5, -0.3, 3.0, -0.7, 2.3, -0.3, -2.9, -0.7, 0.5, -1.1 };
  derr = Diff1(A.GetPointer(), subtract_truth, M*N);
  printf("Test Subtract() : Discrepancy in return matrix %g\n", derr);

  // Test Subtract() ----------------------------------------------------------
  //
  // Reuse A and B
  memcpy(A.GetPointer(), mA, M*N*sizeof(double));
  A.Subtract(B, C);
  derr = Diff1(C.GetPointer(), subtract_truth, M*N);
  printf("Test Subtract() : Discrepancy in return matrix %g\n", derr);

  // Test Add() ---------------------------------------------------------------
  //
  //           [  1.0  0.0  0.3  0.4 ]
  // Reuse A = [  0.0  0.0  1.0  0.2 ],   b = -2.0
  //           [ -0.1  1.0  0.0  0.0 ]
  //
  //     [ -1.0  -2.0  -1.7  -1.6 ]
  // C = [ -2.0  -2.0  -1.0  -1.8 ]
  //     [ -2.1  -1.0  -2.0  -2.0 ]
  memcpy(A.GetPointer(), mA, M*N*sizeof(double));
  double bval = -2.0;
  A.Add(bval);
  double add_truth2[] = { -1.0, -2.0, -2.1, -2.0, -2.0, -1.0, -1.7, -1.0, -2.0, -1.6, -1.8, -2.0 };
  derr = Diff1(A.GetPointer(), add_truth2, M*N);
  printf("Test Add() : Discrepancy in return matrix %g\n", derr);

  // Test Add() ---------------------------------------------------------------
  //
  // Reuse A and b
  memcpy(A.GetPointer(), mA, M*N*sizeof(double));
  A.Add(bval, C);
  derr = Diff1(C.GetPointer(), add_truth2, M*N);
  printf("Test Add() : Discrepancy in return matrix %g\n", derr);

  // Test Subtract() ----------------------------------------------------------
  //
  //           [  1.0  0.0  0.3  0.4 ]
  // Reuse A = [  0.0  0.0  1.0  0.2 ] and b = -2.0
  //           [ -0.1  1.0  0.0  0.0 ]
  //
  //     [  3.0  2.0  2.3  2.4 ]
  // C = [  2.0  2.0  3.0  2.2 ]
  //     [  1.9  3.0  2.0  2.0 ]
  memcpy(A.GetPointer(), mA, M*N*sizeof(double));
  A.Subtract(bval);
  double subtract_truth2[] = { 3.0, 2.0, 1.9, 2.0, 2.0, 3.0, 2.3, 3.0, 2.0, 2.4, 2.2, 2.0 };
  derr = Diff1(A.GetPointer(), subtract_truth2, M*N);
  printf("Test Subtract() : Discrepancy in return matrix %g\n", derr);

  // Test Subtract() ----------------------------------------------------------
  //
  // Reuse A and b
  memcpy(A.GetPointer(), mA, M*N*sizeof(double));
  A.Subtract(bval, C);
  derr = Diff1(C.GetPointer(), subtract_truth2, M*N);
  printf("Test Subtract() : Discrepancy in return matrix %g\n", derr);

  // Test Multiply() ----------------------------------------------------------
  //
  //           [  1.0  0.0  0.3  0.4 ]        [ -0.5   0.3  -2.0   1.1 ]
  // Reuse A = [  0.0  0.0  1.0  0.2 ],   B = [  0.4  -3.0   1.3  -0.3 ]
  //           [ -0.1  1.0  0.0  0.0 ]        [ -0.6   1.7   2.9   1.1 ]
  //
  // C =
  // [ -0.50  0.00 -0.60  0.44 ]
  // [  0.00  0.00  1.30 -0.06 ]
  // [  0.06  1.70  0.00  0.00 ]
  memcpy(A.GetPointer(), mA, M*N*sizeof(double));
  A.Multiply(B);
  double multiply_truth[] = { -0.50, 0.00, 0.06, 0.00, 0.00, 1.70, -0.60, 1.30, 0.00, 0.44, -0.06, 0.00 };
  derr = Diff1(A.GetPointer(), multiply_truth, M*N);
  printf("Test Multiply() : Discrepancy in return matrix %g\n", derr);

  // Test Multiply() ----------------------------------------------------------
  //
  // Reuse A and B
  memcpy(A.GetPointer(), mA, M*N*sizeof(double));
  A.Multiply(B, C);
  derr = Diff1(C.GetPointer(), multiply_truth, M*N);
  printf("Test Multiply() : Discrepancy in return matrix %g\n", derr);

  // Test Divide() ------------------------------------------------------------
  //
  //           [  1.0  0.0  0.3  0.4 ]        [ -0.5   0.3  -2.0   1.1 ]
  // Reuse A = [  0.0  0.0  1.0  0.2 ],   B = [  0.4  -3.0   1.3  -0.3 ]
  //           [ -0.1  1.0  0.0  0.0 ]        [ -0.6   1.7   2.9   1.1 ]
  //
  // C =
  //[-2.000000000000000 0.000000000000000 -0.150000000000000  0.363636363636364]
  //[ 0.000000000000000 0.000000000000000  0.769230769230769 -0.666666666666667]
  //[ 0.166666666666667 0.588235294117647  0.000000000000000  0.000000000000000]
  memcpy(A.GetPointer(), mA, M*N*sizeof(double));
  A.Divide(B);
  double divide_truth[] = { -2.000000000000000,
                             0.000000000000000,
                             0.166666666666667,
                             0.000000000000000,
                             0.000000000000000,
                             0.588235294117647,
                            -0.150000000000000,
                             0.769230769230769,
                             0.000000000000000,
                             0.363636363636364,
                            -0.666666666666667,
                             0.000000000000000 };
  derr = Diff1(A.GetPointer(), divide_truth, M*N);
  printf("Test Divide() : Discrepancy in return matrix %g\n", derr);

  // Test Divide() ------------------------------------------------------------
  //
  // Reuse A and B
  memcpy(A.GetPointer(), mA, M*N*sizeof(double));
  A.Divide(B, C);
  derr = Diff1(C.GetPointer(), divide_truth, M*N);
  printf("Test Divide() : Discrepancy in return matrix %g\n", derr);

  // Test Multiply() ----------------------------------------------------------
  //
  //           [  1.0  0.0  0.3  0.4 ]
  // Reuse A = [  0.0  0.0  1.0  0.2 ] and b = -2.0
  //           [ -0.1  1.0  0.0  0.0 ]
  //
  //     [ -2.0   0.0  -0.6  -0.8 ]
  // C = [  0.0   0.0  -2.0  -0.4 ]
  //     [  0.2  -2.0   0.0   0.0 ]
  memcpy(A.GetPointer(), mA, M*N*sizeof(double));
  A.Multiply(bval);
  double multiply_truth2[] = { -2.0, 0.0, 0.2, 0.0, 0.0, -2.0, -0.6, -2.0, 0.0, -0.8, -0.4, 0.0 };
  derr = Diff1(A.GetPointer(), multiply_truth2, M*N);
  printf("Test Multiply() : Discrepancy in return matrix %g\n", derr);

  // Test Multiply() ----------------------------------------------------------
  //
  // Reuse A and b
  memcpy(A.GetPointer(), mA, M*N*sizeof(double));
  A.Multiply(bval, C);
  derr = Diff1(C.GetPointer(), multiply_truth2, M*N);
  printf("Test Multiply() : Discrepancy in return matrix %g\n", derr);

  // Test Divide() ------------------------------------------------------------
  //
  //           [  1.0  0.0  0.3  0.4 ]
  // Reuse A = [  0.0  0.0  1.0  0.2 ] and b = -2.0
  //           [ -0.1  1.0  0.0  0.0 ]
  //
  //     [ -0.5    0.0  -0.15  -0.2 ]
  // C = [  0.0    0.0  -0.5   -0.1 ]
  //     [  0.05  -0.5   0.0    0.0 ]
  memcpy(A.GetPointer(), mA, M*N*sizeof(double));
  A.Divide(bval);
  double divide_truth2[] = { -0.5, 0.0, 0.05, 0.0, 0.0, -0.5, -0.15, -0.5, 0.0, -0.2, -0.1, 0.0 };
  derr = Diff1(A.GetPointer(), divide_truth2, M*N);
  printf("Test Divide() : Discrepancy in return matrix %g\n", derr);

  // Test Divide() ------------------------------------------------------------
  //
  // Reuse A and b
  memcpy(A.GetPointer(), mA, M*N*sizeof(double));
  A.Divide(bval, C);
  derr = Diff1(C.GetPointer(), divide_truth2, M*N);
  printf("Test Divide() : Discrepancy in return matrix %g\n", derr);

  // Test AddDiagonal() -------------------------------------------------------
  //
  //            [ 0.3  -1.0   0.2 ]
  // Reuse A2 = [ 0.4   0.3   0.0 ] and lambda = 2.3
  //            [ 0.0  -1.0  -2.1 ]
  memcpy(A2.GetPointer(), mA2, N2*N2*sizeof(double));
  A2.AddDiagonal(lambda);
  double add_diag_truth[] = { 2.6, 0.4, 0.0, -1.0, 2.6, -1.0, 0.2, 0.0, 0.2 };
  derr = Diff1(A2.GetPointer(), add_diag_truth, N2*N2);
  printf("Test AddDiagonal() : Discrepancy in return matrix %g\n", derr);

  // Test AddDiagonal() --------------------------------------------------
  //
  // Reuse A2 and lambda
  memcpy(A2.GetPointer(), mA2, N2*N2*sizeof(double));
  DMatrix B2;
  A2.AddDiagonal(lambda, B2);
  derr = Diff1(B2.GetPointer(), add_diag_truth, N2*N2);
  printf("Test AddDiagonal() : Discrepancy in return matrix %g\n", derr);

  // Test SubtractDiagonal() --------------------------------------------------
  //
  //            [ 0.3  -1.0   0.2 ]
  // Reuse A2 = [ 0.4   0.3   0.0 ] and lambda = 2.3
  //            [ 0.0  -1.0  -2.1 ]
  memcpy(A2.GetPointer(), mA2, N2*N2*sizeof(double));
  A2.SubtractDiagonal(lambda);
  double subtract_diag_truth[] = { -2.0, 0.4, 0.0, -1.0, -2.0, -1.0, 0.2, 0.0, -4.4 };
  derr = Diff1(A2.GetPointer(), subtract_diag_truth, N2*N2);
  printf("Test SubtractDiagonal() : Discrepancy in return matrix %g\n", derr);

  // Test SubtractDiagonal() --------------------------------------------------
  //
  // Reuse A2 and lambda
  memcpy(A2.GetPointer(), mA2, N2*N2*sizeof(double));
  A2.SubtractDiagonal(lambda, B2);
  derr = Diff1(B2.GetPointer(), subtract_diag_truth, N2*N2);
  printf("Test SubtractDiagonal() : Discrepancy in return matrix %g\n", derr);

  // Test Cos() ---------------------------------------------------------------
  //
  //           [  1.0  0.0  0.3  0.4 ]
  // Reuse A = [  0.0  0.0  1.0  0.2 ]
  //           [ -0.1  1.0  0.0  0.0 ]
  //
  // B =
  // [ 0.540302305868140 1.000000000000000 0.955336489125606 0.921060994002885 ]
  // [ 1.000000000000000 1.000000000000000 0.540302305868140 0.980066577841242 ]
  // [ 0.995004165278026 0.540302305868140 1.000000000000000 1.000000000000000 ]
  memcpy(A.GetPointer(), mA, M*N*sizeof(double));
  A.Cos();
  double cos_truth[] = { 0.540302305868140,
                         1.000000000000000,
                         0.995004165278026,
                         1.000000000000000,
                         1.000000000000000,
                         0.540302305868140,
                         0.955336489125606,
                         0.540302305868140,
                         1.000000000000000,
                         0.921060994002885,
                         0.980066577841242,
                         1.000000000000000 };
  derr = Diff1(A.GetPointer(), cos_truth, M*N);
  printf("Test Cos() : Discrepancy in return matrix %g\n", derr);

  // Test Cos() ---------------------------------------------------------------
  //
  // Reuse A
  memcpy(A.GetPointer(), mA, M*N*sizeof(double));
  A.Cos(B);
  derr = Diff1(B.GetPointer(), cos_truth, M*N);
  printf("Test Cos() : Discrepancy in return matrix %g\n", derr);

  // Test Sqrt() --------------------------------------------------------------
  //
  //     [ 1.0  0.0  0.3  0.4 ]
  // A = [ 0.0  0.0  1.0  0.2 ]
  //     [ 0.1  1.0  0.0  0.0 ]
  //
  // B =
  // [ 1.000000000000000 0.000000000000000 0.547722557505166 0.632455532033676 ]
  // [ 0.000000000000000 0.000000000000000 1.000000000000000 0.447213595499958 ]
  // [ 0.316227766016838 1.000000000000000 0.000000000000000 0.000000000000000 ]
  double mA11[] = { 1.0, 0.0, 0.1, 0.0, 0.0, 1.0, 0.3, 1.0, 0.0, 0.4, 0.2, 0.0};
  memcpy(A.GetPointer(), mA11, M*N*sizeof(double));
  A.Sqrt();
  double sqrt_truth[] = { 1.000000000000000,
                          0.000000000000000,
                          0.316227766016838,
                          0.000000000000000,
                          0.000000000000000,
                          1.000000000000000,
                          0.547722557505166,
                          1.000000000000000,
                          0.000000000000000,
                          0.632455532033676,
                          0.447213595499958,
                          0.000000000000000 };
  derr = Diff1(A.GetPointer(), sqrt_truth, M*N);
  printf("Test Sqrt() : Discrepancy in return matrix %g\n", derr);

  // Test Sqrt() --------------------------------------------------------------
  //
  // Reuse A
  memcpy(A.GetPointer(), mA11, M*N*sizeof(double));
  A.Sqrt(B);
  derr = Diff1(B.GetPointer(), sqrt_truth, M*N);
  printf("Test Sqrt() : Discrepancy in return matrix %g\n", derr);

  // Test Log() ---------------------------------------------------------------
  //
  //     [ 1.0  2.0  2.3  2.4 ]
  // A = [ 2.9  2.8  1.0  2.2 ]
  //     [ 2.1  1.0  2.0  2.0 ]
  //
  // B =
  // [ 0.000000000000000 0.693147180559945 0.832909122935104 0.875468737353900 ]
  // [ 1.064710736992428 1.029619417181158 0.000000000000000 0.788457360364270 ]
  // [ 0.741937344729377 0.000000000000000 0.693147180559945 0.693147180559945 ]
  double mA12[] = { 1.0, 2.9, 2.1, 2.0, 2.8, 1.0, 2.3, 1.0, 2.0, 2.4, 2.2, 2.0};
  memcpy(A.GetPointer(), mA12, M*N*sizeof(double));
  A.Log();
  double log_truth[] = { 0.000000000000000,
                         1.064710736992428,
                         0.741937344729377,
                         0.693147180559945,
                         1.029619417181158,
                         0.000000000000000,
                         0.832909122935104,
                         0.000000000000000,
                         0.693147180559945,
                         0.875468737353900,
                         0.788457360364270,
                         0.693147180559945 };
  derr = Diff1(A.GetPointer(), log_truth, M*N);
  printf("Test Log() : Discrepancy in return matrix %g\n", derr);

  // Test Log() ---------------------------------------------------------------
  //
  // Reuse A
  memcpy(A.GetPointer(), mA12, M*N*sizeof(double));
  A.Log(B);
  derr = Diff1(B.GetPointer(), log_truth, M*N);
  printf("Test Log() : Discrepancy in return matrix %g\n", derr);

  // Test Sum() ---------------------------------------------------------------
  //
  //           [ 1.0  2.0  2.3  2.4 ]
  // Reuse A = [ 2.9  2.8  1.0  2.2 ]
  //           [ 2.1  1.0  2.0  2.0 ]
  //
  // Sum(A,1) = [ 6.0  5.8  5.3  6.6 ]
  // Sum(A,2) = [ 7.7; 8.9; 7.1 ]
  memcpy(A.GetPointer(), mA12, M*N*sizeof(double));
  A.Sum(b,1);
  double sum_truth1[] = { 6.0, 5.8, 5.3, 6.6 };
  derr = Diff1(b.GetPointer(), sum_truth1, N);
  printf("Test Sum() : Discrepancy in return vector %g\n", derr);

  A.Sum(b,2);
  double sum_truth2[] = { 7.7, 8.9, 7.1 };
  derr = Diff1(b.GetPointer(), sum_truth2, M);
  printf("Test Sum() : Discrepancy in return vector %g\n", derr);

  // Test Prod() --------------------------------------------------------------
  //
  //           [ 1.0  2.0  2.3  2.4 ]
  // Reuse A = [ 2.9  2.8  1.0  2.2 ]
  //           [ 2.1  1.0  2.0  2.0 ]
  //
  // Prod(A,1) = [ 6.09  5.60  4.60  10.56 ]
  // Prod(A,2) = [ 11.040; 17.864; 8.400 ]
  memcpy(A.GetPointer(), mA12, M*N*sizeof(double));
  A.Prod(b,1);
  double prod_truth1[] = { 6.09, 5.60, 4.60, 10.56 };
  derr = Diff1(b.GetPointer(), prod_truth1, N);
  printf("Test Prod() : Discrepancy in return vector %g\n", derr);

  A.Prod(b,2);
  double prod_truth2[] = { 11.040, 17.864, 8.400 };
  derr = Diff1(b.GetPointer(), prod_truth2, M);
  printf("Test Prod() : Discrepancy in return vector %g\n", derr);

  // Test OuterProduct() ------------------------------------------------------
  //
  // b = [ 1.0; 0.0; -0.1 ]
  // c = [  1.0  0.0  0.3  0.4 ]
  //
  //     [  1.0  0.0  0.3   0.4  ]
  // A = [  0.0  0.0  0.0   0.0  ]
  //     [ -0.1  0.0 -0.03 -0.04 ]
  double mbb[] = { 1.0, 0.0, -0.1 };
  double mcc[] = { 1.0, 0.0, 0.3, 0.4 };
  DVector bb(3), cc(4);
  memcpy(bb.GetPointer(), mbb, 3*sizeof(double));
  memcpy(cc.GetPointer(), mcc, 4*sizeof(double));
  A.OuterProduct(bb, cc);
  double outprod[] = { 1.0, 0.0, -0.1, 0.0, 0.0, 0.0, 0.3, 0.0, -0.03, 0.4, 0.0, -0.04 };
  derr = Diff1(A.GetPointer(), outprod, 3*4);
  printf("Test OuterProduct() : Discrepancy in return matrix %g\n", derr);

  //---------------------------------------------------------------------------
  // Beginning from here, major matrices and vectors have a suffix
  // starting from 3
  //---------------------------------------------------------------------------

  // Test MatVec() ------------------------------------------------------------
  //
  // [  0.29 ]   [  1.0  0.0  0.3  0.4 ]   [  0.3 ]
  // [  0.08 ] = [  0.0  0.0  1.0  0.2 ] * [  0.0 ]
  // [ -0.03 ]   [ -0.1  1.0  0.0  0.0 ]   [  0.1 ]
  //                                       [ -0.1 ]
  //                        A3
  double mA3[] = { 1.0, 0.0, -0.1, 0.0, 0.0, 1.0, 0.3, 1.0, 0.0, 0.4, 0.2, 0.0 };
  INTEGER MA = 3, NA = 4;
  DMatrix A3(MA, NA);
  memcpy(A3.GetPointer(), mA3, MA*NA*sizeof(double));
  double mb3[] = { 0.3, 0.0, 0.1, -0.1 };
  DVector b3(NA);
  memcpy(b3.GetPointer(), mb3, NA*sizeof(double));
  DVector y3;
  A3.MatVec(b3, y3, NORMAL);
  double matvec_truth[] = { 0.29, 0.08, -0.03 };
  derr = Diff1(y3.GetPointer(), matvec_truth, MA);
  printf("Test MatVec() : Discrepancy in return vector %g\n", derr);

  double mA3T[] = { 1.0, 0.0, 0.3, 0.4, 0.0, 0.0, 1.0, 0.2, -0.1, 1.0, 0.0, 0.0 };
  DMatrix A3T(NA, MA);
  memcpy(A3T.GetPointer(), mA3T, MA*NA*sizeof(double));
  A3T.MatVec(b3, y3, TRANSPOSE);
  derr = Diff1(y3.GetPointer(), matvec_truth, MA);
  printf("Test MatVec() : Discrepancy in return vector %g\n", derr);

  // Test DGEMV() -------------------------------------------------------------
  //
  // [  2.94 ]          [  1.0  0.0  0.3  0.4 ]   [  0.3 ]         [  0.2 ]
  // [  0.2  ] = 10.0 * [  0.0  0.0  1.0  0.2 ] * [  0.0 ] + 0.2 * [ -3.0 ]
  // [ -0.32 ]          [ -0.1  1.0  0.0  0.0 ]   [  0.1 ]         [ -0.1 ]
  //                                              [ -0.1 ]
  double y3_init[] = { 0.2, -3.0, -0.1 };
  memcpy(y3.GetPointer(), y3_init, MA*sizeof(double));
  double alpha = 10.0, beta = 0.2;
  A3.DGEMV(b3, y3, alpha, beta, NORMAL);
  double dgemv_truth[] = { 2.94, 0.2, -0.32 };
  derr = Diff1(y3.GetPointer(), dgemv_truth, MA);
  printf("Test DGEMV() : Discrepancy in return vector %g\n", derr);

  memcpy(y3.GetPointer(), y3_init, MA*sizeof(double));
  A3T.DGEMV(b3, y3, alpha, beta, TRANSPOSE);
  derr = Diff1(y3.GetPointer(), dgemv_truth, MA);
  printf("Test DGEMV() : Discrepancy in return vector %g\n", derr);

  // Test MatMat() ------------------------------------------------------------
  //
  // [  0.03   0.60 ]   [  1.0  0.0  0.3  0.4 ]   [  0.0   1.0 ]
  // [  0.10  -0.20 ] = [  0.0  0.0  1.0  0.2 ] * [ -0.4   0.3 ]
  // [ -0.40   0.20 ]   [ -0.1  1.0  0.0  0.0 ]   [  0.1   0.0 ]
  //                                              [  0.0  -1.0 ]
  //                               A3
  INTEGER MB = 4, NB = 2;
  double mB3[] = { 0.0, -0.4, 0.1, 0.0, 1.0, 0.3, 0.0, -1.0 };
  DMatrix B3(MB, NB);
  memcpy(B3.GetPointer(), mB3, MB*NB*sizeof(double));
  DMatrix C3;
  A3.MatMat(B3, C3, NORMAL, NORMAL);
  double matmat_truth[] = { 0.03, 0.10, -0.40, 0.60, -0.20, 0.20 };
  derr = Diff1(C3.GetPointer(), matmat_truth, MA*NB);
  printf("Test MatMat() : Discrepancy in return matrix %g\n", derr);

  double mB3T[] = { 0.0, 1.0, -0.4, 0.3, 0.1, 0.0, 0.0, -1.0 };
  DMatrix B3T(NB, MB);
  memcpy(B3T.GetPointer(), mB3T, MB*NB*sizeof(double));
  A3.MatMat(B3T, C3, NORMAL, TRANSPOSE);
  derr = Diff1(C3.GetPointer(), matmat_truth, MA*NB);
  printf("Test MatMat() : Discrepancy in return matrix %g\n", derr);

  A3T.MatMat(B3, C3, TRANSPOSE, NORMAL);
  derr = Diff1(C3.GetPointer(), matmat_truth, MA*NB);
  printf("Test MatMat() : Discrepancy in return matrix %g\n", derr);

  A3T.MatMat(B3T, C3, TRANSPOSE, TRANSPOSE);
  derr = Diff1(C3.GetPointer(), matmat_truth, MA*NB);
  printf("Test MatMat() : Discrepancy in return matrix %g\n", derr);

  // Test DGEMM() -------------------------------------------------------------
  //
  // [  0.38   6.2  ]          [  1.0  0.0  0.3  0.4 ]   [  0.0   1.0 ]
  // [  1.46  -2.02 ] = 10.0 * [  0.0  0.0  1.0  0.2 ] * [ -0.4   0.3 ]
  // [ -4.6    2.24 ]          [ -0.1  1.0  0.0  0.0 ]   [  0.1   0.0 ]
  //                                                     [  0.0  -1.0 ]
  //                          [  0.4   1.0 ]
  //                  + 0.2 * [  2.3  -0.1 ]
  //                          [ -3.0   1.2 ]
  double C3_init[] = { 0.4, 2.3, -3.0, 1.0, -0.1, 1.2 };
  memcpy(C3.GetPointer(), C3_init, MA*NB*sizeof(double));
  A3.DGEMM(B3, C3, alpha, beta, NORMAL, NORMAL);
  double dgemm_truth[] = { 0.38, 1.46, -4.6, 6.2, -2.02, 2.24 };
  derr = Diff1(C3.GetPointer(), dgemm_truth, MA*NB);
  printf("Test DGEMM() : Discrepancy in return matrix %g\n", derr);

  memcpy(C3.GetPointer(), C3_init, MA*NB*sizeof(double));
  A3.DGEMM(B3T, C3, alpha, beta, NORMAL, TRANSPOSE);
  derr = Diff1(C3.GetPointer(), dgemm_truth, MA*NB);
  printf("Test DGEMM() : Discrepancy in return matrix %g\n", derr);

  memcpy(C3.GetPointer(), C3_init, MA*NB*sizeof(double));
  A3T.DGEMM(B3, C3, alpha, beta, TRANSPOSE, NORMAL);
  derr = Diff1(C3.GetPointer(), dgemm_truth, MA*NB);
  printf("Test DGEMM() : Discrepancy in return matrix %g\n", derr);

  memcpy(C3.GetPointer(), C3_init, MA*NB*sizeof(double));
  A3T.DGEMM(B3T, C3, alpha, beta, TRANSPOSE, TRANSPOSE);
  derr = Diff1(C3.GetPointer(), dgemm_truth, MA*NB);
  printf("Test DGEMM() : Discrepancy in return matrix %g\n", derr);

  // Test Mldivide() ----------------------------------------------------------
  //
  // General case, single right-hand side
  // [ -0.66 ]   [ 0.3  -1.0   0.2 ]   [ -0.4 ]
  // [  0.14 ] = [ 0.4   0.3   0.0 ] * [  1.0 ]
  // [ -5.83 ]   [ 0.0  -1.0  -2.1 ]   [  2.3 ]
  //                      A4
  NA = 3;
  double mA4[] = { 0.3, 0.4, 0.0, -1.0, 0.3, -1.0, 0.2, 0.0, -2.1 };
  DMatrix A4(NA);
  memcpy(A4.GetPointer(), mA4, NA*NA*sizeof(double));
  double mb4[] = { -0.66, 0.14, -5.83 };
  DVector b4(NA);
  memcpy(b4.GetPointer(), mb4, NA*sizeof(double));
  DVector x4;
  A4.Mldivide(b4, x4, NORMAL, GENERAL);
  DVector r4;
  A4.MatVec(x4, r4, NORMAL);
  r4.Subtract(b4);
  double relres = r4.Norm2() / b4.Norm2();
  printf("Test Mldivide() : Relative residual %g\n", relres);

  double mA4T[] = { 0.3, -1.0, 0.2, 0.4, 0.3, 0.0, 0.0, -1.0, -2.1 };
  DMatrix A4T(NA);
  memcpy(A4T.GetPointer(), mA4T, NA*NA*sizeof(double));
  A4T.Mldivide(b4, x4, TRANSPOSE, GENERAL);
  A4.MatVec(x4, r4, NORMAL);
  r4.Subtract(b4);
  relres = r4.Norm2() / b4.Norm2();
  printf("Test Mldivide() : Relative residual %g\n", relres);

  // Test Mldivide() ----------------------------------------------------------
  //
  // General case, multiple right-hand sides
  // [ -0.66  -1.6  ]   [ 0.3  -1.0   0.2 ]   [ -0.4   1.4 ]
  // [  0.14   1.16 ] = [ 0.4   0.3   0.0 ] * [  1.0   2.0 ]
  // [ -5.83  -1.79 ]   [ 0.0  -1.0  -2.1 ]   [  2.3  -0.1 ]
  INTEGER NRHS = 2;
  double mB4[] = { -0.66, 0.14, -5.83, -1.6, 1.16, -1.79 };
  DMatrix B4(NA, NRHS);
  memcpy(B4.GetPointer(), mB4, NA*NRHS*sizeof(double));
  DMatrix X4;
  A4.Mldivide(B4, X4, NORMAL, GENERAL);
  DMatrix R4;
  A4.MatMat(X4, R4, NORMAL, NORMAL);
  R4.Subtract(B4);
  R4.GetColumn(0, r4); B4.GetColumn(0, b4);
  double relres0 = r4.Norm2() / b4.Norm2();
  R4.GetColumn(1, r4); B4.GetColumn(1, b4);
  double relres1 = r4.Norm2() / b4.Norm2();
  printf("Test Mldivide() : Relative residuals %g %g\n", relres0, relres1);

  A4T.Mldivide(B4, X4, TRANSPOSE, GENERAL);
  A4.MatMat(X4, R4, NORMAL, NORMAL);
  R4.Subtract(B4);
  R4.GetColumn(0, r4); B4.GetColumn(0, b4);
  relres0 = r4.Norm2() / b4.Norm2();
  R4.GetColumn(1, r4); B4.GetColumn(1, b4);
  relres1 = r4.Norm2() / b4.Norm2();
  printf("Test Mldivide() : Relative residuals %g %g\n", relres0, relres1);

  // Test Mldivide() ----------------------------------------------------------
  //
  // SPD case, single right-hand side
  // [ -1.06 ]   [  1.3  -1.0  0.2 ]   [ -0.4 ]
  // [  1.7  ] = [ -1.0   1.3  0.0 ] * [  1.0 ]
  // [  4.75 ]   [  0.2   0.0  2.1 ]   [  2.3 ]
  //                       A5
  double mA5[] = { 1.3, -1.0, 0.2, -1.0, 1.3, 0.0, 0.2, 0.0, 2.1 };
  DMatrix A5(NA);
  memcpy(A5.GetPointer(), mA5, NA*NA*sizeof(double));
  double mb5[] = { -1.06, 1.7, 4.75 };
  DVector b5(NA);
  memcpy(b5.GetPointer(), mb5, NA*sizeof(double));
  DVector x5;
  A5.Mldivide(b5, x5, NORMAL, SPD);
  DVector r5;
  A5.MatVec(x5, r5, NORMAL);
  r5.Subtract(b5);
  relres = r5.Norm2() / b5.Norm2();
  printf("Test Mldivide() : Relative residual %g\n", relres);

  // Test Mldivide() ----------------------------------------------------------
  //
  // SPD case, multiple right-hand sides
  // [ -1.06  -0.2  ]   [  1.3  -1.0  0.2 ]   [ -0.4   1.4 ]
  // [  1.7    1.2  ] = [ -1.0   1.3  0.0 ] * [  1.0   2.0 ]
  // [  4.75   0.07 ]   [  0.2   0.0  2.1 ]   [  2.3  -0.1 ]
  double mB5[] = { -1.06, 1.7, 4.75, -0.2, 1.2, 0.07 };
  DMatrix B5(NA, NRHS);
  memcpy(B5.GetPointer(), mB5, NA*NRHS*sizeof(double));
  DMatrix X5;
  A5.Mldivide(B5, X5, NORMAL, SPD);
  DMatrix R5;
  A5.MatMat(X5, R5, NORMAL, NORMAL);
  R5.Subtract(B5);
  R5.GetColumn(0, r5); B5.GetColumn(0, b5);
  relres0 = r5.Norm2() / b5.Norm2();
  R5.GetColumn(1, r5); B5.GetColumn(1, b5);
  relres1 = r5.Norm2() / b5.Norm2();
  printf("Test Mldivide() : Relative residuals %g %g\n", relres0, relres1);

  // Test Mldivide() ----------------------------------------------------------
  //
  // Least squares case, single right-hand side
  // [ 0.3  -1.0   0.2 ]              [ -0.4 ]
  // [ 0.4   0.3   0.0 ] * x \approx= [  1.0 ]
  // [ 0.0  -1.0  -2.1 ]              [  2.3 ]
  // [ 1.0  -0.4   1.7 ]              [ -2.5 ]
  //         A10
  double mA10[] = { 0.3, 0.4, 0.0, 1.0, -1.0, 0.3, -1.0, -0.4, 0.2, 0.0, -2.1, 1.7 };
  MA = 4;
  DMatrix A10(MA, NA);
  memcpy(A10.GetPointer(), mA10, MA*NA*sizeof(double));
  double mb10[] = { -0.4, 1.0, 2.3, -2.5 };
  DVector b10(MA);
  memcpy(b10.GetPointer(), mb10, MA*sizeof(double));
  DVector x10;
  double Res[2];
  A10.Mldivide(b10, x10, Res);
  DVector r10;
  A10.MatVec(x10, r10, NORMAL);
  r10.Subtract(b10);
  derr = fabs(r10.Norm2() - Res[0]) /  Res[0];
  printf("Test Mldivide() : Discrepancy in relative residual %g\n", derr);

  // Test Mldivide() ----------------------------------------------------------
  //
  // Least squares case, multiple right-hand sides
  // [ 0.3  -1.0   0.2 ]              [ -0.4   1.4 ]
  // [ 0.4   0.3   0.0 ] * x \approx= [  1.0   2.0 ]
  // [ 0.0  -1.0  -2.1 ]              [  2.3  -0.1 ]
  // [ 1.0  -0.4   1.7 ]              [ -2.5   0.3 ]
  double mB10[] = { -0.4, 1.0, 2.3, -2.5, 1.4, 2.0, -0.1, 0.3 };
  DMatrix B10(MA, NRHS);
  memcpy(B10.GetPointer(), mB10, MA*NRHS*sizeof(double));
  DMatrix X10;
  A10.Mldivide(B10, X10, Res);
  DMatrix R10;
  A10.MatMat(X10, R10, NORMAL, NORMAL);
  R10.Subtract(B10);
  R10.GetColumn(0, r10);
  double derr0 = fabs(r10.Norm2() - Res[0]) /  Res[0];
  R10.GetColumn(1, r10);
  double derr1 = fabs(r10.Norm2() - Res[1]) /  Res[1];
  printf("Test Mldivide() : Discrepancy in relative residuals %g %g\n", derr0, derr1);

  // Test DGETRS() ------------------------------------------------------------
  //
  // General case, reuse A4 as A6
  DMatrix A6(NA);
  memcpy(A6.GetPointer(), mA4, NA*NA*sizeof(double));
  A6.DGETRF();
  A6.DGETRS(b4, x4, NORMAL);
  A4.MatVec(x4, r4, NORMAL);
  r4.Subtract(b4);
  relres = r4.Norm2() / b4.Norm2();
  printf("Test DGETRS() : Relative residual %g\n", relres);

  A6.DGETRS(B4, X4, NORMAL);
  A4.MatMat(X4, R4, NORMAL, NORMAL);
  R4.Subtract(B4);
  R4.GetColumn(0, r4); B4.GetColumn(0, b4);
  relres0 = r4.Norm2() / b4.Norm2();
  R4.GetColumn(1, r4); B4.GetColumn(1, b4);
  relres1 = r4.Norm2() / b4.Norm2();
  printf("Test DGETRS() : Relative residuals %g %g\n", relres0, relres1);

  memcpy(A6.GetPointer(), mA4, NA*NA*sizeof(double));
  A6.DGETRF();
  A6.DGETRS(b4, x4, TRANSPOSE);
  A4.MatVec(x4, r4, TRANSPOSE);
  r4.Subtract(b4);
  relres = r4.Norm2() / b4.Norm2();
  printf("Test DGETRS() : Relative residual %g\n", relres);

  A6.DGETRS(B4, X4, TRANSPOSE);
  A4.MatMat(X4, R4, TRANSPOSE, NORMAL);
  R4.Subtract(B4);
  R4.GetColumn(0, r4); B4.GetColumn(0, b4);
  relres0 = r4.Norm2() / b4.Norm2();
  R4.GetColumn(1, r4); B4.GetColumn(1, b4);
  relres1 = r4.Norm2() / b4.Norm2();
  printf("Test DGETRS() : Relative residuals %g %g\n", relres0, relres1);

  // Test DPOTRS() ------------------------------------------------------------
  //
  // SPD case, reuse A5 as A6
  memcpy(A6.GetPointer(), mA5, NA*NA*sizeof(double));
  A6.DPOTRF(LOWER);
  A6.DPOTRS(b4, x4, LOWER);
  A5.MatVec(x4, r4, NORMAL);
  r4.Subtract(b4);
  relres = r4.Norm2() / b4.Norm2();
  printf("Test DPOTRS() : Relative residual %g\n", relres);

  A6.DPOTRS(B4, X4, LOWER);
  A5.MatMat(X4, R4, NORMAL, NORMAL);
  R4.Subtract(B4);
  R4.GetColumn(0, r4); B4.GetColumn(0, b4);
  relres0 = r4.Norm2() / b4.Norm2();
  R4.GetColumn(1, r4); B4.GetColumn(1, b4);
  relres1 = r4.Norm2() / b4.Norm2();
  printf("Test DPOTRS() : Relative residuals %g %g\n", relres0, relres1);

  // Test DTRSV() -------------------------------------------------------------
  //
  // SPD case, reuse A6 above, and set A7 to be the lower tri part of A6
  DMatrix A7 = A6;
  A7.SetEntry(0, 1, 0.0); A7.SetEntry(0, 2, 0.0); A7.SetEntry(1, 2, 0.0);
  A6.DTRSV(b4, x4, NORMAL, LOWER);
  A7.MatVec(x4, r4, NORMAL);
  r4.Subtract(b4);
  relres = r4.Norm2() / b4.Norm2();
  printf("Test DTRSV() : Relative residual %g\n", relres);

  // Test DTRSM() -------------------------------------------------------------
  //
  // SPD case, reuse A6 and A7 above
  A6.DTRSM(B4, X4, NORMAL, LOWER);
  A7.MatMat(X4, R4, NORMAL, NORMAL);
  R4.Subtract(B4);
  R4.GetColumn(0, r4); B4.GetColumn(0, b4);
  relres0 = r4.Norm2() / b4.Norm2();
  R4.GetColumn(1, r4); B4.GetColumn(1, b4);
  relres1 = r4.Norm2() / b4.Norm2();
  printf("Test DTRSM() : Relative residuals %g %g\n", relres0, relres1);

  // Test Chol() --------------------------------------------------------------
  //
  // Reuse A5
  DMatrix G;
  A5.Chol(G, LOWER);
  G.MatMat(G, A6, NORMAL, TRANSPOSE);
  double res = Diff1(A5.GetPointer(), A6.GetPointer(), NA*NA);
  printf("Test Chol() : Residual %g\n", res);

  // Test SymEig() ------------------------------------------------------------
  //
  //      [  0.3  -1.0  0.2 ]
  // A8 = [ -1.0   1.3  0.0 ]
  //      [  0.2   0.0  2.1 ]
  double mA8[] = { 0.3, -1.0, 0.2, -1.0, 1.3, 0.0, 0.2, 0.0, 2.1 };
  DMatrix A8(NA);
  memcpy(A8.GetPointer(), mA8, NA*NA*sizeof(double));
  DVector eigval;
  DMatrix V;
  A8.SymEig(eigval, V);
  DVector v, Av, r;
  double s;
  V.GetColumn(0, v); s = eigval.GetEntry(0);
  A8.MatVec(v, Av, NORMAL); v.Multiply(s); v.Subtract(Av, r);
  double res0 = r.Norm2();
  V.GetColumn(1, v); s = eigval.GetEntry(1);
  A8.MatVec(v, Av, NORMAL); v.Multiply(s); v.Subtract(Av, r);
  double res1 = r.Norm2();
  V.GetColumn(2, v); s = eigval.GetEntry(2);
  A8.MatVec(v, Av, NORMAL); v.Multiply(s); v.Subtract(Av, r);
  double res2 = r.Norm2();
  printf("Test SymEig() : Residuals %g %g %g\n", res0, res1, res2);

  // Test SymEigBIndef() ------------------------------------------------------
  //
  //            [  0.3  -1.0  0.2 ]        [  1.3  -1.0  0.2 ]
  // Reuse A8 = [ -1.0   1.3  0.0 ],  A9 = [ -1.0   2.3  0.0 ]
  //            [  0.2   0.0  2.1 ]        [  0.2   0.0  2.1 ]
  double mA9[] = { 1.3, -1.0, 0.2, -1.0, 2.3, 0.0, 0.2, 0.0, 2.1 };
  DMatrix A9(NA);
  memcpy(A9.GetPointer(), mA9, NA*NA*sizeof(double));
  A8.SymEigBIndef(A9, eigval, V);
  DVector Bv;
  V.GetColumn(0, v); s = eigval.GetEntry(0);
  A8.MatVec(v, Av, NORMAL); A9.MatVec(v, Bv, NORMAL);
  Bv.Multiply(s); Bv.Subtract(Av, r);
  res0 = r.Norm2();
  V.GetColumn(1, v); s = eigval.GetEntry(1);
  A8.MatVec(v, Av, NORMAL); A9.MatVec(v, Bv, NORMAL);
  Bv.Multiply(s); Bv.Subtract(Av, r);
  res1 = r.Norm2();
  V.GetColumn(2, v); s = eigval.GetEntry(2);
  A8.MatVec(v, Av, NORMAL); A9.MatVec(v, Bv, NORMAL);
  Bv.Multiply(s); Bv.Subtract(Av, r);
  res2 = r.Norm2();
  printf("Test SymEigBIndef() : Residuals %g %g %g\n", res0, res1, res2);

  // Test RealSchur() ---------------------------------------------------------
  //
  // Reuse A4
  DMatrix T, U, UT, UTU;
  A4.RealSchur(T, U, LHP);
  U.MatMat(T, UT, NORMAL, NORMAL);
  UT.MatMat(U, UTU, NORMAL, TRANSPOSE);
  res = Diff1(A4.GetPointer(), UTU.GetPointer(), NA*NA);
  printf("Test RealSchur() : Residual %g\n", res);

  // Test Orth() --------------------------------------------------------------
  //
  // Reuse A4
  DMatrix Q, QQ;
  A4.Orth(Q);
  Q.MatMat(Q, QQ, TRANSPOSE, NORMAL);
  double QQ_truth[] = { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };
  derr = Diff1(QQ.GetPointer(), QQ_truth, NA*NA);
  printf("Test Orth() : Residual %g\n", derr);

  // Test Rank() --------------------------------------------------------------
  //
  // Reuse A10. Rank = 3
  INTEGER rank = A10.Rank();
  int rank_truth = 3;
  int ierr = abs((int)rank - rank_truth);
  printf("Test Rank() : Discrepancy in return value %ld\n", (long)ierr);

  // Test Cond2() -------------------------------------------------------------
  //
  // Reuse A10. Cond2 = 4.828956246106313
  double cond2 = A10.Cond2();
  double cond2_truth = 4.828956246106313;
  derr = fabs(cond2 - cond2_truth);
  printf("Test Cond2() : Discrepancy in return value %g\n", derr);

  // Test Norm2() -------------------------------------------------------------
  //
  // Reuse A10. Norm2 = 2.828313508109951
  double norm2 = A10.Norm2();
  double norm2_truth = 2.828313508109951;
  derr = fabs(norm2 - norm2_truth);
  printf("Test Norm2() : Discrepancy in return value %g\n", derr);

  // Test NormF() -------------------------------------------------------------
  //
  // Reuse A10. NormF = 3.292415526630866
  double normf = A10.NormF();
  double normf_truth = 3.292415526630866;
  derr = fabs(normf - normf_truth);
  printf("Test NormF() : Discrepancy in return value %g\n", derr);

  // Test Diag() --------------------------------------------------------------
  //
  // Reuse A10. Diag = [ 0.3; 0.3; -2.1 ]
  DVector diag;
  A10.Diag(diag);
  double diag_truth2[] = { 0.3, 0.3, -2.1 };
  derr = Diff1(diag.GetPointer(), diag_truth2, 3);
  printf("Test Diag() : Discrepancy in return value %g\n", derr);

  // Test Tr() ----------------------------------------------------------------
  //
  // Reuse A10. Tr = -1.5
  double tr = A10.Tr();
  double tr_truth = -1.5;
  derr = fabs(tr - tr_truth);
  printf("Test Tr() : Discrepancy in return value %g\n", derr);

  // Test Det() ---------------------------------------------------------------
  //
  // Reuse A4. Might be factorized into A6. Det = -1.109
  LogDet LD = A4.Det(GENERAL, UNFACT);
  double det = exp(LD.LogAbsDet)*LD.Sign;
  double det_truth = -1.109;
  derr = fabs(det - det_truth);
  printf("Test Det() : Discrepancy in return value %g\n", derr);

  A6 = A4;
  A6.DGETRF();
  LD = A6.Det(GENERAL, LU_FACT);
  det = exp(LD.LogAbsDet)*LD.Sign;
  det_truth = -1.109;
  derr = fabs(det - det_truth);
  printf("Test Det() : Discrepancy in return value %g\n", derr);

  // Test Det() ---------------------------------------------------------------
  //
  // Reuse A5. Might be factorized into A6. Det = 1.397
  LD = A5.Det(SPD, UNFACT);
  det = exp(LD.LogAbsDet)*LD.Sign;
  det_truth = 1.397;
  derr = fabs(det - det_truth);
  printf("Test Det() : Discrepancy in return value %g\n", derr);

  A6 = A5;
  A6.DPOTRF(LOWER);
  LD = A6.Det(SPD, CHOL_FACT);
  det = exp(LD.LogAbsDet)*LD.Sign;
  det_truth = 1.397;
  derr = fabs(det - det_truth);
  printf("Test Det() : Discrepancy in return value %g\n", derr);

  // Test BuildKernelMatrix() -------------------------------------------------
  //
  //     [ -1.06  -0.2  ]        [ -0.4   1.4 ]
  // X = [  1.7    1.2  ],   Y = [  1.0   2.0 ]
  //     [  4.75   0.07 ]        [  2.3  -0.1 ]
  //                             [ -2.5   0.3 ]
  //
  // Each row is a point; but the matrix is stored in column major order.
  // Kernel: Gaussian. s = 1.2. sigma = 1.4. lambda = 0.4.
  //
  // K(X,X) =
  // [ 1.600000000000000 0.104252050510782 0.000214430392665 ]
  // [ 0.104252050510782 1.600000000000000 0.080741307170112 ]
  // [ 0.000214430392665 0.080741307170112 1.600000000000000 ]
  //
  // K(X,Y) =
  // [ 0.558856740596001 0.118253655216597 0.067190093289071 0.663362590086124 ]
  // [ 0.385627838011378 0.899475242739326 0.711316623533559 0.010842181083260 ]
  // [ 0.000880630595997 0.012838837990185 0.257611951002296 0.000001778083673 ]
  INTEGER NX = 3;
  INTEGER DIM = 2;
  double mPX[] = { -1.06, 1.7, 4.75, -0.2, 1.2, 0.07 };
  DPointArray PX(NX, DIM);
  memcpy(PX.GetPointer(), mPX, NX*DIM*sizeof(double));

  INTEGER NY = 4;
  double mPY[] = { -0.4, 1.0, 2.3, -2.5, 1.4, 2.0, -0.1, 0.3 };
  DPointArray PY(NY, DIM);
  memcpy(PY.GetPointer(), mPY, NY*DIM*sizeof(double));

  s = 1.2;
  double sigma = 1.4;
  lambda = 0.4;
  IsotropicGaussian mKernel(s, sigma);

  DMatrix KXX;
  KXX.BuildKernelMatrix<IsotropicGaussian, DPoint, DPointArray>
    (mKernel, PX, lambda);
  double KXX_truth[] = { 1.600000000000000,
                         0.104252050510782,
                         0.000214430392665,
                         0.104252050510782,
                         1.600000000000000,
                         0.080741307170112,
                         0.000214430392665,
                         0.080741307170112,
                         1.600000000000000 };
  derr = Diff1(KXX.GetPointer(), KXX_truth, NX*NX);
  printf("Test BuildKernelMatrix() : Discrepancy in return matrix %g\n", derr);

  DMatrix KXY;
  KXY.BuildKernelMatrix<IsotropicGaussian, DPoint, DPointArray>
    (mKernel, PX, PY, lambda);
  double KXY_truth[] = { 0.558856740596001,
                         0.385627838011378,
                         0.000880630595997,
                         0.118253655216597,
                         0.899475242739326,
                         0.012838837990185,
                         0.067190093289071,
                         0.711316623533559,
                         0.257611951002296,
                         0.663362590086124,
                         0.010842181083260,
                         0.000001778083673 };
  derr = Diff1(KXY.GetPointer(), KXY_truth, NX*NY);
  printf("Test BuildKernelMatrix() : Discrepancy in return matrix %g\n", derr);

  return 0;

}
