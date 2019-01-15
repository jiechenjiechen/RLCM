// Unit test of the class Lanczos.

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
  printf("|  Test_Lanczos.cpp                                             |\n");
  printf("+---------------------------------------------------------------+\n");

  // Test Run() : Full reorth, no early stopping ------------------------------
  //
  // A is the Standard Laplacian tridiag(-1, 2, -1)
  // v = ones(n,1)
  INTEGER N = 1000;
  INTEGER nnz = 3*N-2;
  
  SMatrix A(N, nnz);
  INTEGER *A_start = A.GetPointerStart();
  A_start[0] = 0;
  A_start[1] = 2;
  for (INTEGER i = 2; i < N; i++) {
    A_start[i] = A_start[i-1] + 3;
  }
  A_start[N] = A_start[N-1] + 2;

  INTEGER *A_idx = A.GetPointerIdx();
  double *A_val = A.GetPointerA();
  INTEGER j = 0;
  for (INTEGER i = 0; i < N; i++) {
    if (i != 0) {
      A_idx[j] = i-1;
      A_val[j++] = -1.0;
    }
    A_idx[j] = i;
    A_val[j++] = 2.0;
    if (i != N-1) {
      A_idx[j] = i+1;
      A_val[j++] = -1.0;
    }
  }

  DVector v(N);
  v.SetConstVal(1.0);

  bool PartialReorth = false;
  INTEGER MaxIt = 20;
  bool EarlyStop = false;
  INTEGER k = 0;
  double RTol = 0;

  Lanczos mLanczos;
  mLanczos.Run<SMatrix>(A, UNFACT, v, PartialReorth, MaxIt, EarlyStop, k, RTol);

  double NormV = mLanczos.GetNormV();
  double NormV_truth = sqrt((double)N);
  double err_NormV = fabs(NormV - NormV_truth);

  INTEGER Iter = mLanczos.GetIter();
  INTEGER Iter_truth = MaxIt;
  double err_Iter = fabs(Iter - Iter_truth);

  DMatrix T, V;
  mLanczos.GetT(T);
  mLanczos.GetV(V);

  DMatrix AV, VT, E, E2;
  A.MatMat(V, AV, NORMAL, NORMAL);
  V.MatMat(T, VT, NORMAL, NORMAL);
  AV.Subtract(VT, E);
  E.GetBlock(0, N, 0, MaxIt-1, E2);
  double err_AV_VT = E2.NormF() / (double)(MaxIt-1);

  DMatrix VV, Eye(MaxIt);
  V.MatMat(V, VV, TRANSPOSE, NORMAL);
  Eye.SetIdentity();
  VV.Subtract(Eye, E);
  double err_VV_Eye = E.NormF() / (double)MaxIt;

  printf("Test Run() : Full reorth, no early stopping\n");
  printf("Error in calculating norm(v): %g\n", err_NormV);
  printf("Error in the actual number of iterations: %g\n", err_Iter);
  printf("Difference between AV and VT: %g\n", err_AV_VT);
  printf("Difference between VV and Eye: %g\n", err_VV_Eye);

  // Test Run() : Full reorth, early stopping ---------------------------------

  EarlyStop = true;
  k = 5;
  RTol = 1e-8;

  mLanczos.Run<SMatrix>(A, UNFACT, v, PartialReorth, MaxIt, EarlyStop, k, RTol);

  Iter = mLanczos.GetIter();
  const double *ResHistory = mLanczos.GetRes();

  DVector S;
  DMatrix U;
  mLanczos.GetRitzValues(S);
  mLanczos.GetRitzVectors(U);

  double *err_Au_su = NULL;
  New_1D_Array<double, INTEGER>(&err_Au_su, MaxIt);
  for (INTEGER i = 0; i < Iter; i++) {
    double s;
    DVector u, Au, su;
    U.GetColumn(i, u);
    s = S.GetEntry(i);
    A.MatVec(u, Au, NORMAL);
    u.Multiply(s, su);
    su.Subtract(Au);
    err_Au_su[i] = su.Norm2();
  }

  printf("Test Run() : Full reorth, early stopping\n");
  printf("Difference between Au and su: truth vs computed\n");
  for (INTEGER i = 0; i < Iter; i++) {
    printf("[%ld]  %.7e  %.7e\n", (long)i, err_Au_su[i], ResHistory[i]);
  }

  // Test Run() : Partial reorth, no early stopping ---------------------------

  PartialReorth = true;
  EarlyStop = false;

  mLanczos.Run<SMatrix>(A, UNFACT, v, PartialReorth, MaxIt, EarlyStop, k, RTol);

  mLanczos.GetT(T);
  mLanczos.GetV(V);

  A.MatMat(V, AV, NORMAL, NORMAL);
  V.MatMat(T, VT, NORMAL, NORMAL);
  AV.Subtract(VT, E);
  E.GetBlock(0, N, 0, MaxIt-1, E2);
  err_AV_VT = E2.NormF() / (double)(MaxIt-1);

  V.MatMat(V, VV, TRANSPOSE, NORMAL);
  VV.Subtract(Eye, E);
  err_VV_Eye = E.NormF() / (double)MaxIt;

  printf("Test Run() : Partial reorth, no early stopping\n");
  printf("Difference between AV and VT: %g\n", err_AV_VT);
  printf("Difference between VV and Eye: %g\n", err_VV_Eye);

  // Test Run() : Partial reorth, early stopping ------------------------------

  EarlyStop = true;
  k = 5;
  RTol = 1e-8;

  mLanczos.Run<SMatrix>(A, UNFACT, v, PartialReorth, MaxIt, EarlyStop, k, RTol);

  Iter = mLanczos.GetIter();
  ResHistory = mLanczos.GetRes();

  mLanczos.GetRitzValues(S);
  mLanczos.GetRitzVectors(U);

  for (INTEGER i = 0; i < Iter; i++) {
    double s;
    DVector u, Au, su;
    U.GetColumn(i, u);
    s = S.GetEntry(i);
    A.MatVec(u, Au, NORMAL);
    u.Multiply(s, su);
    su.Subtract(Au);
    err_Au_su[i] = su.Norm2();
  }

  printf("Test Run() : Partial reorth, early stopping\n");
  printf("Difference between Au and su: truth vs computed\n");
  for (INTEGER i = 0; i < Iter; i++) {
    printf("[%ld]  %.7e  %.7e\n", (long)i, err_Au_su[i], ResHistory[i]);
  }

  // Clean up -----------------------------------------------------------------

  Delete_1D_Array<double>(&err_Au_su);

  return 0;

}
