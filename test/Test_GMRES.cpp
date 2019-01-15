// Unit test of the class GMRES.

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
  printf("|  Test_GMRES.cpp                                               |\n");
  printf("+---------------------------------------------------------------+\n");

  INTEGER MaxIt = String2Integer(argv[2]);

  // Test Solve() : No preconditioner -----------------------------------------
  //
  // A is the Standard Laplacian tridiag(-1, 2, -1)
  // x = ones(n,1)
  INTEGER N = 100;
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

  DVector x_truth(N);
  x_truth.SetConstVal(1.0);

  DVector b;
  A.MatVec(x_truth, b, NORMAL);

  DVector x0(N);
  x0.SetConstVal(0.0);

  SMatrix Eye(N, N);
  INTEGER *Eye_start = Eye.GetPointerStart();
  INTEGER *Eye_idx = Eye.GetPointerIdx();
  double *Eye_val = Eye.GetPointerA();
  for (INTEGER i = 0; i < N; i++) {
    Eye_start[i] = i;
    Eye_idx[i] = i;
    Eye_val[i] = 1.0;
  }
  Eye_start[N] = N;

  INTEGER m = 100;
  double RTol = 1e-8;

  GMRES mGMRES;
  mGMRES.Solve<SMatrix, SMatrix>(A, UNFACT, b, x0, Eye, UNFACT, m, MaxIt, RTol);

  double NormB = mGMRES.GetNormRHS();

  DVector x;
  mGMRES.GetSolution(x);

  INTEGER Iter = -1;
  const double *ResHistory = mGMRES.GetResHistory(Iter);

  DVector r;
  A.MatVec(x, r, NORMAL);
  r.Subtract(b);
  double NormR = r.Norm2();

  printf("Test Solve() : No preconditioner\n");
  printf("Restart cycle: %ld\n", (long)m);
  printf("Number of iterations: %ld\n", (long)Iter);
  printf("Actural relative residual: %g\n", NormR/NormB);
  printf("Relative residual reported by the solver: %g\n", ResHistory[Iter-1]/NormB);

  const char *out_file1 = "Test_GMRES.relres1";
  FILE *fid = fopen(out_file1, "w");
  for (INTEGER i = 0; i < Iter; i++) {
    fprintf(fid, "%g\n", ResHistory[i]/NormB);
  }
  fclose(fid);
  printf("Relative residual history is written to file %s. Use the Matlab program Test_GMRES.m to plot it out!\n", out_file1);

  // Test Solve() : With preconditioner ---------------------------------------
  //
  // M = tridiag(2/9, 6/9, 2/9)
  SMatrix M(N, nnz);
  INTEGER *M_start = M.GetPointerStart();
  M_start[0] = 0;
  M_start[1] = 2;
  for (INTEGER i = 2; i < N; i++) {
    M_start[i] = M_start[i-1] + 3;
  }
  M_start[N] = M_start[N-1] + 2;

  INTEGER *M_idx = M.GetPointerIdx();
  double *M_val = M.GetPointerA();
  j = 0;
  for (INTEGER i = 0; i < N; i++) {
    if (i != 0) {
      M_idx[j] = i-1;
      M_val[j++] = 2.0/9;
    }
    M_idx[j] = i;
    M_val[j++] = 6.0/9;
    if (i != N-1) {
      M_idx[j] = i+1;
      M_val[j++] = 2.0/9;
    }
  }

  mGMRES.Solve<SMatrix, SMatrix>(A, UNFACT, b, x0, M, UNFACT, m, MaxIt, RTol);

  mGMRES.GetSolution(x);

  ResHistory = mGMRES.GetResHistory(Iter);

  A.MatVec(x, r, NORMAL);
  r.Subtract(b);
  NormR = r.Norm2();

  printf("Test Solve() : With preconditioner\n");
  printf("Restart cycle: %ld\n", (long)m);
  printf("Number of iterations: %ld\n", (long)Iter);
  printf("Actural relative residual: %g\n", NormR/NormB);
  printf("Relative residual reported by the solver: %g\n", ResHistory[Iter-1]/NormB);

  const char *out_file2 = "Test_GMRES.relres2";
  fid = fopen(out_file2, "w");
  for (INTEGER i = 0; i < Iter; i++) {
    fprintf(fid, "%g\n", ResHistory[i]/NormB);
  }
  fclose(fid);
  printf("Relative residual history is written to file %s. Use the Matlab program Test_GMRES.m to plot it out!\n", out_file2);

  // Test Solve() : No preconditioner -----------------------------------------
  //
  // Small restart cycle m = 10
  m = 10;

  mGMRES.Solve<SMatrix, SMatrix>(A, UNFACT, b, x0, Eye, UNFACT, m, MaxIt, RTol);

  mGMRES.GetSolution(x);

  ResHistory = mGMRES.GetResHistory(Iter);

  A.MatVec(x, r, NORMAL);
  r.Subtract(b);
  NormR = r.Norm2();

  printf("Test Solve() : No preconditioner\n");
  printf("Restart cycle: %ld\n", (long)m);
  printf("Number of iterations: %ld\n", (long)Iter);
  printf("Actural relative residual: %g\n", NormR/NormB);
  printf("Relative residual reported by the solver: %g\n", ResHistory[Iter-1]/NormB);

  const char *out_file3 = "Test_GMRES.relres3";
  fid = fopen(out_file3, "w");
  for (INTEGER i = 0; i < Iter; i++) {
    fprintf(fid, "%g\n", ResHistory[i]/NormB);
  }
  fclose(fid);
  printf("Relative residual history is written to file %s. Use the Matlab program Test_GMRES.m to plot it out!\n", out_file3);

  // Test Solve() : With preconditioner ---------------------------------------
  //
  // Small restart cycle m = 10
  mGMRES.Solve<SMatrix, SMatrix>(A, UNFACT, b, x0, M, UNFACT, m, MaxIt, RTol);

  mGMRES.GetSolution(x);

  ResHistory = mGMRES.GetResHistory(Iter);

  A.MatVec(x, r, NORMAL);
  r.Subtract(b);
  NormR = r.Norm2();

  printf("Test Solve() : With preconditioner\n");
  printf("Restart cycle: %ld\n", (long)m);
  printf("Number of iterations: %ld\n", (long)Iter);
  printf("Actural relative residual: %g\n", NormR/NormB);
  printf("Relative residual reported by the solver: %g\n", ResHistory[Iter-1]/NormB);

  const char *out_file4 = "Test_GMRES.relres4";
  fid = fopen(out_file4, "w");
  for (INTEGER i = 0; i < Iter; i++) {
    fprintf(fid, "%g\n", ResHistory[i]/NormB);
  }
  fclose(fid);
  printf("Relative residual history is written to file %s. Use the Matlab program Test_GMRES.m to plot it out!\n", out_file4);

  return 0;

}
