// Unit test of the class SPoint.

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
  printf("|  Test_SPoint.cpp                                              |\n");
  printf("+---------------------------------------------------------------+\n");

  // Test InProd() ------------------------------------------------------------
  //
  // x = [ 0.0; 0.4; 0.3; 0.0 ]
  // y = [ 0.1; 0.0; -0.2; 1.0 ]
  // inprod = -0.06
  INTEGER d = 4;
  INTEGER nnzx = 2;
  INTEGER idxx[] = { 1, 2 };
  double sx[] = { 0.4, 0.3 };
  SPoint SPx(d, nnzx);
  memcpy(SPx.GetPointerIdx(), idxx, nnzx*sizeof(INTEGER));
  memcpy(SPx.GetPointerX(), sx, nnzx*sizeof(double));
  double dy[] = { 0.1, 0.0, -0.2, 1.0 };
  DPoint DPy(d);
  memcpy(DPy.GetPointer(), dy, d*sizeof(double));
  double inprod = SPx.InProd(DPy);
  double inprod_truth = -0.06;
  double derr = fabs(inprod-inprod_truth);
  printf("Test InProd() : Discrepancy in return value %g\n", derr);

  INTEGER nnzy = 3;
  INTEGER idxy[] = { 0, 2, 3 };
  double sy[] = { 0.1, -0.2, 1.0 };
  SPoint SPy(d, nnzy);
  memcpy(SPy.GetPointerIdx(), idxy, nnzy*sizeof(INTEGER));
  memcpy(SPy.GetPointerX(), sy, nnzy*sizeof(double));
  inprod = SPx.InProd(SPy);
  derr = fabs(inprod-inprod_truth);
  printf("Test InProd() : Discrepancy in return value %g\n", derr);

  // Test Dist1() -------------------------------------------------------------
  //
  // Use the above x and y
  // dist1 = 2.0
  double dist1 = SPx.Dist1(DPy);
  double dist1_truth = 2.0;
  derr = fabs(dist1-dist1_truth);
  printf("Test Dist1() : Discrepancy in return value %g\n", derr);

  dist1 = SPx.Dist1(SPy);
  derr = fabs(dist1-dist1_truth);
  printf("Test Dist1() : Discrepancy in return value %g\n", derr);

  // Test Dist1() -------------------------------------------------------------
  //
  // Use the above x and y
  // sigma = [ 2.0; 1.0; 0.5; 2.0 ]
  // dist1 = 1.95
  double sigma[] = { 2.0, 1.0, 0.5, 2.0 };
  dist1 = SPx.Dist1(DPy, sigma);
  dist1_truth = 1.95;
  derr = fabs(dist1-dist1_truth);
  printf("Test Dist1() : Discrepancy in return value %g\n", derr);

  dist1 = SPx.Dist1(SPy, sigma);
  derr = fabs(dist1-dist1_truth);
  printf("Test Dist1() : Discrepancy in return value %g\n", derr);

  // Test Dist2() -------------------------------------------------------------
  //
  // Use the above x and y
  // dist1 = 1.42
  double dist2 = SPx.Dist2(DPy);
  double dist2_truth = 1.42;
  derr = fabs(dist2-dist2_truth);
  printf("Test Dist2() : Discrepancy in return value %g\n", derr);

  dist2 = SPx.Dist2(SPy);
  derr = fabs(dist2-dist2_truth);
  printf("Test Dist2() : Discrepancy in return value %g\n", derr);
 
  // Test Dist2() -------------------------------------------------------------
  //
  // Use the above x, y and sigma
  // dist1 = 1.4125
  dist2 = SPx.Dist2(DPy, sigma);
  dist2_truth = 1.4125;
  derr = fabs(dist2-dist2_truth);
  printf("Test Dist2() : Discrepancy in return value %g\n", derr);

  dist2 = SPx.Dist2(SPy, sigma);
  derr = fabs(dist2-dist2_truth);
  printf("Test Dist2() : Discrepancy in return value %g\n", derr);
 
  // Test KernelFuncChi2() ----------------------------------------------------
  //
  // x = [ 0.0; 0.0; 0.3; 0.0 ]
  // y = [ 0.1; 0.0; 1.2; 0.1 ]
  // kval = 0.48
  nnzx = 1;
  idxx[0] = 2;
  sx[0] = 0.3;
  SPx.Init(d, nnzx);
  memcpy(SPx.GetPointerIdx(), idxx, nnzx*sizeof(INTEGER));
  memcpy(SPx.GetPointerX(), sx, nnzx*sizeof(double));
  dy[0] = 0.1; dy[1] = 0.0; dy[2] = 1.2; dy[3] = 0.1;
  DPy.Init(d);
  memcpy(DPy.GetPointer(), dy, d*sizeof(double));
  double kval = SPx.KernelFuncChi2(DPy);
  double kval_truth = 0.48;
  derr = fabs(kval-kval_truth);
  printf("Test KernelFuncChi2() : Discrepancy in return value %g\n", derr);

  nnzy = 3;
  idxy[0] = 0; idxy[1] = 2; idxy[2] = 3;
  sy[0] = 0.1; sy[1] = 1.2; sy[2] = 0.1;
  SPy.Init(d, nnzy);
  memcpy(SPy.GetPointerIdx(), idxy, nnzy*sizeof(INTEGER));
  memcpy(SPy.GetPointerX(), sy, nnzy*sizeof(double));
  kval = SPx.KernelFuncChi2(SPy);
  derr = fabs(kval-kval_truth);
  printf("Test KernelFuncChi2() : Discrepancy in return value %g\n", derr);

  return 0;

}
