// This program performs kernel ridge regression by using the
// BlockDiag method. The actual machine learning task can be either
// regression, binary classification, or multiclass classification.
//
// The BlockDiag method requires a parameter N0, which is the block
// size. Meanwhile, for large data sets, this program may not be able
// to handle the test data all at once. Hence, a parameter Budget is
// used to split the test data into batches, each handled at a
// time. The number of test data in one batch is approximately
// Budget/N0.
//
// The BlockDiag method also offers the option Par, which specifies
// the method for partitioning the data set when building the
// blocks. RAND is much more efficient than PCA.
//
// The current implementation supports the following kernels:
// isotropic Gaussian, isotropic Laplace, product Laplace, and inverse
// multiquadric. For their definitions, see the corresponding .hpp
// files under ${CMATRIX_DIR}/src/Kernels/. The required paramters are
// sigma (bandwidth) and lambda (regularization). The global scaling s
// is always 1.
//
// The current implementation supports two data types for storing
// points: dense (DPoint) and sparse (SPoint). The supported data file
// format is LibSVM format. Even though LibSVM format always treats
// the data sparse, some data are in fact (almost) fully dense. The
// linear algebra is different for dense and sparse points. From a
// practical stand point, using the dense format is often more
// efficient, unless memory explodes.
//
// A note on LibSVM format: The LibSVM reader implemented in this
// software treats the attribute index starting from 1. For binary
// classifications, the labels must be +1/-1. For multiclass
// classifications, the labels must start from 0 and must be
// consecutive integers.
//
// The current implementation supports parallelism. One may set the
// USE_OPENBLAS flag at compile-time so that all the BLAS and LAPACK
// routines are threaded. Additionally, one may set the USE_OPENMP
// flag so that other parts are threaded.
//
// Compile-time macros:
//
//   KernelType:       One of IsotropicGaussian, IsotropicLaplace,
//                     ProductLaplace, InvMultiquadric
//   PointFormat:      One of DPoint, SPoint
//   PointArrayFormat: One of DPointArray, SPointArray. Must be
//                     consistent with PointFormat
//   USE_OPENBLAS:     Either use this macro or not (no value)
//   USE_OPENMP:       Either use this macro or not (no value)
//
// Usage:
//
//   KRR_BlockDiag.ex NumThreads FileTrain FileTest NumClasses d Seed
//   N0 Budget Num_sigma List_sigma Num_lambda List_lambda Par
//
//   NumThreads:  Number of threads
//   FileTrain:   File name (including path) of the data file for
//                training
//   FileTest:    File name (including path) of the data file for
//                testing
//   NumClasses:  1 for regression, 2 for classification, and intergers
//                larger than 2 for multiclass classification. In the
//                last case, the integer denotes the number of classes
//   d:           Dimension of the data
//   Seed:        If >= 0, will use this value as the seed for RNG;
//                otherwise, use the current time to seed the RNG.
//   N0:          Block size
//   Budget:      Budget. If <=0, test data are handled all at once
//   Num_sigma:   Length of list of sigma
//   List_sigma:  List of sigma
//   Num_lambba:  Length of list of lambda
//   List_lambda: List of lambda
//   Par:         Partitioning method. One of RAND, PCA

#include "KRR_Common.hpp"

int main(int argc, char **argv) {

  //---------- Parameters from command line --------------------

  INTEGER idx = 1;
  int NumThreads = atoi(argv[idx++]);               // Number of threads
  char *FileTrain = argv[idx++];                    // Training data
  char *FileTest = argv[idx++];                     // Testing data
  INTEGER NumClasses = String2Integer(argv[idx++]); // Number of classes
  INTEGER d = String2Integer(argv[idx++]);          // Data dimension
  INTEGER sSeed = String2Integer(argv[idx++]);      // Seed for randomization
  unsigned Seed;
  if (sSeed < 0) {
    Seed = (unsigned)time(NULL);
  }
  else {
    Seed = (unsigned)sSeed;
  }
  INTEGER N0 = String2Integer(argv[idx++]);         // Block size
  INTEGER Budget = String2Integer(argv[idx++]);     // Budget
  INTEGER Num_sigma = String2Integer(argv[idx++]);  // Length of list of sigma
  double *List_sigma = NULL;                        // List of sigma
  New_1D_Array<double, INTEGER>(&List_sigma, Num_sigma);
  for (INTEGER i = 0; i < Num_sigma; i++) {
    List_sigma[i] = atof(argv[idx++]);
  }
  INTEGER Num_lambda = String2Integer(argv[idx++]); // Length of list of lambda
  double *List_lambda = NULL;                       // List of lambda
  New_1D_Array<double, INTEGER>(&List_lambda, Num_lambda);
  for (INTEGER i = 0; i < Num_lambda; i++) {
    List_lambda[i] = atof(argv[idx++]);
  }
  char *ParString = argv[idx++];                    // Partitioning method
  PartMethod Par;
  if (strcmp(ParString, "RAND") == 0) {
    Par = RAND;
  }
  else if (strcmp(ParString, "PCA") == 0) {
    Par = PCA;
  }
  else {
    printf("KRR_BlockDiag. Error: Unkown partitioning method!\n");
    return 0;
  }

  //---------- Read in data --------------------

  PointArrayFormat Xtrain; // Training points
  PointArrayFormat Xtest;  // Testing points
  DVector ytrain;          // Training labels
  DVector ytest;           // Testing labels (ground truth)
  DVector ytest_predict;   // Predictions

  if (LibSVM_IO::ReadData(FileTrain, Xtrain, ytrain, d) == false) {
    printf("KRR_BlockDiag: Error reading input data!\n");
    return 0;
  }
  if (LibSVM_IO::ReadData(FileTest, Xtest, ytest, d) == false) {
    printf("KRR_BlockDiag: Error reading input data!\n");
    return 0;
  }

  // For multiclass classification, need to convert a single vector
  // ytrain to a matrix Ytrain. The "predictions" are stored in the
  // corresponding matrix Ytest_predict. The vector ytest_predict is
  // unused.
  DMatrix Ytrain;
  DMatrix Ytest_predict;
  if (NumClasses > 2) {
    ConvertYtrain(ytrain, Ytrain, NumClasses);
  }

  //---------- Threading --------------------

#ifdef USE_OPENBLAS
  openblas_set_num_threads(NumThreads);
#elif defined USE_OPENMP
  omp_set_num_threads(NumThreads);
#else
  NumThreads = 1; // To avoid compiler warining of unused variable
#endif

  //---------- Main computation --------------------

  PREPARE_CLOCK(true);

  KRR_BlockDiag<KernelType, PointFormat, PointArrayFormat> mKRR_BlockDiag;
  INTEGER *Perm = NULL, *iPerm = NULL;
  INTEGER N = Xtrain.GetN();
  New_1D_Array<INTEGER, INTEGER>(&Perm, N);
  New_1D_Array<INTEGER, INTEGER>(&iPerm, N);

  // Pre-training
  START_CLOCK;
  double MemEst;
  if (NumClasses <= 2) {
    MemEst = mKRR_BlockDiag.PreTrain(Xtrain, ytrain, Perm, iPerm, N0, Seed,
                                     Par);
  }
  else {
    MemEst = mKRR_BlockDiag.PreTrain(Xtrain, Ytrain, Perm, iPerm, N0, Seed,
                                     Par);
  }
  END_CLOCK;
  double TimePreTrain = ELAPSED_TIME;
  printf("KRR_BlockDiag: pretrain time = %g\n", TimePreTrain); fflush(stdout);

  double s = 1.0; // s is always 1
  for (INTEGER j = 0; j < Num_sigma; j++) {
    double sigma = List_sigma[j];
    for (INTEGER k = 0; k < Num_lambda; k++) {
      double lambda = List_lambda[k];

      KernelType mKernel(s, sigma);

      // Training
      START_CLOCK;
      mKRR_BlockDiag.Train(Xtrain, mKernel, lambda);
      END_CLOCK;
      double TimeTrain = ELAPSED_TIME;

      // Testing
      START_CLOCK;
      if (NumClasses <= 2) {
        if (Budget <= 0) {
          mKRR_BlockDiag.Test(Xtrain, Xtest, ytrain, mKernel, ytest_predict);
        }
        else {
          ytest_predict.Init(Xtest.GetN());
          INTEGER BatchSize = (INTEGER)ceil((double)Budget/N0);
          INTEGER RemainSize = Xtest.GetN();
          INTEGER ThisBatchSize;
          INTEGER ThisBatchStart = 0;
          while (RemainSize > 0) {
            ThisBatchSize = BatchSize < RemainSize ? BatchSize : RemainSize;
            PointArrayFormat Xtest_batch;
            DVector ytest_predict_batch;
            Xtest.GetSubset(ThisBatchStart, ThisBatchSize, Xtest_batch);
            mKRR_BlockDiag.Test(Xtrain, Xtest_batch, ytrain, mKernel,
                                ytest_predict_batch);
            ytest_predict.SetBlock(ThisBatchStart, ThisBatchSize,
                                   ytest_predict_batch);
            RemainSize -= ThisBatchSize;
            ThisBatchStart += ThisBatchSize;
          }
        }
      }
      else {
        if (Budget <= 0) {
          mKRR_BlockDiag.Test(Xtrain, Xtest, Ytrain, mKernel, Ytest_predict);
        }
        else {
          Ytest_predict.Init(Xtest.GetN(), NumClasses);
          INTEGER BatchSize = (INTEGER)ceil((double)Budget/N0);
          INTEGER RemainSize = Xtest.GetN();
          INTEGER ThisBatchSize;
          INTEGER ThisBatchStart = 0;
          while (RemainSize > 0) {
            ThisBatchSize = BatchSize < RemainSize ? BatchSize : RemainSize;
            PointArrayFormat Xtest_batch;
            DMatrix Ytest_predict_batch;
            Xtest.GetSubset(ThisBatchStart, ThisBatchSize, Xtest_batch);
            mKRR_BlockDiag.Test(Xtrain, Xtest_batch, Ytrain, mKernel,
                                Ytest_predict_batch);
            Ytest_predict.SetBlock(ThisBatchStart, ThisBatchSize, 0, NumClasses,
                                   Ytest_predict_batch);
            RemainSize -= ThisBatchSize;
            ThisBatchStart += ThisBatchSize;
          }
        }
      }
      END_CLOCK;
      double TimeTest = ELAPSED_TIME;

      // Performance
      double Perf = 0.0;
      if (NumClasses <= 2) {
        Perf = Performance(ytest, ytest_predict, NumClasses);
      }
      else {
        Perf = Performance(ytest, Ytest_predict, NumClasses);
      }
      printf("KRR_BlockDiag: N0 = %ld, param = %g %g, perf = %g, time = %g %g, mem = %g\n", (long)N0, sigma, lambda, Perf, TimeTrain, TimeTest, MemEst); fflush(stdout);

    }
  }

  //---------- Clean up --------------------

  Delete_1D_Array<INTEGER>(&Perm);
  Delete_1D_Array<INTEGER>(&iPerm);
  Delete_1D_Array<double>(&List_sigma);
  Delete_1D_Array<double>(&List_lambda);

  return 0;

}
