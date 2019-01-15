// This program performs GP analysis (MLE, kriging) by using the RLCM
// method.
//
// The current implementation uses a 2-dimensional grid (DPointArray),
// a chordal Matern kernel, and 3 parameters for the kernel: alpha,
// ell, and tau. The actual nugget is lambda = 10^tau. The global
// scaling s = 10^alpha. The smoothness nu is always 1/2.
//
// The code is written for the 500 mb temperature data set, which has
// been demeaned and each time step is an independent sample. The
// ordering of the data is lat * lon * time, where the first dimension
// lat varies the fastest.
//
// The current implementation supports parallelism. One may set the
// USE_OPENBLAS flag at compile-time so that all the BLAS and LAPACK
// routines are threaded. Additionally, one may set the USE_OPENMP
// flag so that other parts are threaded.
//
// Compile-time macros:
//
//   USE_OPENBLAS:     Either use this macro or not (no value)
//   USE_OPENMP:       Either use this macro or not (no value)
//
// Usage:
//
//   GP_RLCM_Multiple_Samples.ex NumThreads DataFile Dim Lower Upper
//   NumSamples Num_alpha List_alpha Num_ell List_ell Num_tau List_tau r
//   DiagCorrect Seed IsCheckFiniteDiff DiffStepSize OutputRandomField
//   [RandomFieldFileBasename] OutputLogLik [LogLikFileName]
//   IsComputeFisher OutputFisher [FisherFileName]
//   OutputKrigedRandomField [KrigedRandomFieldFileBasename]
//   OutputPredictions [PredictionsFileName]
//
//   NumThreads:   Number of threads
//   DataFile:     File name (including path) of the data file
//   Dim:          (Grid) Size of each dimension of the grid. Length d = 2
//   Lower:        (Grid) Lower end of the physical grid. Length d = 2
//   Upper:        (Grid) Upper end of the physical grid. Length d = 2
//   NumSamples:   Number of samples (time steps)
//   Num_alpha:    (Param grid search) Length of list of alpha
//   List_alpha:   (Param grid search) List of alpha
//   Num_ell:      (Param grid search) Length of list of ell
//   List_ell:     (Param grid search) List of ell
//   Num_tau:      (Param grid search) Length of list of tau
//   List_tau:     (Param grid search) List of tau
//   r:            (Matrix structure) Rank
//   DiagCorrect:  (Matrix structure) Diagonal correction, e.g., 1e-8
//   Seed:         If >= 0, will use this value as the seed for RNG;
//                 otherwise, use the current time to seed the RNG.
//   IsCheckFiniteDiff:       If > 0, inspect what diff step size is appropriate
//   DiffStepSize:            An array of NUM_PARAM numbers
//   OutputRandomField:       If > 0, output random field to file
//   RandomFieldFileBasename:       If above > 0, file basename (no ext)
//   OutputLogLik:            If > 0, output calculated loglik's to file
//   LogLikFileName:                If above > 0, file name (including ext)
//   IsComputeFisher:         If > 0, compute Fisher and related quantities
//   OutputFisher:            If > 0 and above > 0, output Fisher to file
//   FisherFileName:                If above two > 0, file name (including ext)
//   OutputKrigedRandomField: If > 0, output kriged random field to file
//   KrigedRandomFieldFileBasename: If above > 0, file basename (no ext)
//   OutputPredictions:       If > 0, output predictions to file
//   PredictionsFileName:           If above > 0, file name (including ext)
//
// For file formats of the output files, see GP_Common.hpp.

#include "GP_Common.hpp"

#define NUM_PARAM 3

int main(int argc, char **argv) {

  //---------- Parameters from command line --------------------

  INTEGER idx = 1;
  int NumThreads = atoi(argv[idx++]);

  // Data
  char *DataFile = argv[idx++];

  // Grid
  INTEGER d = 2; // d must be 2
  INTEGER *Dim = NULL;
  New_1D_Array<INTEGER, INTEGER>(&Dim, d);
  for (INTEGER i = 0; i < d; i++) {
    Dim[i] = String2Integer(argv[idx++]);
  }
  // Dim[0] is the number of latitudes and Dim[1] is the number of
  // longitudes. Dim[1] must be even to ensure that after
  // downsampling, the longitudes are still evenly spaced.
  if (Dim[1] % 2 != 0) {
    printf("GP_RLCM_Multiple_Samples: Error. Dim[1] (number of longitudes) must be even!\n"); fflush(stdout);
    return 0;
  }
  double *Lower = NULL;
  New_1D_Array<double, INTEGER>(&Lower, d);
  for (INTEGER i = 0; i < d; i++) {
    Lower[i] = atof(argv[idx++]);
  }
  // Lower[1] is the lower end of longitude. It must be -pi+2*pi/Dim[1]
  double LowerLon = -M_PI + 2*M_PI/Dim[1];
  if (fabs(Lower[1] - LowerLon) > EPSx10) {
    printf("GP_RLCM_Multiple_Samples: Warning. Lower end of longitude %.16e is reset to %.16e\n", Lower[1], LowerLon); fflush(stdout);
    Lower[1] = LowerLon;
  }
  double *Upper = NULL;
  New_1D_Array<double, INTEGER>(&Upper, d);
  for (INTEGER i = 0; i < d; i++) {
    Upper[i] = atof(argv[idx++]);
  }
  // Upper[1] is the upper end of longitude. It must be pi
  double UpperLon = M_PI;
  if (fabs(Upper[1] - UpperLon) > EPSx10) {
    printf("GP_RLCM_Multiple_Samples: Warning. Upper end of longitude %.16e is reset to %.16e\n", Upper[1], UpperLon); fflush(stdout);
    Upper[1] = UpperLon;
  }
  INTEGER NumSamples = String2Integer(argv[idx++]);

  // Kernel
  double nu = 1.5; // nu is always 1.5

  // Param grid search
  INTEGER Num_alpha = String2Integer(argv[idx++]);
  double *List_alpha = NULL;
  New_1D_Array<double, INTEGER>(&List_alpha, Num_alpha);
  for (INTEGER i = 0; i < Num_alpha; i++) {
    List_alpha[i] = atof(argv[idx++]);
  }
  INTEGER Num_ell = String2Integer(argv[idx++]);
  double *List_ell = NULL;
  New_1D_Array<double, INTEGER>(&List_ell, Num_ell);
  for (INTEGER i = 0; i < Num_ell; i++) {
    List_ell[i] = atof(argv[idx++]);
  }
  INTEGER Num_tau = String2Integer(argv[idx++]);
  double *List_tau = NULL;
  New_1D_Array<double, INTEGER>(&List_tau, Num_tau);
  for (INTEGER i = 0; i < Num_tau; i++) {
    List_tau[i] = atof(argv[idx++]);
  }

  // Matrix structure
  INTEGER r = String2Integer(argv[idx++]);
  double mDiagCorrect = atof(argv[idx++]);

  // RNG
  INTEGER sSeed = String2Integer(argv[idx++]);
  unsigned Seed;
  if (sSeed < 0) {
    Seed = (unsigned)time(NULL);
  }
  else {
    Seed = (unsigned)sSeed;
  }

  // Finite difference
  bool IsCheckFiniteDiff = atoi(argv[idx++]) ? true : false;
  double DiffStepSize[NUM_PARAM];
  for (INTEGER i = 0; i < NUM_PARAM; i++) {
    DiffStepSize[i] = atof(argv[idx++]);
  }

  // Diagnostics
  bool OutputRandomField = atoi(argv[idx++]) ? true : false;
  char *RandomFieldFileBasename = NULL;
  if (OutputRandomField) {
    RandomFieldFileBasename = argv[idx++];
  }
  bool OutputLogLik = atoi(argv[idx++]) ? true : false;
  char *LogLikFileName = NULL;
  if (OutputLogLik) {
    LogLikFileName = argv[idx++];
  }
  bool IsComputeFisher = atoi(argv[idx++]) ? true : false;
  bool OutputFisher = atoi(argv[idx++]) ? true : false;
  char *FisherFileName = NULL;
  if (OutputFisher) {
    FisherFileName = argv[idx++];
  }
  bool OutputKrigedRandomField = atoi(argv[idx++]) ? true : false;
  char *KrigedRandomFieldFileBasename = NULL;
  if (OutputKrigedRandomField) {
    KrigedRandomFieldFileBasename = argv[idx++];
  }
  bool OutputPredictions = atoi(argv[idx++]) ? true : false;
  char *PredictionsFileName = NULL;
  if (OutputPredictions) {
    PredictionsFileName = argv[idx++];
  }

  //---------- Threading --------------------

#ifdef USE_OPENBLAS
  openblas_set_num_threads(NumThreads);
#elif defined USE_OPENMP
  omp_set_num_threads(NumThreads);
#else
  NumThreads = 1; // To avoid compiler warining of unused variable
#endif

  //---------- Main Computation --------------------

  // Timing
  PREPARE_CLOCK(true);

  // Seed the RNG
  srandom(Seed);

  // Generate regular grid X
  DPointArray X;
  X.SetRegularGrid(d, Dim, Lower, Upper);

  // Load in data sample y. Treat it as a long vector (whose length
  // should be Dim[0]*Dim[1]*NumSamples) and ignore the time steps.
  DVector y;
  if (Raw_binary::ReadData(DataFile, y) == false) {
    printf("GP_RLCM_Multiple_Samples: Error reading input data!\n");
    return 0;
  }
  if (y.GetN() != Dim[0]*Dim[1]*NumSamples) {
    printf("GP_RLCM_Multiple_Samples: Error. Size of data is not Dim[0] * Dim[1] * NumSamples!\n");
    return 0;
  }

  // Output random field to file. For the sake of practical purposes,
  // there is no need to output the whole data. We already have the
  // data. We output only the first sample for debugging purposes.
  if (OutputRandomField) {
    DVector y1;
    y.GetBlock(0, Dim[0]*Dim[1], y1);
    WriteRandomFieldToFile(y1, d, Dim, Lower, Upper, 0, NULL,
                           RandomFieldFileBasename, "GP_RLCM_Multiple_Samples");
  }

  // Train/test split. Ytrain is a matrix, where each column is an
  // independent sample. We krige only the first sample. Therefore,
  // ytest is a vector. IdxTrain and IdxTest are vectors, which take
  // care of only the first time step (because for different time
  // steps each vector is the same).
  DPointArray Xtrain, Xtest;
  DMatrix Ytrain;
  DVector ytest;
  INTEGER *IdxTrain = NULL, *IdxTest = NULL;
  TrainTestSplitMultipleSamples(X, y, Dim, NumSamples, Xtrain, Xtest,
                                Ytrain, ytest, &IdxTrain, &IdxTest);

  // Save some memory
  X.ReleaseAllMemory();
  y.ReleaseAllMemory();

  // Kernel matrix
  CMatrix Ktrain;
  INTEGER Ntrain = Xtrain.GetN();
  INTEGER *PermXtrain = NULL;
  New_1D_Array<INTEGER, INTEGER>(&PermXtrain, Ntrain);
  INTEGER NumLevel = (INTEGER)log2((double)Ntrain/r); //excluding the root level
  Ktrain.BuildTree<DPoint, DPointArray>(Xtrain, PermXtrain, NULL, r,
                                        NumLevel, mDiagCorrect, Seed, BBOX);

  // Xtrain is permuted in Ktrain.BuildTree(). Accordingly permute
  // Ytrain
  for (INTEGER i = 0; i < NumSamples; i++) {
    DVector ytrain;
    Ytrain.GetColumn(i, ytrain);
    ytrain.Permute(PermXtrain, Ntrain);
    Ytrain.SetColumn(i, ytrain);
  }

  // MLE through grid search
  INTEGER ListLength = Num_alpha * Num_ell * Num_tau;
  double *mLogLik = NULL;
  New_1D_Array<double, INTEGER>(&mLogLik, ListLength);
  double **ListParam = NULL;
  New_2D_Array<double, INTEGER, INTEGER>(&ListParam, ListLength, NUM_PARAM);
  MLE_RLCM<ChordalMatern, DPoint, DPointArray> mMLE;
  START_CLOCK;
  INTEGER t = 0;
  for (INTEGER i = 0; i < Num_alpha; i++) {
    double alpha = List_alpha[i];
    for (INTEGER j = 0; j < Num_ell; j++) {
      double ell = List_ell[j];
      for (INTEGER k = 0; k < Num_tau; k++) {
        double tau = List_tau[k];
        ListParam[t][0] = alpha;
        ListParam[t][1] = ell;
        ListParam[t][2] = tau;
        double s        = pow(10.0, alpha);
        double lambda   = pow(10.0, tau);
        ChordalMatern mKernel(s, nu, ell);
        mLogLik[t] = mMLE.LogLik(Ktrain, Xtrain, Ytrain, mKernel, lambda);
        printf("MLE_RLCM: Grid search alpha = %g, ell = %g, tau = %g, loglik = %.16e\n", alpha, ell, tau, mLogLik[t]); fflush(stdout);
        t++;
      }
    }
  }
  END_CLOCK;
  double TimeMLE = ELAPSED_TIME;

  double HatParam[NUM_PARAM];
  double MaxLogLik;
  EstimatedParam(NUM_PARAM, ListLength, ListParam, mLogLik,
                 HatParam, MaxLogLik);
  double hat_alpha = HatParam[0];
  double hat_ell   = HatParam[1];
  double hat_tau   = HatParam[2];
  printf("GP_RLCM_Multiple_Samples: Estimated alpha = %g, ell = %g, tau = %g, max loglik = %.16e, MLE time = %gs\n", hat_alpha, hat_ell, hat_tau, MaxLogLik, TimeMLE); fflush(stdout);

  // Output all loglik's to file
  if (OutputLogLik) {
    WriteLogLikToFile(NUM_PARAM, ListLength, ListParam, mLogLik,
                      LogLikFileName, "GP_RLCM_Multiple_Samples");
  }

  // Finite difference
  if (IsCheckFiniteDiff) {

    double delta = 1e-3;
    double fac = 2.0;
    INTEGER NumStep = 10;
    for (INTEGER i = 0; i < NUM_PARAM; i++) {
      printf("\nNumerical differentiation for param #%ld:\n", (long)i); fflush(stdout);
      printf("Step size                First difference          Second difference\n"); fflush(stdout);
      for (INTEGER j = 0; j < NumStep; j++) {

        double CenterParam[NUM_PARAM];
        memcpy(CenterParam, HatParam, NUM_PARAM*sizeof(double));
        double epsilon = delta / pow(fac, j);
        CenterParam[i] += epsilon;
        double alpha   = CenterParam[0];
        double ell     = CenterParam[1];
        double tau     = CenterParam[2];
        double s       = pow(10.0, alpha);
        double lambda  = pow(10.0, tau);
        ChordalMatern mKernel(s, nu, ell);
        double mLogLikp = mMLE.LogLik(Ktrain, Xtrain, Ytrain, mKernel, lambda);

        memcpy(CenterParam, HatParam, NUM_PARAM*sizeof(double));
        CenterParam[i] -= epsilon;
        alpha  = CenterParam[0];
        ell    = CenterParam[1];
        tau    = CenterParam[2];
        s      = pow(10.0, alpha);
        lambda = pow(10.0, tau);
        mKernel.Init(s, nu, ell);
        double mLogLikm = mMLE.LogLik(Ktrain, Xtrain, Ytrain, mKernel, lambda);

        double FirstDiff = (mLogLikp - mLogLikm) / (2.0*epsilon);
        double SecondDiff = (mLogLikp + mLogLikm - 2.0*MaxLogLik) /
          (epsilon*epsilon);

        printf("%.16e   %+.16e   %+.16e\n", epsilon, FirstDiff, SecondDiff); fflush(stdout);
      }
    }

  }

  // Fisher information
  if (IsComputeFisher) {

    INTEGER ListLength2 = 2 * NUM_PARAM + 4 * NUM_PARAM*(NUM_PARAM-1)/2;
    double *mLogLik2 = NULL;
    New_1D_Array<double, INTEGER>(&mLogLik2, ListLength2);
    double **ListParam2 = NULL;
    New_2D_Array<double, INTEGER, INTEGER>(&ListParam2, ListLength2, NUM_PARAM);
    PrepareListParamForFisher(NUM_PARAM, ListLength2, ListParam2, HatParam,
                              DiffStepSize);

    for (INTEGER i = 0; i < ListLength2; i++) {
      double alpha  = ListParam2[i][0];
      double ell    = ListParam2[i][1];
      double tau    = ListParam2[i][2];
      double s      = pow(10.0, alpha);
      double lambda = pow(10.0, tau);
      ChordalMatern mKernel(s, nu, ell);
      mLogLik2[i] = mMLE.LogLik(Ktrain, Xtrain, Ytrain, mKernel, lambda);
      printf("MLE_RLCM: Fisher info alpha = %g, ell = %g, tau = %g, loglik = %.16e\n", alpha, ell, tau, mLogLik2[i]); fflush(stdout);
    }

    DMatrix Fisher, Cov;
    DVector Stderr;
    ComputeFisher(NUM_PARAM, ListLength2, ListParam2, mLogLik2,
                  HatParam, MaxLogLik, DiffStepSize, Fisher, Cov, Stderr);
    double *mStderr = Stderr.GetPointer();
    printf("GP_RLCM_Multiple_Samples: Stderr alpha = %g, ell = %g, tau = %g\n", mStderr[0], mStderr[1], mStderr[2]); fflush(stdout);

    Delete_1D_Array<double>(&mLogLik2);
    Delete_2D_Array<double, INTEGER>(&ListParam2, ListLength2);

    // Output Fisher to file
    if (OutputFisher) {
      WriteFisherToFile(Fisher, FisherFileName, "GP_RLCM_Multiple_Samples");
    }

  }

  // Use estimated parameters to do kriging. Krige only the first time
  // step.
  Kriging_RLCM<ChordalMatern, DPoint, DPointArray> mKriging;
  double hat_s      = pow(10.0, hat_alpha);
  double hat_lambda = pow(10.0, hat_tau);
  ChordalMatern mKernel2(hat_s, nu, hat_ell);
  DVector ytest_predict(Xtest.GetN());
  DVector ytest_stddev(Xtest.GetN());
  START_CLOCK;
  mKriging.Train(Ktrain, Xtrain, mKernel2, hat_lambda);
  DVector ytrain;
  Ytrain.GetColumn(0, ytrain);
  mKriging.Test(Ktrain, Xtrain, Xtest, ytrain, mKernel2, hat_lambda,
                ytest_predict, ytest_stddev);
  END_CLOCK;
  double TimeKriging = ELAPSED_TIME;
  printf("GP_RLCM_Multiple_Samples: Kriging time = %gs\n", TimeKriging); fflush(stdout);

  // Output kriged random field to file
  if (OutputKrigedRandomField) {
    DVector ytrain_recover_order;
    Ytrain.GetColumn(0, ytrain_recover_order);
    ytrain_recover_order.iPermute(PermXtrain, Ntrain);
    DVector y_all;
    AssembleY(ytrain_recover_order, ytest_predict, IdxTrain, IdxTest, y_all);
    WriteRandomFieldToFile(y_all, d, Dim, Lower, Upper, NUM_PARAM, HatParam,
                           KrigedRandomFieldFileBasename,
                           "GP_RLCM_Multiple_Samples");
  }

  // Output predictions to file
  if (OutputPredictions) {
    WritePredictionsToFile(ytest, ytest_predict, ytest_stddev,
                           PredictionsFileName, "GP_RLCM_Multiple_Samples");
  }

  //---------- Clean up --------------------

  Delete_1D_Array<INTEGER>(&Dim);
  Delete_1D_Array<double>(&Lower);
  Delete_1D_Array<double>(&Upper);
  Delete_1D_Array<double>(&List_alpha);
  Delete_1D_Array<double>(&List_ell);
  Delete_1D_Array<double>(&List_tau);
  Delete_1D_Array<INTEGER>(&PermXtrain);
  Delete_1D_Array<INTEGER>(&IdxTrain);
  Delete_1D_Array<INTEGER>(&IdxTest);
  Delete_1D_Array<double>(&mLogLik);
  Delete_2D_Array<double, INTEGER>(&ListParam, ListLength);

  return 0;

}
