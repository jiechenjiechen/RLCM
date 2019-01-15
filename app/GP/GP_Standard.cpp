// This program performs GP analysis (sampling, MLE, kriging) by using
// the Standard method.
//
// The current implementation uses a d-dimensional grid (DPointArray),
// an isotropic Matern kernel, and 3 parameters for the kernel: ell,
// nu, and tau. The actual nugget is lambda = 10^tau. The global
// scaling s is always 1.
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
//   GP_Standard.ex NumThreads d Dim Lower Upper ell nu tau Num_ell
//   List_ell Num_nu List_nu Num_tau List_tau Seed PortionTrain
//   IsCheckFiniteDiff DiffStepSize OutputRandomField
//   [RandomFieldFileBasename] OutputLogLik [LogLikFileName]
//   IsComputeFisher OutputFisher [FisherFileName]
//   OutputKrigedRandomField [KrigedRandomFieldFileBasename]
//   OutputPredictions [PredictionsFileName]
//
//   NumThreads:   Number of threads
//   d:            (Grid) Dimension of the data
//   Dim:          (Grid) Size of each dimension of the grid. Length d
//   Lower:        (Grid) Lower end of the physical grid. Length d
//   Upper:        (Grid) Upper end of the physical grid. Length d
//   ell:          (True param) Param of Matern kernel
//   nu:           (True param) Param of Matern kernel
//   tau:          (True param) Param of Matern kernel
//   Num_ell:      (Param grid search) Length of list of ell
//   List_ell:     (Param grid search) List of ell
//   Num_nu:       (Param grid search) Length of list of nu
//   List_nu:      (Param grid search) List of nu
//   Num_tau:      (Param grid search) Length of list of tau
//   List_tau:     (Param grid search) List of tau
//   Seed:         If >= 0, will use this value as the seed for RNG;
//                 otherwise, use the current time to seed the RNG.
//   PortionTrain:            Portion of grid for training
//   IsCheckFiniteDiff:       If > 0, inspect what diff step size is appropriate
//   DiffStepSize:            An array of NUM_PARAM numbers
//   TimingEachLogLik:        If > 0, record time of each loglik calculation
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

  // Grid
  INTEGER d = String2Integer(argv[idx++]);
  INTEGER *Dim = NULL;
  New_1D_Array<INTEGER, INTEGER>(&Dim, d);
  for (INTEGER i = 0; i < d; i++) {
    Dim[i] = String2Integer(argv[idx++]);
  }
  double *Lower = NULL;
  New_1D_Array<double, INTEGER>(&Lower, d);
  for (INTEGER i = 0; i < d; i++) {
    Lower[i] = atof(argv[idx++]);
  }
  double *Upper = NULL;
  New_1D_Array<double, INTEGER>(&Upper, d);
  for (INTEGER i = 0; i < d; i++) {
    Upper[i] = atof(argv[idx++]);
  }

  // Kernel
  double s = 1.0; // s is always 1
  double ell = atof(argv[idx++]);
  double nu = atof(argv[idx++]);
  double tau = atof(argv[idx++]);
  double Param[NUM_PARAM] = {ell, nu, tau};

  // Param grid search
  INTEGER Num_ell = String2Integer(argv[idx++]);
  double *List_ell = NULL;
  New_1D_Array<double, INTEGER>(&List_ell, Num_ell);
  for (INTEGER i = 0; i < Num_ell; i++) {
    List_ell[i] = atof(argv[idx++]);
  }
  INTEGER Num_nu = String2Integer(argv[idx++]);
  double *List_nu = NULL;
  New_1D_Array<double, INTEGER>(&List_nu, Num_nu);
  for (INTEGER i = 0; i < Num_nu; i++) {
    List_nu[i] = atof(argv[idx++]);
  }
  INTEGER Num_tau = String2Integer(argv[idx++]);
  double *List_tau = NULL;
  New_1D_Array<double, INTEGER>(&List_tau, Num_tau);
  for (INTEGER i = 0; i < Num_tau; i++) {
    List_tau[i] = atof(argv[idx++]);
  }

  // RNG
  INTEGER sSeed = String2Integer(argv[idx++]);
  unsigned Seed;
  if (sSeed < 0) {
    Seed = (unsigned)time(NULL);
  }
  else {
    Seed = (unsigned)sSeed;
  }

  // Train/test split
  double PortionTrain = atof(argv[idx++]);

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

  // Instantiate kernel
  IsotropicMatern mKernel(s, nu, ell);

  // Simulate a random field y
  DVector y;
  Sampling_Standard<IsotropicMatern, DPoint, DPointArray> mSampler;
  double lambda = pow(10.0, tau);
  START_CLOCK;
  mSampler.GenRandomField(X, y, mKernel, lambda);
  END_CLOCK;
  double TimeSampling = ELAPSED_TIME;
  printf("GP_Standard: Sampling time = %gs\n", TimeSampling); fflush(stdout);

  // Output random field to file
  if (OutputRandomField) {
    WriteRandomFieldToFile(y, d, Dim, Lower, Upper, NUM_PARAM, Param,
                           RandomFieldFileBasename, "GP_Standard");
  }

  // Train/test split
  DPointArray Xtrain, Xtest;
  DVector ytrain, ytest;
  INTEGER *IdxTrain = NULL, *IdxTest = NULL;
  TrainTestSplit(X, y, Xtrain, Xtest, ytrain, ytest, &IdxTrain, &IdxTest,
                 PortionTrain);

  // Save some memory
  X.ReleaseAllMemory();
  y.ReleaseAllMemory();

  // MLE through grid search
  INTEGER ListLength = Num_ell * Num_nu * Num_tau;
  double *mLogLik = NULL;
  New_1D_Array<double, INTEGER>(&mLogLik, ListLength);
  double **ListParam = NULL;
  New_2D_Array<double, INTEGER, INTEGER>(&ListParam, ListLength, NUM_PARAM);
  MLE_Standard<IsotropicMatern, DPoint, DPointArray> mMLE;
  START_CLOCK;
  INTEGER t = 0;
  for (INTEGER i = 0; i < Num_ell; i++) {
    double ell = List_ell[i];
    for (INTEGER j = 0; j < Num_nu; j++) {
      double nu = List_nu[j];
      for (INTEGER k = 0; k < Num_tau; k++) {
        double tau = List_tau[k];
        ListParam[t][0] = ell;
        ListParam[t][1] = nu;
        ListParam[t][2] = tau;
        double lambda   = pow(10.0, tau);
        IsotropicMatern mKernel(s, nu, ell);
        mLogLik[t] = mMLE.LogLik(Xtrain, ytrain, mKernel, lambda);
        printf("MLE_Standard: Grid search ell = %g, nu = %g, tau = %g, loglik = %.16e\n", ell, nu, tau, mLogLik[t]); fflush(stdout);
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
  double hat_ell = HatParam[0];
  double hat_nu  = HatParam[1];
  double hat_tau = HatParam[2];
  printf("GP_Standard: Truth     ell = %g, nu = %g, tau = %g\n", ell, nu, tau); fflush(stdout);
  printf("GP_Standard: Estimated ell = %g, nu = %g, tau = %g, max loglik = %.16e, MLE time = %gs\n", hat_ell, hat_nu, hat_tau, MaxLogLik, TimeMLE); fflush(stdout);

  // Output all loglik's to file
  if (OutputLogLik) {
    WriteLogLikToFile(NUM_PARAM, ListLength, ListParam, mLogLik,
                      LogLikFileName, "GP_Standard");
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
        double ell    = CenterParam[0];
        double nu     = CenterParam[1];
        double tau    = CenterParam[2];
        double lambda = pow(10.0, tau);
        IsotropicMatern mKernel(s, nu, ell);
        double mLogLikp = mMLE.LogLik(Xtrain, ytrain, mKernel, lambda);

        memcpy(CenterParam, HatParam, NUM_PARAM*sizeof(double));
        CenterParam[i] -= epsilon;
        ell    = CenterParam[0];
        nu     = CenterParam[1];
        tau    = CenterParam[2];
        lambda = pow(10.0, tau);
        mKernel.Init(s, nu, ell);
        double mLogLikm = mMLE.LogLik(Xtrain, ytrain, mKernel, lambda);

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
      double ell    = ListParam2[i][0];
      double nu     = ListParam2[i][1];
      double tau    = ListParam2[i][2];
      double lambda = pow(10.0, tau);
      IsotropicMatern mKernel(s, nu, ell);
      mLogLik2[i] = mMLE.LogLik(Xtrain, ytrain, mKernel, lambda);
      printf("MLE_Standard: Fisher info ell = %g, nu = %g, tau = %g, loglik = %.16e\n", ell, nu, tau, mLogLik2[i]); fflush(stdout);
    }

    DMatrix Fisher, Cov;
    DVector Stderr;
    ComputeFisher(NUM_PARAM, ListLength2, ListParam2, mLogLik2,
                  HatParam, MaxLogLik, DiffStepSize, Fisher, Cov, Stderr);
    double *mStderr = Stderr.GetPointer();
    printf("GP_Standard: Stderr ell = %g, nu = %g, tau = %g\n", mStderr[0], mStderr[1], mStderr[2]); fflush(stdout);

    Delete_1D_Array<double>(&mLogLik2);
    Delete_2D_Array<double, INTEGER>(&ListParam2, ListLength2);

    // Output Fisher to file
    if (OutputFisher) {
      WriteFisherToFile(Fisher, FisherFileName, "GP_Standard");
    }

  }

  // Use estimated parameters to do kriging
  Kriging_Standard<IsotropicMatern, DPoint, DPointArray> mKriging;
  DVector ytest_predict, ytest_stddev;
  double hat_lambda = pow(10.0, hat_tau);
  IsotropicMatern mKernel2(s, hat_nu, hat_ell);
  START_CLOCK;
  mKriging.Train(Xtrain, mKernel2, hat_lambda);
  mKriging.Test(Xtrain, Xtest, ytrain, mKernel2, hat_lambda,
                ytest_predict, ytest_stddev);
  END_CLOCK;
  double TimeKriging = ELAPSED_TIME;
  printf("GP_Standard: Kriging time = %gs\n", TimeKriging); fflush(stdout);

  // Output kriged random field to file
  if (OutputKrigedRandomField) {
    DVector y_all;
    AssembleY(ytrain, ytest_predict, IdxTrain, IdxTest, y_all);
    WriteRandomFieldToFile(y_all, d, Dim, Lower, Upper, NUM_PARAM, HatParam,
                           KrigedRandomFieldFileBasename, "GP_Standard");
    y_all.ReleaseAllMemory();
  }

  // Output predictions to file
  if (OutputPredictions) {
    WritePredictionsToFile(ytest, ytest_predict, ytest_stddev,
                           PredictionsFileName, "GP_Standard");
  }

  //---------- Clean up --------------------

  Delete_1D_Array<INTEGER>(&Dim);
  Delete_1D_Array<double>(&Lower);
  Delete_1D_Array<double>(&Upper);
  Delete_1D_Array<double>(&List_ell);
  Delete_1D_Array<double>(&List_nu);
  Delete_1D_Array<double>(&List_tau);
  Delete_1D_Array<INTEGER>(&IdxTrain);
  Delete_1D_Array<INTEGER>(&IdxTest);
  Delete_1D_Array<double>(&mLogLik);
  Delete_2D_Array<double, INTEGER>(&ListParam, ListLength);

  return 0;

}
