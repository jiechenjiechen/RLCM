#include "GP_Common.hpp"


//--------------------------------------------------------------------------
void WriteRandomFieldToFile(const DVector &y, INTEGER d, const INTEGER *Dim,
                            const double *Lower, const double *Upper,
                            INTEGER NumParam, const double *Param,
                            const char *FileBasename, const char *PrintString) {

  char FileName[1000];

  sprintf(FileName, "%s.val", FileBasename);
  FILE *fp = fopen(FileName, "w");
  double *my = y.GetPointer();
  INTEGER N = y.GetN();
  for (INTEGER i = 0; i < N; i++) {
    fprintf(fp, "%.15e\n", my[i]);
  }
  fclose(fp);

  printf("%s: Written random field to file %s\n", PrintString, FileName); fflush(stdout);

  sprintf(FileName, "%s.info", FileBasename);
  fp = fopen(FileName, "w");
  fprintf(fp, "%ld\n", (long)d);
  for (INTEGER i = 0; i < d; i++) {
    fprintf(fp, "%ld ", (long)(Dim[i]));
  }
  fprintf(fp, "\n");
  for (INTEGER i = 0; i < d; i++) {
    fprintf(fp, "%.15e ", Lower[i]);
  }
  fprintf(fp, "\n");
  for (INTEGER i = 0; i < d; i++) {
    fprintf(fp, "%.15e ", Upper[i]);
  }
  fprintf(fp, "\n");
  fprintf(fp, "%ld\n", (long)NumParam);
  for (INTEGER i = 0; i < NumParam; i++) {
    fprintf(fp, "%.15e\n", Param[i]);
  }
  fclose(fp);

}


//--------------------------------------------------------------------------
void ReadRandomFieldInfo(INTEGER &d, INTEGER **Dim, INTEGER &N,
                         double **Lower, double **Upper,
                         INTEGER &NumParam, double **Param,
                         const char *FileBasename, const char *PrintString) {

  char FileName[1000];
  sprintf(FileName, "%s.info", FileBasename);
  FILE *fp = fopen(FileName, "r");

  int some_int = 0;
  fscanf(fp, "%d", &some_int);
  d = (INTEGER)some_int;
  New_1D_Array<INTEGER, INTEGER>(Dim, d);
  New_1D_Array<double, INTEGER>(Lower, d);
  New_1D_Array<double, INTEGER>(Upper, d);

  N = 1;
  for (INTEGER i = 0; i < d; i++) {
    fscanf(fp, "%d", &some_int);
    (*Dim)[i] = (INTEGER)some_int;
    N *= (*Dim)[i];
  }
  for (INTEGER i = 0; i < d; i++) {
    fscanf(fp, "%lf", (*Lower)+i);
  }
  for (INTEGER i = 0; i < d; i++) {
    fscanf(fp, "%lf", (*Upper)+i);
  }

  fscanf(fp, "%d", &some_int);
  NumParam = (INTEGER)some_int;
  New_1D_Array<double, INTEGER>(Param, NumParam);

  for (INTEGER i = 0; i < NumParam; i++) {
    fscanf(fp, "%lf", (*Param)+i);
  }

  fclose(fp);

  printf("%s: Loaded random field info from file %s\n", PrintString, FileName); fflush(stdout);

}


//--------------------------------------------------------------------------
void ReadRandomField(DVector &y, const char *FileBasename,
                     const char *PrintString) {

  char FileName[1000];
  sprintf(FileName, "%s.val", FileBasename);
  double *A = NULL;
  INTEGER N = ReadArrayFromFile(&A, FileName);
  if (N <= 0) {
    return;
  }
  y.Init(N);
  memcpy(y.GetPointer(), A, N*sizeof(double));
  Delete_1D_Array<double>(&A);

  printf("%s: Loaded random field from file %s\n", PrintString, FileName); fflush(stdout);

}


//--------------------------------------------------------------------------
void ReadTrainTestSplit(INTEGER **IdxTrain, INTEGER **IdxTest,
                        INTEGER &Ntrain, INTEGER &Ntest,
                        const char *FileBasename, const char *PrintString) {

  char FileName[1000];
  sprintf(FileName, "%s.idxTrain", FileBasename);
  Ntrain = ReadArrayFromFile(IdxTrain, FileName);
  printf("%s: Loaded training indices from file %s\n", PrintString, FileName); fflush(stdout);

  sprintf(FileName, "%s.idxTest", FileBasename);
  Ntest = ReadArrayFromFile(IdxTest, FileName);
  printf("%s: Loaded testing indices from file %s\n", PrintString, FileName); fflush(stdout);

}


//--------------------------------------------------------------------------
void TrainTestSplit(const DPointArray &X, const DVector &y, DPointArray &Xtrain,
                    DPointArray &Xtest, DVector &ytrain, DVector &ytest,
                    INTEGER **IdxTrain, INTEGER **IdxTest,
                    double PortionTrain) {

  INTEGER N = X.GetN();
  INTEGER Ntrain = (int)(N * PortionTrain);
  INTEGER Ntest = N - Ntrain;

  INTEGER *IdxTmp = NULL;
  New_1D_Array<INTEGER, INTEGER>(IdxTrain, Ntrain);
  New_1D_Array<INTEGER, INTEGER>(IdxTest, Ntest);
  New_1D_Array<INTEGER, INTEGER>(&IdxTmp, N); // Init zero

  RandPerm(N, Ntrain, *IdxTrain);
  for (INTEGER i = 0; i < Ntrain; i++) {
    IdxTmp[(*IdxTrain)[i]] = 1;
  }
  INTEGER j = 0;
  for (INTEGER i = 0; i < N; i++) {
    if (IdxTmp[i] == 0) {
      (*IdxTest)[j++] = i;
    }
  }

  X.GetSubset(*IdxTrain, Ntrain, Xtrain);
  X.GetSubset(*IdxTest, Ntest, Xtest);
  y.GetBlock(*IdxTrain, Ntrain, ytrain);
  y.GetBlock(*IdxTest, Ntest, ytest);

  Delete_1D_Array(&IdxTmp);

}


//--------------------------------------------------------------------------
void WriteTrainTestSplitToFile(const DVector &ytrain, const DVector &ytest,
                               const INTEGER *IdxTrain, const INTEGER *IdxTest,
                               const char *FileBasename,
                               const char *PrintString) {

  char FileName[1000];

  sprintf(FileName, "%s.idxTrain", FileBasename);
  FILE *fp = fopen(FileName, "w");
  INTEGER Ntrain = ytrain.GetN();
  for (INTEGER i = 0; i < Ntrain; i++) {
    fprintf(fp, "%d\n", IdxTrain[i]);
  }
  fclose(fp);

  printf("%s: Written training indices to file %s\n", PrintString, FileName); fflush(stdout);

  sprintf(FileName, "%s.idxTest", FileBasename);
  fp = fopen(FileName, "w");
  INTEGER Ntest = ytest.GetN();
  for (INTEGER i = 0; i < Ntest; i++) {
    fprintf(fp, "%d\n", IdxTest[i]);
  }
  fclose(fp);

  printf("%s: Written testing indices to file %s\n", PrintString, FileName); fflush(stdout);

}


//--------------------------------------------------------------------------
void TrainTestSplitMultipleSamples(const DPointArray &X, const DVector &y,
                                   const INTEGER *Dim, INTEGER NumSamples,
                                   DPointArray &Xtrain, DPointArray &Xtest,
                                   DMatrix &Ytrain, DVector &ytest,
                                   INTEGER **IdxTrain, INTEGER **IdxTest) {

  INTEGER N = X.GetN(); // Should = Dim[0]*Dim[1]
  INTEGER Ntrain = ( Dim[0]/2 + Dim[0]%2 ) * ( Dim[1]/2 ); // Dim[1] is even
  INTEGER Ntest = N - Ntrain;

  New_1D_Array<INTEGER, INTEGER>(IdxTrain, Ntrain);
  New_1D_Array<INTEGER, INTEGER>(IdxTest, Ntest);

  INTEGER cur_train = 0, cur_test = 0;
  for (INTEGER j = 0; j < Dim[1]; j++) {
    bool j_odd = j%2;
    for (INTEGER i = 0; i < Dim[0]; i++) {
      INTEGER idx = i + j*Dim[0];
      if ( !j_odd && !(i%2) ) {
        (*IdxTrain)[cur_train++] = idx;
      }
      else {
        (*IdxTest)[cur_test++] = idx;
      }
    }
  } // cur_train should = Ntrain at the end. cur_test = Ntest

  X.GetSubset(*IdxTrain, Ntrain, Xtrain);
  X.GetSubset(*IdxTest, Ntest, Xtest);

  INTEGER *IdxTrainCopy = NULL, *IdxTestCopy = NULL;
  New_1D_Array<INTEGER, INTEGER>(&IdxTrainCopy, Ntrain);
  New_1D_Array<INTEGER, INTEGER>(&IdxTestCopy, Ntest);
  memcpy(IdxTrainCopy, *IdxTrain, Ntrain*sizeof(INTEGER));
  memcpy(IdxTestCopy, *IdxTest, Ntest*sizeof(INTEGER));
  DVector y2;
  INTEGER Offset = Dim[0] * Dim[1];
  Ytrain.Init(Ntrain, NumSamples);
  ytest.Init(Ntest);
  for (INTEGER i = 0; i < NumSamples; i++) {
    y.GetBlock(IdxTrainCopy, Ntrain, y2);
    Ytrain.SetColumn(i, y2);
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (INTEGER j = 0; j < Ntrain; j++) {
      IdxTrainCopy[j] += Offset;
    }
    if (i == 0) {
      y.GetBlock(IdxTestCopy, Ntest, ytest);
    }
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (INTEGER j = 0; j < Ntest; j++) {
      IdxTestCopy[j] += Offset;
    }
  }
  Delete_1D_Array<INTEGER>(&IdxTrainCopy);
  Delete_1D_Array<INTEGER>(&IdxTestCopy);

}


//--------------------------------------------------------------------------
void EstimatedParam(INTEGER NumParam, INTEGER ListLength,
                    const double * const *ListParam, const double *mLogLik,
                    double *HatParam, double &MaxLogLik) {

  MaxLogLik = -DBL_MAX;
  for (INTEGER i = 0; i < NumParam; i++) {
    HatParam[i] = ListParam[0][i];
  }
  for (INTEGER i = 0; i < ListLength; i++) {
    if (MaxLogLik < mLogLik[i]) {
      MaxLogLik = mLogLik[i];
      for (INTEGER j = 0; j < NumParam; j++) {
        HatParam[j] = ListParam[i][j];
      }
    }
  }

}


//--------------------------------------------------------------------------
void PrepareListParamForFisher(INTEGER NumParam, INTEGER ListLength,
                               double **ListParam, const double *Param,
                               const double *DiffStepSize) {

  // For example, if NumParam = 3,
  //  Param1:  + - | 0 0 | 0 0 | + - + - | + - + - | 0 0 0 0
  //  Param2:  0 0 | + - | 0 0 | + - - + | 0 0 0 0 | + - + -
  //  Param3:  0 0 | 0 0 | + - | 0 0 0 0 | + - - + | + - - +

  for (INTEGER i = 0; i < ListLength; i++) {
    for (INTEGER j = 0; j < NumParam; j++) {
      ListParam[i][j] = Param[j];
    }
  }

  for (INTEGER i = 0; i < NumParam; i++) {
    ListParam[2*i][i]   += DiffStepSize[i];
    ListParam[2*i+1][i] -= DiffStepSize[i];
  }

  INTEGER off = 2*NumParam;
  for (INTEGER i = 0; i < NumParam; i++) {
    for (INTEGER j = i+1; j < NumParam; j++) {
      ListParam[off][i]   += DiffStepSize[i];
      ListParam[off][j]   += DiffStepSize[j];
      ListParam[off+1][i] -= DiffStepSize[i];
      ListParam[off+1][j] -= DiffStepSize[j];
      ListParam[off+2][i] += DiffStepSize[i];
      ListParam[off+2][j] -= DiffStepSize[j];
      ListParam[off+3][i] -= DiffStepSize[i];
      ListParam[off+3][j] += DiffStepSize[j];
      off += 4;
    }
  }

}


//--------------------------------------------------------------------------
void ComputeFisher(INTEGER NumParam, INTEGER ListLength,
                   const double * const *ListParam, double *mLogLik,
                   const double *Param, double MaxLogLik,
                   const double *DiffStepSize,
                   DMatrix &Fisher, DMatrix &Cov, DVector &Stderr) {

  // Shorter variable name
  double *b = mLogLik;
  double a = MaxLogLik;

  // Fisher matrix
  Fisher.Init(NumParam);
  double *mFisher = Fisher.GetPointer();
  INTEGER off = 2*NumParam;
  for (INTEGER i = 0; i < NumParam; i++) {
    mFisher[i+i*NumParam] = ( b[2*i] + b[2*i+1] - 2*a ) /
      ( DiffStepSize[i] * DiffStepSize[i] );
    for (INTEGER j = i+1; j < NumParam; j++) {
      mFisher[i+j*NumParam] = ( (b[off] + b[off+1]) - (b[off+2] + b[off+3]) ) /
        ( 4.0 * DiffStepSize[i] * DiffStepSize[j] );
      mFisher[j+i*NumParam] = mFisher[i+j*NumParam];
      off += 4;
    }
  }
  Fisher.Negate(); // Negative Hessian

  Fisher.PrintMatrixMatlabForm("Fisher");

  // Variance matrix. The Fisher matrix should be positive definite;
  // otherwise there is something wrong.
  DMatrix Eye(NumParam);
  Eye.SetIdentity();
  Fisher.Mldivide(Eye, Cov, NORMAL, SPD);

  Cov.PrintMatrixMatlabForm("Cov");

  // Standard error
  Cov.Diag(Stderr);
  Stderr.Sqrt();

  Stderr.PrintVectorMatlabForm("Stderr");

}


//--------------------------------------------------------------------------
void WriteLogLikToFile(INTEGER NumParam, INTEGER ListLength,
                       const double * const *ListParam, const double *mLogLik,
                       const char *FileName, const char *PrintString) {

  FILE *fp = fopen(FileName, "w");
  for (INTEGER i = 0; i < ListLength; i++) {
    for (INTEGER j = 0; j < NumParam; j++) {
      fprintf(fp, "%g ", ListParam[i][j]);
    }
    fprintf(fp, "%.16e", mLogLik[i]);
    fprintf(fp, "\n");
  }
  fclose(fp);

  printf("%s: Written loglik to file %s\n", PrintString, FileName); fflush(stdout);

}


//--------------------------------------------------------------------------
void WriteFisherToFile(const DMatrix &Fisher,
                       const char *FileName, const char *PrintString) {

  FILE *fp = fopen(FileName, "w");
  INTEGER NumParam = Fisher.GetN();
  double *mFisher = Fisher.GetPointer();
  INTEGER t = 0;
  for (INTEGER i = 0; i < NumParam; i++) {
    for (INTEGER j = 0; j < NumParam; j++) {
      fprintf(fp, "%g ", mFisher[t]);
      t++;
    }
    fprintf(fp, "\n");
  }
  fclose(fp);

  printf("%s: Written Fisher to file %s\n", PrintString, FileName); fflush(stdout);

}


//--------------------------------------------------------------------------
void AssembleY(const DVector &ytrain, const DVector &ytest,
               const INTEGER *IdxTrain, const INTEGER *IdxTest, DVector &y) {

  INTEGER Ntrain = ytrain.GetN();
  INTEGER Ntest = ytest.GetN();
  INTEGER N = Ntrain + Ntest;
  y.Init(N);

  double *my = y.GetPointer();
  double *mytrain = ytrain.GetPointer();
  double *mytest = ytest.GetPointer();

  for (INTEGER i = 0; i < Ntrain; i++) {
    my[IdxTrain[i]] = mytrain[i];
  }
  for (INTEGER i = 0; i < Ntest; i++) {
    my[IdxTest[i]] = mytest[i];
  }

}


//--------------------------------------------------------------------------
void WritePredictionsToFile(const DVector &ytest,
                            const DVector &ytest_predict,
                            const DVector &ytest_stddev,
                            const char *FileName,
                            const char *PrintString) {

  FILE *fp = fopen(FileName, "w");
  double *mytest = ytest.GetPointer();
  double *mytest_predict = ytest_predict.GetPointer();
  double *mytest_stddev = ytest_stddev.GetPointer();
  INTEGER N = ytest.GetN();
  for (INTEGER i = 0; i < N; i++) {
    fprintf(fp, "%g %g %g\n", mytest[i], mytest_predict[i], mytest_stddev[i]);
  }
  fclose(fp);

  printf("%s: Written predictions to file %s\n", PrintString, FileName); fflush(stdout);

}


//--------------------------------------------------------------------------
INTEGER ReadArrayFromFile(INTEGER **A, const char *FileName) {

  FILE *fp = fopen(FileName, "r");

  // Read once and get count
  INTEGER N = 0;
  int some_int = 0, ret = 0;
  while ((ret = fscanf(fp, "%d", &some_int)) == 1) {
    N++;
  }
  fclose(fp);

  if (N <= 0) {
    return N;
  }
  else {
    New_1D_Array<INTEGER, INTEGER>(A, N);
  }

  // Actual read
  fp = fopen(FileName, "r");
  for (INTEGER i = 0; i < N; i++) {
    fscanf(fp, "%d", &some_int);
    (*A)[i] = (INTEGER)some_int;
  }
  fclose(fp);

  return N;

}


//--------------------------------------------------------------------------
INTEGER ReadArrayFromFile(double **A, const char *FileName) {

  FILE *fp = fopen(FileName, "r");

  // Read once and get count
  INTEGER N = 0;
  int ret = 0;
  double some_double;
  while ((ret = fscanf(fp, "%lf", &some_double)) == 1) {
    N++;
  }
  fclose(fp);

  if (N <= 0) {
    return N;
  }
  else {
    New_1D_Array<double, INTEGER>(A, N);
  }

  // Actual read
  fp = fopen(FileName, "r");
  for (INTEGER i = 0; i < N; i++) {
    fscanf(fp, "%lf", (*A)+i);
  }
  fclose(fp);

  return N;

}
