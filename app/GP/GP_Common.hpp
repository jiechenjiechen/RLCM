// This file contains miscellaneous stuffs supporting the programs
// implemented in this directory.
//
// The implementation is parallelized.

#ifndef _GP_COMMON_
#define _GP_COMMON_

#include "LibCMatrix.hpp"

// File format:
// .val contains the y vector, one number per line
// .info contains 5 + NumParam rows:
//    d
//    Dim[0], ... Dim[d-1]
//    Lower[0], ... Lower[d-1]
//    Upper[0], ... Upper[d-1]
//    NumParam
//    Param[0]
//    ...
//    Param[NumParam-1]
void WriteRandomFieldToFile(const DVector &y, INTEGER d, const INTEGER *Dim,
                            const double *Lower, const double *Upper,
                            INTEGER NumParam, const double *Param,
                            const char *FileBasename, const char *PrintString);

void ReadRandomFieldInfo(INTEGER &d, INTEGER **Dim, INTEGER &N,
                         double **Lower, double **Upper,
                         INTEGER &NumParam, double **Param,
                         const char *FileBasename, const char *PrintString);
void ReadRandomField(DVector &y, const char *FileBasename,
                     const char *PrintString);
void ReadTrainTestSplit(INTEGER **IdxTrain, INTEGER **IdxTest,
                        INTEGER &Ntrain, INTEGER &Ntest,
                        const char *FileBasename, const char *PrintString);

// Arrays IdxTrain and IdxTest should be empty before this routine is
// called. This routine will allocate memory for the arrays. The user
// is responsible to free the memory when the arrays are no longer
// used.
void TrainTestSplit(const DPointArray &X, const DVector &y, DPointArray &Xtrain,
                    DPointArray &Xtest, DVector &ytrain, DVector &ytest,
                    INTEGER **IdxTrain, INTEGER **IdxTest,
                    double PortionTrain);

void WriteTrainTestSplitToFile(const DVector &ytrain, const DVector &ytest,
                               const INTEGER *IdxTrain, const INTEGER *IdxTest,
                               const char *FileBasename,
                               const char *PrintString);

// Arrays IdxTrain and IdxTest should be empty before this routine is
// called. This routine will allocate memory for the arrays. The user
// is responsible to free the memory when the arrays are no longer
// used.
//
// x : train    o : test
//
//   x o x o x o x o x o x o x o x o
//   o o o o o o o o o o o o o o o o
//   x o x o x o x o x o x o x o x o
//   o o o o o o o o o o o o o o o o
//   x o x o x o x o x o x o x o x o
//   o o o o o o o o o o o o o o o o
//   x o x o x o x o x o x o x o x o
void TrainTestSplitMultipleSamples(const DPointArray &X, const DVector &y,
                                   const INTEGER *Dim, INTEGER NumSamples,
                                   DPointArray &Xtrain, DPointArray &Xtest,
                                   DMatrix &Ytrain, DVector &ytest,
                                   INTEGER **IdxTrain, INTEGER **IdxTest);

// Find the max loglik
void EstimatedParam(INTEGER NumParam, INTEGER ListLength,
                    const double * const *ListParam, const double *mLogLik,
                    double *HatParam, double &MaxLogLik);

// Compute Fisher information
void PrepareListParamForFisher(INTEGER NumParam, INTEGER ListLength,
                               double **ListParam, const double *Param,
                               const double *DiffStepSize);
void ComputeFisher(INTEGER NumParam, INTEGER ListLength,
                   const double * const *ListParam, double *mLogLik,
                   const double *Param, double MaxLogLik,
                   const double *DiffStepSize,
                   DMatrix &Fisher, DMatrix &Cov, DVector &Stderr);

// File format:
//  This file contains ListLength rows. Each row has NumParam+1 numbers:
//    ...
//    ListParam[i][0]  ... ListPararam[i][NumParam-1]  mLogLik[i]
//    ...
void WriteLogLikToFile(INTEGER NumParam, INTEGER ListLength,
                       const double * const *ListParam, const double *mLogLik,
                       const char *FileName, const char *PrintString);

// File format:
//  This file contains one matrix of size NumParam * NumParam
void WriteFisherToFile(const DMatrix &Fisher,
                       const char *FileName, const char *PrintString);

// Assemble ytrain and ytest to y
void AssembleY(const DVector &ytrain, const DVector &ytest,
               const INTEGER *IdxTrain, const INTEGER *IdxTest, DVector &y);

// File format:
//  This file contains ytest.GetN() rows. Each row contains three numbers:
//    ...
//    ytest[i] ytest_predict[i] ytest_stddev[i]
//    ...
void WritePredictionsToFile(const DVector &ytest,
                            const DVector &ytest_predict,
                            const DVector &ytest_stddev,
                            const char *FileName,
                            const char *PrintString);

// Auxilary function.
// Return value:
//  -1  Fail
//   0  Empty array
//  >0  Normal array
INTEGER ReadArrayFromFile(INTEGER **A, const char *FileName);
INTEGER ReadArrayFromFile(double **A, const char *FileName);


#endif
