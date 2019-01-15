// This file contains miscellaneous stuffs supporting the programs
// implemented in this directory.
//
// The implementation is parallelized.

#ifndef _KRR_COMMON_
#define _KRR_COMMON_

#include "LibCMatrix.hpp"

// For multiclass classification, need to convert a single vector
// ytest to a matrix Ytest.
void ConvertYtrain(const DVector &ytrain, DMatrix &Ytrain, INTEGER NumClasses);

// If NumClasses = 1 (regression), return relative error. If
// NumClasses = 2 (binary classification), return accuracy% (between 0
// and 100).
double Performance(const DVector &ytest_truth, const DVector &ytest_predict,
                   INTEGER NumClasses);

// If Numclasses > 2 (multiclass classification), must call this
// routine. The input Ytest_predict consists of NumClasses
// columns. Return accuracy% (between 0 and 100).
double Performance(const DVector &ytest_truth, const DMatrix &Ytest_predict,
                   INTEGER NumClasses);

// Mean and standard deviation of an array a
void MeanAndStddev(double *a, INTEGER n, double &mean, double &stddev);

#endif
