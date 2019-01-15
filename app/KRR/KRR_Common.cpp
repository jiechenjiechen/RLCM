#include "KRR_Common.hpp"


//--------------------------------------------------------------------------
void ConvertYtrain(const DVector &ytrain, DMatrix &Ytrain, INTEGER NumClasses) {
  Ytrain.Init(ytrain.GetN(), NumClasses);
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < NumClasses; i++) {
    DVector y = ytrain;
    double *my = y.GetPointer();
    for (INTEGER j = 0; j < y.GetN(); j++) {
      my[j] = ((INTEGER)my[j])==i ? 1.0 : -1.0;
    }
    Ytrain.SetColumn(i, y);
  }
}


//--------------------------------------------------------------------------
double Performance(const DVector &ytest_truth, const DVector &ytest_predict,
                   INTEGER NumClasses) {

  INTEGER n = ytest_truth.GetN();
  if (n != ytest_predict.GetN()) {
    printf("Performance. Error: Vector lengths mismatch. Return NAN");
    return NAN;
  }
  if (NumClasses != 1 && NumClasses != 2) {
    printf("Performance. Error: Neither regression nor binary classification. Return NAN");
    return NAN;
  }
  double perf = 0.0;

  if (NumClasses == 1) { // Relative error
    DVector diff;
    ytest_truth.Subtract(ytest_predict, diff);
    perf = diff.Norm2()/ytest_truth.Norm2();
  }
  else if (NumClasses == 2) { // Accuracy
    double *y1 = ytest_truth.GetPointer();
    double *y2 = ytest_predict.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for reduction(+:perf)
#endif
    for (INTEGER i = 0; i < n; i++) {
      perf += y1[i]*y2[i]>0 ? 1.0:0.0;
    }
    perf = perf/n * 100.0;
  }

  return perf;
}


//--------------------------------------------------------------------------
double Performance(const DVector &ytest_truth, const DMatrix &Ytest_predict,
                   INTEGER NumClasses) {

  INTEGER n = ytest_truth.GetN();
  if (n != Ytest_predict.GetM() || Ytest_predict.GetN() != NumClasses ) {
    printf("Performance. Error: Size mismatch. Return NAN");
    return NAN;
  }
  if (NumClasses <= 2) {
    printf("Performance. Error: Not multiclass classification. Return NAN");
    return NAN;
  }
  double perf = 0.0;

  // Compute ytest_predict
  DVector ytest_predict(n);
  double *y = ytest_predict.GetPointer();
  for (INTEGER i = 0; i < n; i++) {
    DVector row(NumClasses);
    Ytest_predict.GetRow(i, row); // GetRow is already parallelized
    INTEGER idx = -1;
    row.Max(idx);
    y[i] = (double)idx;
  }

  // Accuracy
  double *y1 = ytest_truth.GetPointer();
  double *y2 = ytest_predict.GetPointer();
#ifdef USE_OPENMP
#pragma omp parallel for reduction(+:perf)
#endif
  for (INTEGER i = 0; i < n; i++) {
    perf += ((INTEGER)y1[i])==((INTEGER)y2[i]) ? 1.0 : 0.0;
  }
  perf = perf/n * 100.0;

  return perf;
}


//--------------------------------------------------------------------------
void MeanAndStddev(double *a, INTEGER n, double &mean, double &stddev) {

  double m_mean = 0.0;
#ifdef USE_OPENMP
#pragma omp parallel for reduction(+:m_mean)
#endif
  for (INTEGER i = 0; i < n; i++) {
    m_mean += a[i];
  }
  m_mean /= n;

  double m_stddev = 0.0;
#ifdef USE_OPENMP
#pragma omp parallel for reduction(+:m_stddev)
#endif
  for (INTEGER i = 0; i < n; i++) {
    m_stddev += Square(a[i]-m_mean);
  }
  m_stddev = sqrt(m_stddev/n);
  //m_stddev = sqrt(m_stddev/(n-1));

  mean = m_mean;
  stddev = m_stddev;
}
