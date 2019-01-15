#ifndef _DVECTOR_TPP_
#define _DVECTOR_TPP_


//--------------------------------------------------------------------------
template<class TestFunction, class Point, class PointArray>
void DVector::
BuildResponseVector(const TestFunction &mTestFunction, const PointArray &X) {

  Init(X.GetN());
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < N; i++) {
    Point x;
    X.GetPoint(i, x);
    a[i] = mTestFunction.Eval(x);
  }

}


#endif
