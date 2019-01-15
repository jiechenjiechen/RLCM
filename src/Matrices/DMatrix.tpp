#ifndef _DMATRIX_TPP_
#define _DMATRIX_TPP_


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void DMatrix::
BuildKernelMatrix(const Kernel &mKernel, const PointArray &X, double lambda) {

  if (mKernel.IsSymmetric() == false) {
    BuildKernelMatrix<Kernel, Point, PointArray>(mKernel, X, X, lambda);
    return;
  }

  Init(X.GetN(), X.GetN());
  for (INTEGER i = 0; i < M; i++) {
    Point x;
    X.GetPoint(i, x);
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (INTEGER j = 0; j <= i; j++) {
      Point y;
      X.GetPoint(j, y);
      A[j*M+i] = mKernel.Eval(x, y, lambda);
      A[i*M+j] = A[j*M+i];
    }
  }

}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void DMatrix::
BuildKernelMatrix(const Kernel &mKernel, const PointArray &X,
                  const PointArray &Y, double lambda) {

  Init(X.GetN(), Y.GetN());

  if (M > N) {

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (INTEGER i = 0; i < M; i++) {
      Point x;
      X.GetPoint(i, x);
      for (INTEGER j = 0; j < N; j++) {
        Point y;
        Y.GetPoint(j, y);
        A[j*M+i] = mKernel.Eval(x, y, lambda);
      }
    }

  }
  else {

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (INTEGER j = 0; j < N; j++) {
      Point y;
      Y.GetPoint(j, y);
      for (INTEGER i = 0; i < M; i++) {
        Point x;
        X.GetPoint(i, x);
        A[j*M+i] = mKernel.Eval(x, y, lambda);
      }
    }

  }

}


#endif
