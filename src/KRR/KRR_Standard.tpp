#ifndef _KRR_STANDARD_TPP_
#define _KRR_STANDARD_TPP_


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
KRR_Standard<Kernel, Point, PointArray>::
KRR_Standard() {
  Init();
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void KRR_Standard<Kernel, Point, PointArray>::
Init(void) {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void KRR_Standard<Kernel, Point, PointArray>::
ReleaseAllMemory(void) {
  K.ReleaseAllMemory();
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
KRR_Standard<Kernel, Point, PointArray>::
KRR_Standard(const KRR_Standard &G) {
  Init();
  DeepCopy(G);
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
KRR_Standard<Kernel, Point, PointArray>& KRR_Standard<Kernel, Point, PointArray>::
operator= (const KRR_Standard &G) {
  if (this != &G) {
    DeepCopy(G);
  }
  return *this;
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void KRR_Standard<Kernel, Point, PointArray>::
DeepCopy(const KRR_Standard &G) {
  ReleaseAllMemory();
  K = G.K;
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
KRR_Standard<Kernel, Point, PointArray>::
~KRR_Standard() {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
double KRR_Standard<Kernel, Point, PointArray>::
Train(const PointArray &Xtrain,
      const Kernel &mKernel,
      double lambda) {

  // K = phi(X,X)
  K.BuildKernelMatrix<Kernel, Point, PointArray>(mKernel, Xtrain);

  // Add lambda to diagonal
  K.AddDiagonal(lambda);

  // Cholesky factorization
  K.DPOTRF(UPPER);

  // Memory estimation
  double mem_est = (double)Xtrain.GetN();
  return mem_est;

}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void KRR_Standard<Kernel, Point, PointArray>::
Test(const PointArray &Xtrain,
     const PointArray &Xtest,
     const DVector &ytrain,
     const Kernel &mKernel,
     DVector &ytest_predict) const {

  // yy = y-ymean
  double ymean = ytrain.Mean();
  DVector yy;
  ytrain.Subtract(ymean, yy);

  // c = (K+lambda*I)\(y-ymean)
  DVector c;
  K.DPOTRS(yy, c, UPPER); // c = K\yy

  // K0 = phi(Xtest,Xtrain)
  DMatrix K0;
  K0.BuildKernelMatrix<Kernel, Point, PointArray>(mKernel, Xtest, Xtrain);

  // y = K0*c + ymean
  K0.MatVec(c, ytest_predict, NORMAL);
  ytest_predict.Add(ymean);

}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void KRR_Standard<Kernel, Point, PointArray>::
Test(const PointArray &Xtrain,
     const PointArray &Xtest,
     const DMatrix &Ytrain,
     const Kernel &mKernel,
     DMatrix &Ytest_predict) const {

  // yy = y-ymean
  INTEGER m = Ytrain.GetN();
  double *ymean = NULL;
  New_1D_Array<double, INTEGER>(&ymean, m);
  DMatrix YY(Ytrain.GetM(), m);
  for (INTEGER i = 0; i < m; i++) {
    DVector ytrain;
    Ytrain.GetColumn(i, ytrain);
    ymean[i] = ytrain.Mean();
    ytrain.Subtract(ymean[i]);
    YY.SetColumn(i, ytrain);
  }

  // c = (K+lambda*I)\(y-ymean)
  DMatrix C;
  K.DPOTRS(YY, C, UPPER); // c = K\yy

  // K0 = phi(Xtest,Xtrain)
  DMatrix K0;
  K0.BuildKernelMatrix<Kernel, Point, PointArray>(mKernel, Xtest, Xtrain);

  // y = K0*c + ymean
  K0.MatMat(C, Ytest_predict, NORMAL, NORMAL);
  for (INTEGER i = 0; i < m; i++) {
    DVector ytest_predict;
    Ytest_predict.GetColumn(i, ytest_predict);
    ytest_predict.Add(ymean[i]);
    Ytest_predict.SetColumn(i, ytest_predict);
  }

  // Clean up
  Delete_1D_Array<double>(&ymean);

}


#endif
