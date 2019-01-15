#ifndef _KRIGING_STANDARD_TPP_
#define _KRIGING_STANDARD_TPP_


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
Kriging_Standard<Kernel, Point, PointArray>::
Kriging_Standard() {
  Init();
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void Kriging_Standard<Kernel, Point, PointArray>::
Init(void) {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void Kriging_Standard<Kernel, Point, PointArray>::
ReleaseAllMemory(void) {
  K.ReleaseAllMemory();
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
Kriging_Standard<Kernel, Point, PointArray>::
Kriging_Standard(const Kriging_Standard &G) {
  Init();
  DeepCopy(G);
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
Kriging_Standard<Kernel, Point, PointArray>&
Kriging_Standard<Kernel, Point, PointArray>::
operator= (const Kriging_Standard &G) {
  if (this != &G) {
    DeepCopy(G);
  }
  return *this;
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void Kriging_Standard<Kernel, Point, PointArray>::
DeepCopy(const Kriging_Standard &G) {
  ReleaseAllMemory();
  K = G.K;
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
Kriging_Standard<Kernel, Point, PointArray>::
~Kriging_Standard() {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void Kriging_Standard<Kernel, Point, PointArray>::
Train(const PointArray &Xtrain,
      const Kernel &mKernel,
      double lambda) {

  // K = phi(X,X)
  K.BuildKernelMatrix<Kernel, Point, PointArray>(mKernel, Xtrain);

  // Add lambda to diagonal
  K.AddDiagonal(lambda);

  // Cholesky factorization
  K.DPOTRF(UPPER);

}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void Kriging_Standard<Kernel, Point, PointArray>::
Test(const PointArray &Xtrain,
     const PointArray &Xtest,
     const DVector &ytrain,
     const Kernel &mKernel,
     double lambda,
     DVector &ytest_predict,
     DVector &ytest_stddev) const {

  // K0 = phi(Xtrain,Xtest)
  DMatrix K0;
  K0.BuildKernelMatrix<Kernel, Point, PointArray>(mKernel, Xtrain, Xtest);

  // C = K\K0
  DMatrix C;
  K.DPOTRS(K0, C, UPPER);

  // y0 = C'*y
  C.MatVec(ytrain, ytest_predict, TRANSPOSE);

  // stddev(y0) = sqrt( K(0)+lambda - K0 * inv(K+lambda*I) * K0 )
  DVector c, k0;
  INTEGER m = C.GetM();
  double s = mKernel.GetS();
  ytest_stddev.Init(m);
  double *my = ytest_stddev.GetPointer();
  for (INTEGER i = 0; i < m; i++) {
    C.GetColumn(i, c);
    K0.GetColumn(i, k0);
    my[i] = sqrt( s + lambda - c.InProd(k0) );
  }

}


#endif
