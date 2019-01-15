#ifndef _MLE_STANDARD_TPP_
#define _MLE_STANDARD_TPP_


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
MLE_Standard<Kernel, Point, PointArray>::
MLE_Standard() {
  Init();
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void MLE_Standard<Kernel, Point, PointArray>::
Init(void) {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void MLE_Standard<Kernel, Point, PointArray>::
ReleaseAllMemory(void) {
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
MLE_Standard<Kernel, Point, PointArray>::
MLE_Standard(const MLE_Standard &G) {
  Init();
  DeepCopy(G);
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
MLE_Standard<Kernel, Point, PointArray>&
MLE_Standard<Kernel, Point, PointArray>::
operator= (const MLE_Standard &G) {
  if (this != &G) {
    DeepCopy(G);
  }
  return *this;
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void MLE_Standard<Kernel, Point, PointArray>::
DeepCopy(const MLE_Standard &G) {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
MLE_Standard<Kernel, Point, PointArray>::
~MLE_Standard() {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
double MLE_Standard<Kernel, Point, PointArray>::
LogLik(const PointArray &X,
       const DVector &y,
       const Kernel &mKernel,
       double lambda) const {

  // K = phi(X,X)
  DMatrix K;
  K.BuildKernelMatrix<Kernel, Point, PointArray>(mKernel, X);

  // Add lambda to diagonal
  K.AddDiagonal(lambda);

  // Cholesky factorization
  K.DPOTRF(LOWER);

  // -1/2*y'*inv(K)*y = -1/2*(invG*y)'*(invG*y), if K = GG'
  DVector invGy;
  K.DTRSV(y, invGy, NORMAL, LOWER);
  double Term1 = -0.5 * invGy.InProd(invGy);

  // -1/2*logdet(K)
  LogDet logdetK = K.Det(SPD, CHOL_FACT);
  if (logdetK.Sign < 0) {
    printf("MLE_Standard::LogLik. Error: Kernel matrix is not positive-definite. Return NAN.\n");
    return NAN;
  }
  double Term2 = -0.5 * logdetK.LogAbsDet;

  // -n/2*log(2*pi)
  INTEGER N = X.GetN();
  double Term3 = -0.5 * N * log(PIx2);

  // Total
  double Total = Term1 + Term2 + Term3;
  return Total;

}


#endif
