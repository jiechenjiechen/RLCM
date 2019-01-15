#ifndef _MLE_RLCM_TPP_
#define _MLE_RLCM_TPP_


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
MLE_RLCM<Kernel, Point, PointArray>::
MLE_RLCM() {
  Init();
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void MLE_RLCM<Kernel, Point, PointArray>::
Init(void) {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void MLE_RLCM<Kernel, Point, PointArray>::
ReleaseAllMemory(void) {
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
MLE_RLCM<Kernel, Point, PointArray>::
MLE_RLCM(const MLE_RLCM &G) {
  Init();
  DeepCopy(G);
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
MLE_RLCM<Kernel, Point, PointArray>& MLE_RLCM<Kernel, Point, PointArray>::
operator= (const MLE_RLCM &G) {
  if (this != &G) {
    DeepCopy(G);
  }
  return *this;
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void MLE_RLCM<Kernel, Point, PointArray>::
DeepCopy(const MLE_RLCM &G) {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
MLE_RLCM<Kernel, Point, PointArray>::
~MLE_RLCM() {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
double MLE_RLCM<Kernel, Point, PointArray>::
LogLik(CMatrix &K,           // Must have called CMatrix::BuildTree()
       const PointArray &X,  // Must use permuted X
       const DVector &y,     // Must use permuted y
       const Kernel &mKernel,
       double lambda) const {

  // K = phi(X,X) + lambda*I
  K.BuildKernelMatrix<Kernel, Point, PointArray>(mKernel, X, lambda);

  // invK
  CMatrix invK;
  K.Invert(invK);

  // -1/2*y'*inv(K)*y
  DVector invKy;
  invK.MatVec(y, invKy, NORMAL);
  GMRES mGMRES; // Refinement
  mGMRES.RefineSolve<CMatrix, CMatrix>(K, UNFACT, y, invKy, invK, UNFACT);
  double Term1 = -0.5 * y.InProd(invKy);

  // -1/2*logdet(K)
  LogDet logdetK = K.Det(); // Note that must call K.Invert() eariler
  if (logdetK.Sign < 0) {
    printf("MLE_RLCM::LogLik. Error: Kernel matrix is not positive-definite. Return NAN.\n");
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


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
double MLE_RLCM<Kernel, Point, PointArray>::
LogLik(CMatrix &K,           // Must have called CMatrix::BuildTree()
       const PointArray &X,  // Must use permuted X
       const DMatrix &Y,     // Must use permuted Y
       const Kernel &mKernel,
       double lambda) const {

  // K = phi(X,X) + lambda*I
  K.BuildKernelMatrix<Kernel, Point, PointArray>(mKernel, X, lambda);

  // invK
  CMatrix invK;
  K.Invert(invK);

  // Number of samples
  INTEGER NumSamples = Y.GetN();

  // -1/2 sum_i y_i'*inv(K)*y_i
  double Term1 = 0;
  for (INTEGER i = 0; i < NumSamples; i++) {
    DVector y, invKy;
    Y.GetColumn(i, y);
    invK.MatVec(y, invKy, NORMAL);
    GMRES mGMRES; // Refinement
    mGMRES.RefineSolve<CMatrix, CMatrix>(K, UNFACT, y, invKy, invK, UNFACT);
    Term1 += ( -0.5 * y.InProd(invKy) );
  }

  // -1/2*logdet(K) * NumSamples
  LogDet logdetK = K.Det(); // Note that must call K.Invert() eariler
  if (logdetK.Sign < 0) {
    printf("MLE_RLCM::LogLik. Error: Kernel matrix is not positive-definite. Return NAN.\n");
    return NAN;
  }
  double Term2 = -0.5 * logdetK.LogAbsDet * NumSamples;

  // -n/2*log(2*pi) * NumSamples
  INTEGER N = X.GetN();
  double Term3 = -0.5 * N * log(PIx2) * NumSamples;

  // Total
  double Total = Term1 + Term2 + Term3;
  return Total;

}


#endif
