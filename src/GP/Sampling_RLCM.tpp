#ifndef _SAMPLING_RLCM_TPP_
#define _SAMPLING_RLCM_TPP_


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
Sampling_RLCM<Kernel, Point, PointArray>::
Sampling_RLCM() {
  Init();
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void Sampling_RLCM<Kernel, Point, PointArray>::
Init(void) {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void Sampling_RLCM<Kernel, Point, PointArray>::
ReleaseAllMemory(void) {
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
Sampling_RLCM<Kernel, Point, PointArray>::
Sampling_RLCM(const Sampling_RLCM &G) {
  Init();
  DeepCopy(G);
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
Sampling_RLCM<Kernel, Point, PointArray>&
Sampling_RLCM<Kernel, Point, PointArray>::
operator= (const Sampling_RLCM &G) {
  if (this != &G) {
    DeepCopy(G);
  }
  return *this;
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void Sampling_RLCM<Kernel, Point, PointArray>::
DeepCopy(const Sampling_RLCM &G) {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
Sampling_RLCM<Kernel, Point, PointArray>::
~Sampling_RLCM() {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void Sampling_RLCM<Kernel, Point, PointArray>::
GenRandomField(CMatrix &K,  // Must have called CMatrix::BuildTree()
               const PointArray &X,  // Must use permuted X
               DVector &y,           // In the permuting order of X
               const Kernel &mKernel,
               double lambda) const {

  // K = phi(X,X) + lambda*I
  K.BuildKernelMatrix<Kernel, Point, PointArray>(mKernel, X, lambda);

  // Cholesky factorization K = GG'
  CMatrix G;
  K.Sqrt(G);

  // z = randn(n,1)
  INTEGER N = X.GetN();
  DVector z(N);
  z.SetStandardNormal();

  // y = G * z
  G.MatVec(z, y, NORMAL);

}


#endif
