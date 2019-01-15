#ifndef _SAMPLING_STANDARD_TPP_
#define _SAMPLING_STANDARD_TPP_


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
Sampling_Standard<Kernel, Point, PointArray>::
Sampling_Standard() {
  Init();
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void Sampling_Standard<Kernel, Point, PointArray>::
Init(void) {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void Sampling_Standard<Kernel, Point, PointArray>::
ReleaseAllMemory(void) {
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
Sampling_Standard<Kernel, Point, PointArray>::
Sampling_Standard(const Sampling_Standard &G) {
  Init();
  DeepCopy(G);
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
Sampling_Standard<Kernel, Point, PointArray>&
Sampling_Standard<Kernel, Point, PointArray>::
operator= (const Sampling_Standard &G) {
  if (this != &G) {
    DeepCopy(G);
  }
  return *this;
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void Sampling_Standard<Kernel, Point, PointArray>::
DeepCopy(const Sampling_Standard &G) {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
Sampling_Standard<Kernel, Point, PointArray>::
~Sampling_Standard() {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void Sampling_Standard<Kernel, Point, PointArray>::
GenRandomField(const PointArray &X,
               DVector &y,
               const Kernel &mKernel,
               double lambda) const {

  // K = phi(X,X)
  DMatrix K;
  K.BuildKernelMatrix<Kernel, Point, PointArray>(mKernel, X);

  // Add lambda to diagonal
  K.AddDiagonal(lambda);

  // Cholesky factorization K = GG'
  DMatrix G;
  K.Chol(G, LOWER);

  // z = randn(n,1)
  INTEGER N = X.GetN();
  DVector z(N);
  z.SetStandardNormal();

  // y = G * z
  G.MatVec(z, y, NORMAL);

}


#endif
