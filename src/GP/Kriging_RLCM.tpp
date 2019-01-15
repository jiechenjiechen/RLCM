#ifndef _KRIGING_RLCM_TPP_
#define _KRIGING_RLCM_TPP_


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
Kriging_RLCM<Kernel, Point, PointArray>::
Kriging_RLCM() {
  Init();
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void Kriging_RLCM<Kernel, Point, PointArray>::
Init(void) {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void Kriging_RLCM<Kernel, Point, PointArray>::
ReleaseAllMemory(void) {
  invK.ReleaseAllMemory();
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
Kriging_RLCM<Kernel, Point, PointArray>::
Kriging_RLCM(const Kriging_RLCM &G) {
  Init();
  DeepCopy(G);
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
Kriging_RLCM<Kernel, Point, PointArray>&
Kriging_RLCM<Kernel, Point, PointArray>::
operator= (const Kriging_RLCM &G) {
  if (this != &G) {
    DeepCopy(G);
  }
  return *this;
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void Kriging_RLCM<Kernel, Point, PointArray>::
DeepCopy(const Kriging_RLCM &G) {
  ReleaseAllMemory();
  invK = G.invK;
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
Kriging_RLCM<Kernel, Point, PointArray>::
~Kriging_RLCM() {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void Kriging_RLCM<Kernel, Point, PointArray>::
Train(CMatrix &K,  // Must have called CMatrix::BuildTree()
      const PointArray &Xtrain, // Must use permuted Xtrain
      const Kernel &mKernel,
      double lambda) {

  // K = phi(X,X) + lambda*I
  K.BuildKernelMatrix<Kernel, Point, PointArray>(mKernel, Xtrain, lambda);

  // invK
  K.Invert(invK);

}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void Kriging_RLCM<Kernel, Point, PointArray>::
Test(CMatrix &K,  // Must be the same as that input to Train()
     const PointArray &Xtrain, // Must use permuted Xtrain
     const PointArray &Xtest,
     const DVector &ytrain,    // Must use permuted ytrain
     const Kernel &mKernel,
     double lambda,
     DVector &ytest_predict,
     DVector &ytest_stddev) {

  // w = K\y
  DVector w;
  invK.MatVec(ytrain, w, NORMAL);

  // y0 = phi(Xtest,Xtrain) * w
  K.MatVecImplicit<Kernel, Point, PointArray>
    (Xtrain, Xtest, mKernel, w, ytest_predict);

  // stddev(y0) = sqrt( K(0)+lambda - K0 * inv(K+lambda*I) * K0 )
  DVector z;
  K.BilinearImplicit<Kernel, Point, PointArray>
    (Xtrain, Xtest, mKernel, invK, z);
  z.Negate();
  double s = mKernel.GetS();
  z.Add(s + lambda);
  z.Sqrt(ytest_stddev);

}


#endif
