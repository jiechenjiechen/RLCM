#ifndef _KRR_FOURIER_TPP_
#define _KRR_FOURIER_TPP_


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
KRR_Fourier<Kernel, Point, PointArray>::
KRR_Fourier() {
  Init();
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void KRR_Fourier<Kernel, Point, PointArray>::
Init(void) {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void KRR_Fourier<Kernel, Point, PointArray>::
ReleaseAllMemory(void) {
  w.ReleaseAllMemory();
  b.ReleaseAllMemory();
  Z.ReleaseAllMemory();
  ZZ.ReleaseAllMemory();
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
KRR_Fourier<Kernel, Point, PointArray>::
KRR_Fourier(const KRR_Fourier &G) {
  Init();
  DeepCopy(G);
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
KRR_Fourier<Kernel, Point, PointArray>& KRR_Fourier<Kernel, Point, PointArray>::
operator= (const KRR_Fourier &G) {
  if (this != &G) {
    DeepCopy(G);
  }
  return *this;
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void KRR_Fourier<Kernel, Point, PointArray>::
DeepCopy(const KRR_Fourier &G) {
  ReleaseAllMemory();
  w  = G.w;
  b  = G.b;
  Z  = G.Z;
  ZZ = G.ZZ;
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
KRR_Fourier<Kernel, Point, PointArray>::
~KRR_Fourier() {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
double KRR_Fourier<Kernel, Point, PointArray>::
Train(const PointArray &Xtrain,
      const Kernel &mKernel,
      double lambda,
      INTEGER r,
      unsigned Seed) {

  // Seed the RNG
  srandom(Seed);

  // Dimensions
  INTEGER N = Xtrain.GetN();
  INTEGER d = Xtrain.GetD();

  // Generate w. Currently support only the following kernels:
  //   IsotropicGaussian, IsotropicLaplace, ProductLaplace
  w.Init(d, r);
  if (mKernel.GetKernelName() == "IsotropicGaussian") {

    // Standard normal scaled by 1/sigma
    w.SetStandardNormal();
    w.Divide(mKernel.GetSigma());

  }
  else if (mKernel.GetKernelName() == "IsotropicLaplace") {

    // Multivariate student-t of degree 1, scaled by 1/sigma
    w.SetMultivariateStudentT1();
    w.Divide(mKernel.GetSigma());

  }
  else if (mKernel.GetKernelName() == "ProductLaplace") {

    // Student-t of degree 1, scaled by 1/sigma
    w.SetStudentT1();
    w.Divide(mKernel.GetSigma());

  }
  else {

    printf("KRR_Fourier::Train. Error: The supplied kernel is not supported by the Fourier method. Function call takes no effect.\n");
    return NAN;

  }

  // b = rand(1,r)*2*pi
  b.Init(r);
  b.SetUniformRandom01();
  b.Multiply(2.0*M_PI);

  // Z = cos(X*w+repmat(b,n,1)) * sqrt(2/r) * sqrt(s)
  DMatrix Xw;
  Xtrain.MatMat(w, Xw, NORMAL, NORMAL);
  DVector ones(N);
  ones.SetConstVal(1.0);
  DMatrix oneb;
  oneb.OuterProduct(ones, b);
  Xw.Add(oneb, Z);
  Z.Cos();
  Z.Multiply(sqrt(2.0/r*mKernel.GetS()));

  // ZZ = Z'*Z+lambda*I and factorize
  Z.MatMat(Z, ZZ, TRANSPOSE, NORMAL);
  ZZ.AddDiagonal(lambda);
  ZZ.DPOTRF(UPPER);

  // Memory estimation
  double mem_est = (double)r;
  return mem_est;

}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void KRR_Fourier<Kernel, Point, PointArray>::
Test(const PointArray &Xtest,
     const DVector &ytrain,
     const Kernel &mKernel,
     DVector &ytest_predict) const {

  // Dimensions
  INTEGER n = Xtest.GetN();
  INTEGER r = w.GetN();

  // Z0 = cos(X*w+repmat(b,n,1)) * sqrt(2/r) * sqrt(s)
  DMatrix Xw;
  Xtest.MatMat(w, Xw, NORMAL, NORMAL);
  DVector ones(n);
  ones.SetConstVal(1.0);
  DMatrix oneb;
  oneb.OuterProduct(ones, b);
  DMatrix Z0;
  Xw.Add(oneb, Z0);
  Z0.Cos();
  Z0.Multiply(sqrt(2.0/r*mKernel.GetS()));

  // yy = y-ymean
  double ymean = ytrain.Mean();
  DVector yy;
  ytrain.Subtract(ymean, yy);

  // c = (Z'*Z+lambda*I)\(Z'*(y-ymean))
  DVector z2;
  Z.MatVec(yy, z2, TRANSPOSE); // Z2 = Z'*(y-ymean)
  DVector c;
  ZZ.DPOTRS(z2, c, UPPER); // c = ZZ\Z2

  // y = Z0*c + ymean
  Z0.MatVec(c, ytest_predict, NORMAL);
  ytest_predict.Add(ymean);

}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void KRR_Fourier<Kernel, Point, PointArray>::
Test(const PointArray &Xtest,
     const DMatrix &Ytrain,
     const Kernel &mKernel,
     DMatrix &Ytest_predict) const {

  // Dimensions
  INTEGER n = Xtest.GetN();
  INTEGER r = w.GetN();

  // Z0 = cos(X*w+repmat(b,n,1)) * sqrt(2/r) * sqrt(s)
  DMatrix Xw;
  Xtest.MatMat(w, Xw, NORMAL, NORMAL);
  DVector ones(n);
  ones.SetConstVal(1.0);
  DMatrix oneb;
  oneb.OuterProduct(ones, b);
  DMatrix Z0;
  Xw.Add(oneb, Z0);
  Z0.Cos();
  Z0.Multiply(sqrt(2.0/r*mKernel.GetS()));

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

  // c = (Z'*Z+lambda*I)\(Z'*(y-ymean))
  DMatrix Z2;
  Z.MatMat(YY, Z2, TRANSPOSE, NORMAL); // Z2 = Z'*(y-ymean)
  DMatrix C;
  ZZ.DPOTRS(Z2, C, UPPER); // c = ZZ\Z2

  // y = Z0*c + ymean
  Z0.MatMat(C, Ytest_predict, NORMAL, NORMAL);
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
