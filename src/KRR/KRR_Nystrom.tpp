#ifndef _KRR_NYSTROM_TPP_
#define _KRR_NYSTROM_TPP_


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
KRR_Nystrom<Kernel, Point, PointArray>::
KRR_Nystrom() {
  Init();
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void KRR_Nystrom<Kernel, Point, PointArray>::
Init(void) {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void KRR_Nystrom<Kernel, Point, PointArray>::
ReleaseAllMemory(void) {
  Y.ReleaseAllMemory();
  K3.ReleaseAllMemory();
  Z.ReleaseAllMemory();
  ZZ.ReleaseAllMemory();
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
KRR_Nystrom<Kernel, Point, PointArray>::
KRR_Nystrom(const KRR_Nystrom &G) {
  Init();
  DeepCopy(G);
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
KRR_Nystrom<Kernel, Point, PointArray>& KRR_Nystrom<Kernel, Point, PointArray>::
operator= (const KRR_Nystrom &G) {
  if (this != &G) {
    DeepCopy(G);
  }
  return *this;
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void KRR_Nystrom<Kernel, Point, PointArray>::
DeepCopy(const KRR_Nystrom &G) {
  ReleaseAllMemory();
  Y  = G.Y;
  K3 = G.K3;
  Z  = G.Z;
  ZZ = G.ZZ;
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
KRR_Nystrom<Kernel, Point, PointArray>::
~KRR_Nystrom() {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
double KRR_Nystrom<Kernel, Point, PointArray>::
Train(const PointArray &Xtrain,
      const Kernel &mKernel,
      double lambda,
      INTEGER r,
      unsigned Seed) {

  // Seed the RNG
  srandom(Seed);

  // Centers Y: random subset of X
  INTEGER *idx = NULL;
  New_1D_Array<INTEGER, INTEGER>(&idx, r);
  RandPerm(Xtrain.GetN(), r, idx);
  Xtrain.GetSubset(idx, r, Y);

  // K1 = phi(X,Y)
  DMatrix K1;
  K1.BuildKernelMatrix<Kernel, Point, PointArray>(mKernel, Xtrain, Y);

  // K2 = phi(Y,Y)
  DMatrix K2;
  K2.BuildKernelMatrix<Kernel, Point, PointArray>(mKernel, Y);

  // K3 = some form of sqrt of pinv(K2)
  DVector D;
  DMatrix V;
  K2.SymEig(D, V); // K2 = V*diag(D)*V'
  INTEGER *idx2 = NULL;
  New_1D_Array<INTEGER, INTEGER>(&idx2, r);
  INTEGER num = D.FindLargerThan(D.Max()*EPS, idx2);
  DMatrix V2;
  V.GetColumns(idx2, num, V2); // V2 = V(:,idx2)
  DVector D2;
  D.GetBlock(idx2, num, D2);
  D2.Sqrt();
  D2.Inv(); // D2 = 1/sqrt(D(idx2))
  DMatrix D3;
  D3.MakeDiag(D2); // D3 = diag(D2)
  V2.MatMat(D3, K3, NORMAL, NORMAL); // K3 = V2*D3

  // Z = K1*K3
  K1.MatMat(K3, Z, NORMAL, NORMAL);

  // ZZ = Z'*Z+lambda*I and factorize
  Z.MatMat(Z, ZZ, TRANSPOSE, NORMAL);
  ZZ.AddDiagonal(lambda);
  ZZ.DPOTRF(UPPER);

  // Cleanup
  Delete_1D_Array<INTEGER>(&idx);
  Delete_1D_Array<INTEGER>(&idx2);

  // Memory estimation
  double mem_est = (double)r;
  return mem_est;

}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void KRR_Nystrom<Kernel, Point, PointArray>::
Test(const PointArray &Xtest,
     const DVector &ytrain,
     const Kernel &mKernel,
     DVector &ytest_predict) const {

  // K1 = phi(X,Y)
  DMatrix K1;
  K1.BuildKernelMatrix<Kernel, Point, PointArray>(mKernel, Xtest, Y);

  // Z0 = K1*K3
  DMatrix Z0;
  K1.MatMat(K3, Z0, NORMAL, NORMAL);

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
void KRR_Nystrom<Kernel, Point, PointArray>::
Test(const PointArray &Xtest,
     const DMatrix &Ytrain,
     const Kernel &mKernel,
     DMatrix &Ytest_predict) const {

  // K1 = phi(X,Y)
  DMatrix K1;
  K1.BuildKernelMatrix<Kernel, Point, PointArray>(mKernel, Xtest, Y);

  // Z0 = K1*K3
  DMatrix Z0;
  K1.MatMat(K3, Z0, NORMAL, NORMAL);

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
