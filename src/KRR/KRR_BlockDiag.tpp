#ifndef _KRR_BLOCKDIAG_TPP_
#define _KRR_BLOCKDIAG_TPP_


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
KRR_BlockDiag<Kernel, Point, PointArray>::
KRR_BlockDiag() {
  Init();
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void KRR_BlockDiag<Kernel, Point, PointArray>::
Init(void) {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void KRR_BlockDiag<Kernel, Point, PointArray>::
ReleaseAllMemory(void) {
  K.ReleaseAllMemory();
  invK.ReleaseAllMemory();
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
KRR_BlockDiag<Kernel, Point, PointArray>::
KRR_BlockDiag(const KRR_BlockDiag &G) {
  Init();
  DeepCopy(G);
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
KRR_BlockDiag<Kernel, Point, PointArray>& KRR_BlockDiag<Kernel, Point, PointArray>::
operator= (const KRR_BlockDiag &G) {
  if (this != &G) {
    DeepCopy(G);
  }
  return *this;
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void KRR_BlockDiag<Kernel, Point, PointArray>::
DeepCopy(const KRR_BlockDiag &G) {
  ReleaseAllMemory();
  K    = G.K;
  invK = G.invK;
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
KRR_BlockDiag<Kernel, Point, PointArray>::
~KRR_BlockDiag() {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
double KRR_BlockDiag<Kernel, Point, PointArray>::
PreTrain(PointArray &Xtrain,    // Will be permuted
         DVector &ytrain,       // Will be permuted accordingly
         INTEGER *Perm,         // Must preallocate memory
         INTEGER *iPerm,        // Must preallocate memory
         INTEGER N0,
         unsigned Seed,
         PartMethod mPar        // Partitioning method
         ) {

  // N
  INTEGER N = Xtrain.GetN();

  // Build tree. Xtrain is permuted
  K.BuildTree<Point, PointArray>(Xtrain, Perm, iPerm, N0, Seed, mPar);

  // Permute ytrain accordingly
  ytrain.Permute(Perm, N);

  // Memory estimation
  double mem_est = (double)K.GetMemEst() / N;
  return mem_est;

}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
double KRR_BlockDiag<Kernel, Point, PointArray>::
PreTrain(PointArray &Xtrain,    // Will be permuted
         DMatrix &Ytrain,       // Will be permuted accordingly
         INTEGER *Perm,         // Must preallocate memory
         INTEGER *iPerm,        // Must preallocate memory
         INTEGER N0,
         unsigned Seed,
         PartMethod mPar        // Partitioning method
         ) {

  // N
  INTEGER N = Xtrain.GetN();

  // Build tree. Xtrain is permuted
  K.BuildTree<Point, PointArray>(Xtrain, Perm, iPerm, N0, Seed, mPar);

  // Permute Ytrain accordingly
  INTEGER m = Ytrain.GetN();
  for (INTEGER i = 0; i < m; i++) {
    DVector ytrain;
    Ytrain.GetColumn(i, ytrain);
    ytrain.Permute(Perm, N);
    Ytrain.SetColumn(i, ytrain);
  }

  // Memory estimation
  double mem_est = (double)K.GetMemEst() / N;
  return mem_est;

}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void KRR_BlockDiag<Kernel, Point, PointArray>::
Train(const PointArray &Xtrain, // Must use permuted Xtrain
      const Kernel &mKernel,
      double lambda) {

  // K = phi(X,X) + lambda*I
  K.BuildKernelMatrix<Kernel, Point, PointArray>(mKernel, Xtrain, lambda);

  // invK
  K.Invert(invK, SPD);

}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void KRR_BlockDiag<Kernel, Point, PointArray>::
Test(const PointArray &Xtrain, // Must use permuted Xtrain
     const PointArray &Xtest,
     const DVector &ytrain,    // Must use permuted ytrain
     const Kernel &mKernel,
     DVector &ytest_predict) {

  // yy = y-ymean
  double ymean = ytrain.Mean();
  DVector yy;
  ytrain.Subtract(ymean, yy);

  // c = K\(y-ymean)
  DVector c;
  invK.MatVec(yy, c, NORMAL, CHOL_FACT);

  // K0 = phi(Xtest,Xtrain)
  // y = K0*c
  K.MatVecImplicit<Kernel, Point, PointArray>
    (Xtrain, Xtest, mKernel, c, ytest_predict);

  // y = K0*c + ymean
  ytest_predict.Add(ymean);

}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void KRR_BlockDiag<Kernel, Point, PointArray>::
Test(const PointArray &Xtrain, // Must use permuted Xtrain
     const PointArray &Xtest,
     const DMatrix &Ytrain,    // Must use permuted Ytrain
     const Kernel &mKernel,
     DMatrix &Ytest_predict) {

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

  // c = K\(y-ymean)
  INTEGER N = Xtrain.GetN();
  DMatrix C(N, m);
  invK.MatMat(YY, C, NORMAL, CHOL_FACT);

  // K0 = phi(Xtest,Xtrain)
  // y = K0*c
  K.MatMatImplicit<Kernel, Point, PointArray>
    (Xtrain, Xtest, mKernel, C, Ytest_predict);

  // y = K0*c + ymean
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
