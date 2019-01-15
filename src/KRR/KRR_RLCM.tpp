#ifndef _KRR_RLCM_TPP_
#define _KRR_RLCM_TPP_

#define INITVAL_Refinement false


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
KRR_RLCM<Kernel, Point, PointArray>::
KRR_RLCM() {
  Init();
  Refinement = INITVAL_Refinement;
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void KRR_RLCM<Kernel, Point, PointArray>::
Init(void) {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void KRR_RLCM<Kernel, Point, PointArray>::
ReleaseAllMemory(void) {
  K.ReleaseAllMemory();
  invK.ReleaseAllMemory();
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
KRR_RLCM<Kernel, Point, PointArray>::
KRR_RLCM(const KRR_RLCM &G) {
  Init();
  DeepCopy(G);
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
KRR_RLCM<Kernel, Point, PointArray>& KRR_RLCM<Kernel, Point, PointArray>::
operator= (const KRR_RLCM &G) {
  if (this != &G) {
    DeepCopy(G);
  }
  return *this;
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void KRR_RLCM<Kernel, Point, PointArray>::
DeepCopy(const KRR_RLCM &G) {
  ReleaseAllMemory();
  K          = G.K;
  invK       = G.invK;
  Refinement = G.Refinement;
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
KRR_RLCM<Kernel, Point, PointArray>::
~KRR_RLCM() {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
double KRR_RLCM<Kernel, Point, PointArray>::
PreTrain(PointArray &Xtrain,    // Will be permuted
         DVector &ytrain,       // Will be permuted accordingly
         INTEGER *Perm,         // Must preallocate memory
         INTEGER *iPerm,        // Must preallocate memory
         INTEGER r,
         double mDiagCorrect,   // Diagonal correction
         bool Refinement_,      // Whether to refine the linear solve
         unsigned Seed,
         PartMethod mPar        // Partitioning method
         ) {

  // N
  INTEGER N = Xtrain.GetN();

  // N0
  INTEGER N0 = r;

  // Refinement
  Refinement = Refinement_;

  // Build tree. Xtrain is permuted
  K.BuildTree<Point, PointArray>
    (Xtrain, Perm, iPerm, r, N0, mDiagCorrect, Seed, mPar);

  // Permute ytrain accordingly
  ytrain.Permute(Perm, N);

  // Memory estimation
  double mem_est = (double)K.GetMemEst() / N;
  return mem_est;

}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
double KRR_RLCM<Kernel, Point, PointArray>::
PreTrain(PointArray &Xtrain,    // Will be permuted
         DMatrix &Ytrain,       // Will be permuted accordingly
         INTEGER *Perm,         // Must preallocate memory
         INTEGER *iPerm,        // Must preallocate memory
         INTEGER r,
         double mDiagCorrect,   // Diagonal correction
         bool Refinement_,      // Whether to refine the linear solve
         unsigned Seed,
         PartMethod mPar        // Partitioning method
         ) {

  // N
  INTEGER N = Xtrain.GetN();

  // N0
  INTEGER N0 = r;

  // Refinement
  Refinement = Refinement_;

  // Build tree. Xtrain is permuted
  K.BuildTree<Point, PointArray>
    (Xtrain, Perm, iPerm, r, N0, mDiagCorrect, Seed, mPar);

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
void KRR_RLCM<Kernel, Point, PointArray>::
Train(const PointArray &Xtrain, // Must use permuted Xtrain
      const Kernel &mKernel,
      double lambda) {

  // K = phi(X,X) + lambda*I
  K.BuildKernelMatrix<Kernel, Point, PointArray>(mKernel, Xtrain, lambda);

  // invK
  K.Invert(invK);

}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void KRR_RLCM<Kernel, Point, PointArray>::
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
  invK.MatVec(yy, c, NORMAL);

  // Refine the linear solves
  if (Refinement) {
    GMRES mGMRES;
    mGMRES.RefineSolve<CMatrix, CMatrix>(K, UNFACT, yy, c, invK, UNFACT);
  }

  // K0 = phi(Xtest,Xtrain)
  // y = K0*c
  K.MatVecImplicit<Kernel, Point, PointArray>
    (Xtrain, Xtest, mKernel, c, ytest_predict);

  // y = K0*c + ymean
  ytest_predict.Add(ymean);

}


//--------------------------------------------------------------------------
template<class Kernel, class Point, class PointArray>
void KRR_RLCM<Kernel, Point, PointArray>::
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
  DMatrix C;
  invK.MatMat(YY, C, NORMAL);

  // Refine the linear solves
  if (Refinement) {
    GMRES mGMRES;
    for (INTEGER i = 0; i < m; i++) {
      DVector yy, c;
      YY.GetColumn(i, yy);
      C.GetColumn(i, c);
      mGMRES.RefineSolve<CMatrix, CMatrix>(K, UNFACT, yy, c, invK, UNFACT);
      mGMRES.GetSolution(c);
      C.SetColumn(i, c);
    }
  }

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
