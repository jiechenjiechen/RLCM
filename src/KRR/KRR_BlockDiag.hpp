// The KRR_BlockDiag class implements the kernel ridge regression by
// using block-diagonal approximation.
//
//   y0 = K0 * inv(K+lambda*I) * (y-ymean) + ymean.
//
// The implementation of this class is parallelized.

#ifndef _KRR_BLOCKDIAG_
#define _KRR_BLOCKDIAG_

#include "../Matrices/BMatrix.hpp"

template<class Kernel, class Point, class PointArray>
class KRR_BlockDiag {

public:

  KRR_BlockDiag();
  KRR_BlockDiag(const KRR_BlockDiag &G);
  KRR_BlockDiag& operator= (const KRR_BlockDiag &G);
  ~KRR_BlockDiag();

  void Init(void);
  void ReleaseAllMemory(void);
  void DeepCopy(const KRR_BlockDiag &G);

  // This function builds the hierarchy tree only. The function
  // returns the estimated memory consumption per training data point
  // according to a simplified memory metric model.
  double PreTrain(PointArray &Xtrain,    // Will be permuted
                  DVector &ytrain,       // Will be permuted accordingly
                  INTEGER *Perm,         // Must preallocate memory
                  INTEGER *iPerm,        // Must preallocate memory
                  INTEGER N0,
                  unsigned Seed,
                  PartMethod mPar = RAND // Partitioning method
                  );
  double PreTrain(PointArray &Xtrain,    // Will be permuted
                  DMatrix &Ytrain,       // Will be permuted accordingly
                  INTEGER *Perm,         // Must preallocate memory
                  INTEGER *iPerm,        // Must preallocate memory
                  INTEGER N0,
                  unsigned Seed,
                  PartMethod mPar = RAND // Partitioning method
                  );

  // This function builds the kernel matrix and factorizes it.
  void Train(const PointArray &Xtrain, // Must use permuted Xtrain
             const Kernel &mKernel,
             double lambda);

  // Compute prediction.
  void Test(const PointArray &Xtrain, // Must use permuted Xtrain
            const PointArray &Xtest,
            const DVector &ytrain,    // Must use permuted ytrain
            const Kernel &mKernel,
            DVector &ytest_predict);
  void Test(const PointArray &Xtrain, // Must use permuted Xtrain
            const PointArray &Xtest,
            const DMatrix &Ytrain,    // Must use permuted Ytrain
            const Kernel &mKernel,
            DMatrix &Ytest_predict);

protected:

private:

  BMatrix K;    // K+lambda*I
  BMatrix invK; // Factored form of K+lambda*I

};

#include "KRR_BlockDiag.tpp"

#endif
