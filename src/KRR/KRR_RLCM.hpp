// The KRR_RLCM class implements the kernel ridge regression by using
// recursively low-rank compressed matrices.
//
//   y0 = K0 * inv(K+lambda*I) * (y-ymean) + ymean.
//
// The implementation of this class is parallelized.

#ifndef _KRR_RLCM_
#define _KRR_RLCM_

#include "../Matrices/CMatrix.hpp"
#include "../Solvers/GMRES.hpp"

template<class Kernel, class Point, class PointArray>
class KRR_RLCM {

public:

  KRR_RLCM();
  KRR_RLCM(const KRR_RLCM &G);
  KRR_RLCM& operator= (const KRR_RLCM &G);
  ~KRR_RLCM();

  void Init(void);
  void ReleaseAllMemory(void);
  void DeepCopy(const KRR_RLCM &G);

  // This function builds the hierarchy tree only. The function
  // returns the estimated memory consumption per training data point
  // according to a simplified memory metric model.
  double PreTrain(PointArray &Xtrain,    // Will be permuted
                  DVector &ytrain,       // Will be permuted accordingly
                  INTEGER *Perm,         // Must preallocate memory
                  INTEGER *iPerm,        // Must preallocate memory
                  INTEGER r,
                  double mDiagCorrect,   // Diagonal correction
                  bool Refinement_,      // Whether to refine the linear solve
                  unsigned Seed,
                  PartMethod mPar = RAND // Partitioning method
                  );
  double PreTrain(PointArray &Xtrain,    // Will be permuted
                  DMatrix &Ytrain,       // Will be permuted accordingly
                  INTEGER *Perm,         // Must preallocate memory
                  INTEGER *iPerm,        // Must preallocate memory
                  INTEGER r,
                  double mDiagCorrect,   // Diagonal correction
                  bool Refinement_,      // Whether to refine the linear solve
                  unsigned Seed,
                  PartMethod mPar = RAND // Partitioning method
                  );

  // This function builds the kernel matrix and inverts it.
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

  CMatrix K;       // K+lambda*I
  CMatrix invK;    // inv(K+lambda*I)
  bool Refinement; // Whether to refine the linear solve

};

#include "KRR_RLCM.tpp"

#endif
