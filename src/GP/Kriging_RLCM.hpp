// The Kriging_RLCM class implements kriging (kernel ridge regression)
// by using recursively low-rank compressed matrices.
//
//   y0 = K0 * inv(K+lambda*I) * y
//   var(y0) = K(0)+lambda - K0 * inv(K+lambda*I) * K0
//
// The implementation of this class is parallelized.

#ifndef _KRIGING_RLCM_
#define _KRIGING_RLCM_

#include "../Matrices/CMatrix.hpp"
#include "../Solvers/GMRES.hpp"

template<class Kernel, class Point, class PointArray>
class Kriging_RLCM {

public:

  Kriging_RLCM();
  Kriging_RLCM(const Kriging_RLCM &G);
  Kriging_RLCM& operator= (const Kriging_RLCM &G);
  ~Kriging_RLCM();

  void Init(void);
  void ReleaseAllMemory(void);
  void DeepCopy(const Kriging_RLCM &G);

  // This function builds the kernel matrix and inverts it.
  void Train(CMatrix &K,  // Must have called CMatrix::BuildTree()
             const PointArray &Xtrain, // Must use permuted Xtrain
             const Kernel &mKernel,
             double lambda);

  // Compute prediction and variance.
  void Test(CMatrix &K,  // Must be the same as that input to Train()
            const PointArray &Xtrain, // Must use permuted Xtrain
            const PointArray &Xtest,
            const DVector &ytrain,    // Must use permuted ytrain
            const Kernel &mKernel,
            double lambda,
            DVector &ytest_predict,
            DVector &ytest_stddev);

  // The split of Train() and Test() is used to handle the case when
  // there is a huge set of test points, in which case they need be
  // packed in batches and Test() is called repeatedly.

protected:

private:

  CMatrix invK; // Inverse of K+lambda*I

};

#include "Kriging_RLCM.tpp"

#endif
