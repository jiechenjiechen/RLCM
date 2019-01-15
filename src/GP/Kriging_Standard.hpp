// The Kriging_Standard class implements the standard kriging method
// (kernel ridge regression).
//
//   y0 = K0 * inv(K+lambda*I) * y
//   var(y0) = K(0)+lambda - K0 * inv(K+lambda*I) * K0
//
// The implementation of this class is parallelized.

#ifndef _KRIGING_STANDARD_
#define _KRIGING_STANDARD_

#include "../Matrices/DMatrix.hpp"

template<class Kernel, class Point, class PointArray>
class Kriging_Standard {

public:

  Kriging_Standard();
  Kriging_Standard(const Kriging_Standard &G);
  Kriging_Standard& operator= (const Kriging_Standard &G);
  ~Kriging_Standard();

  void Init(void);
  void ReleaseAllMemory(void);
  void DeepCopy(const Kriging_Standard &G);

  // This function performs factorization of K+lambda*I only.
  void Train(const PointArray &Xtrain,
             const Kernel &mKernel,
             double lambda);

  // Compute prediction and variance.
  void Test(const PointArray &Xtrain,
            const PointArray &Xtest,
            const DVector &ytrain,
            const Kernel &mKernel,
            double lambda,
            DVector &ytest_predict,
            DVector &ytest_stddev) const;

  // The split of Train() and Test() is used to handle the case when
  // there is a huge set of test points, in which case they need be
  // packed in batches and Test() is called repeatedly.

protected:

private:

  DMatrix K; // Factored form of K+lambda*I

};

#include "Kriging_Standard.tpp"

#endif
