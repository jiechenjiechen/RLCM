// The KRR_Standard class implements the standard kernel ridge
// regression.
//
//   y0 = K0 * inv(K+lambda*I) * (y-ymean) + ymean.
//
// The implementation of the methods in this class is parallelized.

#ifndef _KRR_STANDARD_
#define _KRR_STANDARD_

#include "../Matrices/DMatrix.hpp"

template<class Kernel, class Point, class PointArray>
class KRR_Standard {

public:

  KRR_Standard();
  KRR_Standard(const KRR_Standard &G);
  KRR_Standard& operator= (const KRR_Standard &G);
  ~KRR_Standard();

  void Init(void);
  void ReleaseAllMemory(void);
  void DeepCopy(const KRR_Standard &G);

  // This function performs factorization of K+lambda*I only. The
  // function returns the estimated memory consumption per training
  // data point according to a simplified memory metric model.
  double Train(const PointArray &Xtrain,
               const Kernel &mKernel,
               double lambda);

  // Compute prediction.
  void Test(const PointArray &Xtrain,
            const PointArray &Xtest,
            const DVector &ytrain,
            const Kernel &mKernel,
            DVector &ytest_predict) const;
  void Test(const PointArray &Xtrain,
            const PointArray &Xtest,
            const DMatrix &Ytrain,
            const Kernel &mKernel,
            DMatrix &Ytest_predict) const;

protected:

private:

  DMatrix K; // Factored form of K+lambda*I

};

#include "KRR_Standard.tpp"

#endif
