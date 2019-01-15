// The KRR_Nystrom class implements the Nystrom method for performing
// approximate kernel ridge regression.
//
// Given training points X (n*d), labels y (n*1), and rank r, the
// Nystrom feature Z (n*r) is computed in the following way:
//
//   kernel matrix Phi approx= Z*Z' = K1*pinv(K2)*K1,
//
// where K1 is the data*center matrix and K2 is the center*center
// matrix. The r centers are found by using random sampling.
//
// For testing points X0, the output y0 is computed as
//
//   y0 = Z0 * inv(Z'*Z+lambda*I) * Z'*(y-ymean) + ymean.
//
// The implementation of this class is parallelized.

#ifndef _KRR_NYSTROM_
#define _KRR_NYSTROM_

#include "../Matrices/DMatrix.hpp"

template<class Kernel, class Point, class PointArray>
class KRR_Nystrom {

public:

  KRR_Nystrom();
  KRR_Nystrom(const KRR_Nystrom &G);
  KRR_Nystrom& operator= (const KRR_Nystrom &G);
  ~KRR_Nystrom();

  void Init(void);
  void ReleaseAllMemory(void);
  void DeepCopy(const KRR_Nystrom &G);

  // This function performs factorization of Z'*Z+lambda*I only. The
  // function returns the estimated memory consumption per training
  // data point according to a simplified memory metric model.
  double Train(const PointArray &Xtrain,
               const Kernel &mKernel,
               double lambda,
               INTEGER r,
               unsigned Seed);

  // Compute prediction.
  void Test(const PointArray &Xtest,
            const DVector &ytrain,
            const Kernel &mKernel,
            DVector &ytest_predict) const;
  void Test(const PointArray &Xtest,
            const DMatrix &Ytrain,
            const Kernel &mKernel,
            DMatrix &Ytest_predict) const;

protected:

private:

  PointArray Y; // Centers
  DMatrix K3;   // sqrt(pinv(K2))
  DMatrix Z;    // Z
  DMatrix ZZ;   // Factored form of Z'*Z+lambda*I

};

#include "KRR_Nystrom.tpp"

#endif
