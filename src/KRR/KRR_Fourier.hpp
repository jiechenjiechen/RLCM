// The KRR_Fourier class implements the Fourier method for performing
// approximate kernel ridge regression.
//
// Given training points X (n*d), labels y (n*1), and rank r, the
// Fourier feature Z (n*r) is computed in the following way:
//
//   kernel matrix Phi approx= Z*Z',  Z = cos(X*w+b) * sqrt(2/r) * sqrt(s),
//
// where each column of w is iid sampled from the Fourier transform of
// the kernel and each entry of b is sampled from Uniform(0,2*pi).
//
// For testing points X0, the output y0 is computed as
//
//   y0 = Z0 * inv(Z'*Z + lambda) * Z'*(y-ymean) + ymean.
//
// The implementation of this class is parallelized.

#ifndef _KRR_FOURIER_
#define _KRR_FOURIER_

#include "../Matrices/DMatrix.hpp"

template<class Kernel, class Point, class PointArray>
class KRR_Fourier {

public:

  KRR_Fourier();
  KRR_Fourier(const KRR_Fourier &G);
  KRR_Fourier& operator= (const KRR_Fourier &G);
  ~KRR_Fourier();

  void Init(void);
  void ReleaseAllMemory(void);
  void DeepCopy(const KRR_Fourier &G);

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

  DMatrix w;  // Random w
  DVector b;  // Random b
  DMatrix Z;  // Z
  DMatrix ZZ; // Factored form of Z'*Z+lambda*I

};

#include "KRR_Fourier.tpp"

#endif
