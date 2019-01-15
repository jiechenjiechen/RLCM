// The Sampling_RLCM class implements GP sampling based on RLCM sqrt
// factorization.
//
//   y = sqrt(K) * randn(n,1).
//
// The implementation of this class is parallelized.

#ifndef _Sampling_RLCM_
#define _Sampling_RLCM_

#include "../Matrices/CMatrix.hpp"

template<class Kernel, class Point, class PointArray>
class Sampling_RLCM {

public:

  Sampling_RLCM();
  Sampling_RLCM(const Sampling_RLCM &G);
  Sampling_RLCM& operator= (const Sampling_RLCM &G);
  ~Sampling_RLCM();

  void Init(void);
  void ReleaseAllMemory(void);
  void DeepCopy(const Sampling_RLCM &G);

  // Sampling
  void GenRandomField(CMatrix &K,  // Must have called CMatrix::BuildTree()
                      const PointArray &X,  // Must use permuted X
                      DVector &y,           // In the permuting order of X
                      const Kernel &mKernel,
                      double lambda) const;

protected:

private:

};

#include "Sampling_RLCM.tpp"

#endif
