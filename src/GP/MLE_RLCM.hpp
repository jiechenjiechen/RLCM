// The MLE_Standard class implements the likelihood calculation by
// using recursively low-rank compressed matrices.
//
//   Loglik = -1/2*y'*inv(K)*y -1/2*logdet(K) -n/2*log(2*pi).
//
// N samples:
//
//   Loglik = -1/2 sum_i y_i'*inv(K)*y_i -N/2*logdet(K) -Nn/2*log(2*pi).
//
// The implementation of this class is parallelized.

#ifndef _MLE_RLCM_
#define _MLE_RLCM_

#include "../Matrices/CMatrix.hpp"
#include "../Solvers/GMRES.hpp"

template<class Kernel, class Point, class PointArray>
class MLE_RLCM {

public:

  MLE_RLCM();
  MLE_RLCM(const MLE_RLCM &G);
  MLE_RLCM& operator= (const MLE_RLCM &G);
  ~MLE_RLCM();

  void Init(void);
  void ReleaseAllMemory(void);
  void DeepCopy(const MLE_RLCM &G);

  // Loglikelihood
  double LogLik(CMatrix &K,           // Must have called CMatrix::BuildTree()
                const PointArray &X,  // Must use permuted X
                const DVector &y,     // Must use permuted y
                const Kernel &mKernel,
                double lambda) const;
  double LogLik(CMatrix &K,           // Must have called CMatrix::BuildTree()
                const PointArray &X,  // Must use permuted X
                const DMatrix &Y,     // Must use permuted Y
                const Kernel &mKernel,
                double lambda) const;

protected:

private:

};

#include "MLE_RLCM.tpp"

#endif
