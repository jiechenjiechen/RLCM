// The MLE_Standard class implements the standard likelihood
// calculation.
//
//   Loglik = -1/2*y'*inv(K)*y -1/2*logdet(K) -n/2*log(2*pi).
//
// The implementation of this class is parallelized.

#ifndef _MLE_STANDARD_
#define _MLE_STANDARD_

#include "../Matrices/DMatrix.hpp"

template<class Kernel, class Point, class PointArray>
class MLE_Standard {

public:

  MLE_Standard();
  MLE_Standard(const MLE_Standard &G);
  MLE_Standard& operator= (const MLE_Standard &G);
  ~MLE_Standard();

  void Init(void);
  void ReleaseAllMemory(void);
  void DeepCopy(const MLE_Standard &G);

  // Loglikelihood
  double LogLik(const PointArray &X,
                const DVector &y,
                const Kernel &mKernel,
                double lambda) const;

protected:

private:

};

#include "MLE_Standard.tpp"

#endif
