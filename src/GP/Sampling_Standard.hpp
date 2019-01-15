// The Sampling_Standard class implements the standard GP sampling.
//
//   y = chol(K, 'lower') * randn(n,1).
//
// The implementation of this class is parallelized.

#ifndef _Sampling_STANDARD_
#define _Sampling_STANDARD_

#include "../Matrices/DMatrix.hpp"

template<class Kernel, class Point, class PointArray>
class Sampling_Standard {

public:

  Sampling_Standard();
  Sampling_Standard(const Sampling_Standard &G);
  Sampling_Standard& operator= (const Sampling_Standard &G);
  ~Sampling_Standard();

  void Init(void);
  void ReleaseAllMemory(void);
  void DeepCopy(const Sampling_Standard &G);

  // Sampling
  void GenRandomField(const PointArray &X,
                      DVector &y,
                      const Kernel &mKernel,
                      double lambda) const;

protected:

private:

};

#include "Sampling_Standard.tpp"

#endif
