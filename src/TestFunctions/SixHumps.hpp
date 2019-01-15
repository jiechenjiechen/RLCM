// The SixHumps class implements the six-humps test function (2D)
//
//   z = (4-2.1*x.^2+x.^4/3).*x.^2 + x.*y + (-4+4*y.^2).*y.^2.
//
// Note that the x-y ordering is different from the d-dimensional
// ordering. In particular,
//
//   x  <->  p[1]
//   y  <->  p[0].
//
// The implementation of this class is NOT parallelized, because the
// major computation---function evaluation---is more likely to be used
// in the context of creating a vector. The parallelism should be used
// on iterating the vector elements rather than on an individual
// evaluation.

#ifndef _SIX_HUMPS_
#define _SIX_HUMPS_

#include "../Misc/Common.hpp"
#include "../Matrices/DPoint.hpp"

class SixHumps {

public:

  SixHumps();
  SixHumps(const SixHumps &G);
  SixHumps& operator= (const SixHumps &G);
  ~SixHumps();

  void Init(void);
  void ReleaseAllMemory(void);
  void DeepCopy(const SixHumps &G);

  // Properties
  INTEGER GetD(void) const;

  // Evaluate
  double Eval(const DPoint &x) const;

protected:

private:

  INTEGER d; // Dimension (must be 2)

};

#endif
