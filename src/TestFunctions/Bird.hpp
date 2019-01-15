// The Bird class implements Mishra's Bird test function (2D)
//
//   z = sin(y).*exp((1-cos(x)).^2) + 
//       cos(x).*exp((1-sin(y)).^2) + (x-y).^2.
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

#ifndef _BIRD_
#define _BIRD_

#include "../Misc/Common.hpp"
#include "../Matrices/DPoint.hpp"

class Bird {

public:

  Bird();
  Bird(const Bird &G);
  Bird& operator= (const Bird &G);
  ~Bird();

  void Init(void);
  void ReleaseAllMemory(void);
  void DeepCopy(const Bird &G);

  // Properties
  INTEGER GetD(void) const;

  // Evaluate
  double Eval(const DPoint &x) const;

protected:

private:

  INTEGER d; // Dimension (must be 2)

};

#endif
