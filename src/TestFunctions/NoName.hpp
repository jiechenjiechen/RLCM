// The NoName class implements the following test function (2D)
//
//   z = ( exp(1.4*x1) * cos(3.5*pi*x1) ) * ( sin(2*pi*x2) + 0.2*sin(8*pi*x2) ).
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

#ifndef _NO_NAME_
#define _NO_NAME_

#include "../Misc/Common.hpp"
#include "../Matrices/DPoint.hpp"

class NoName {

public:

  NoName();
  NoName(const NoName &G);
  NoName& operator= (const NoName &G);
  ~NoName();

  void Init(void);
  void ReleaseAllMemory(void);
  void DeepCopy(const NoName &G);

  // Properties
  INTEGER GetD(void) const;

  // Evaluate
  double Eval(const DPoint &x) const;

protected:

private:

  INTEGER d; // Dimension (must be 2)

};

#endif
