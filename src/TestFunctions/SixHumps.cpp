#include "SixHumps.hpp"

#define VAL_d 2


//--------------------------------------------------------------------------
SixHumps::
SixHumps() {
  Init();
}


//--------------------------------------------------------------------------
void SixHumps::
Init(void) {
  ReleaseAllMemory();
  d = VAL_d;
}


//--------------------------------------------------------------------------
void SixHumps::
ReleaseAllMemory(void) {
}


//--------------------------------------------------------------------------
SixHumps::
SixHumps(const SixHumps &G) {
  Init();
  DeepCopy(G);
}


//--------------------------------------------------------------------------
SixHumps& SixHumps::
operator= (const SixHumps &G) {
  if (this != &G) {
    DeepCopy(G);
  }
  return *this;
}


//--------------------------------------------------------------------------
void SixHumps::
DeepCopy(const SixHumps &G) {
  d = G.d;
}


//--------------------------------------------------------------------------
SixHumps::
~SixHumps() {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
INTEGER SixHumps::
GetD(void) const {
  return d;
}


//--------------------------------------------------------------------------
double SixHumps::
Eval(const DPoint &x) const {
  if (x.GetD() != d) {
    printf("SixHumps::Eval. Error: Point dimension is not %ld. Return NAN.\n", (long)d);
  }
  double *mx = x.GetPointer();
  double x1 = mx[1], x2 = mx[0];
  double y = (4.0 - 2.1 * Square(x1) + Square(Square(x1))/3.0) * Square(x1) +
    x1 * x2 + (-4.0 + 4.0 * Square(x2)) * Square(x2);
  return y;
}
