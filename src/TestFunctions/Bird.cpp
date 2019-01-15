#include "Bird.hpp"

#define VAL_d 2


//--------------------------------------------------------------------------
Bird::
Bird() {
  Init();
}


//--------------------------------------------------------------------------
void Bird::
Init(void) {
  ReleaseAllMemory();
  d = VAL_d;
}


//--------------------------------------------------------------------------
void Bird::
ReleaseAllMemory(void) {
}


//--------------------------------------------------------------------------
Bird::
Bird(const Bird &G) {
  Init();
  DeepCopy(G);
}


//--------------------------------------------------------------------------
Bird& Bird::
operator= (const Bird &G) {
  if (this != &G) {
    DeepCopy(G);
  }
  return *this;
}


//--------------------------------------------------------------------------
void Bird::
DeepCopy(const Bird &G) {
  d = G.d;
}


//--------------------------------------------------------------------------
Bird::
~Bird() {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
INTEGER Bird::
GetD(void) const {
  return d;
}


//--------------------------------------------------------------------------
double Bird::
Eval(const DPoint &x) const {
  if (x.GetD() != d) {
    printf("Bird::Eval. Error: Point dimension is not %ld. Return NAN.\n", (long)d);
  }
  double *mx = x.GetPointer();
  double x1 = mx[1], x2 = mx[0];
  double y = sin(x2) * exp(Square(1-cos(x1))) + 
    cos(x1) * exp(Square(1-sin(x2))) + Square(x1-x2);

  return y;
}
