#include "NoName.hpp"

#define VAL_d 2


//--------------------------------------------------------------------------
NoName::
NoName() {
  Init();
}


//--------------------------------------------------------------------------
void NoName::
Init(void) {
  ReleaseAllMemory();
  d = VAL_d;
}


//--------------------------------------------------------------------------
void NoName::
ReleaseAllMemory(void) {
}


//--------------------------------------------------------------------------
NoName::
NoName(const NoName &G) {
  Init();
  DeepCopy(G);
}


//--------------------------------------------------------------------------
NoName& NoName::
operator= (const NoName &G) {
  if (this != &G) {
    DeepCopy(G);
  }
  return *this;
}


//--------------------------------------------------------------------------
void NoName::
DeepCopy(const NoName &G) {
  d = G.d;
}


//--------------------------------------------------------------------------
NoName::
~NoName() {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
INTEGER NoName::
GetD(void) const {
  return d;
}


//--------------------------------------------------------------------------
double NoName::
Eval(const DPoint &x) const {
  if (x.GetD() != d) {
    printf("NoName::Eval. Error: Point dimension is not %ld. Return NAN.\n", (long)d);
  }
  double *mx = x.GetPointer();
  double x1 = mx[1], x2 = mx[0];
  double y = ( exp(1.4*x1) * cos(3.5*M_PI*x1) ) *
    ( sin(2*M_PI*x2) + 0.2*sin(8*M_PI*x2) );
  return y;
}
