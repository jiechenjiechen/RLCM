#include "Polynomial.hpp"

#define INITVAL_s   DBL_MAX
#define INITVAL_a   DBL_MAX
#define INITVAL_c   DBL_MAX
#define INITVAL_deg DBL_MAX


//--------------------------------------------------------------------------
Polynomial::
Polynomial() {
  Init();
}


//--------------------------------------------------------------------------
void Polynomial::
Init(void) {
  ReleaseAllMemory();
  s = INITVAL_s;
  a = INITVAL_a;
  c = INITVAL_c;
  deg = INITVAL_deg;
}


//--------------------------------------------------------------------------
void Polynomial::
ReleaseAllMemory(void) {
}


//--------------------------------------------------------------------------
void Polynomial::
Init(double s_, double a_, double c_, double deg_) {
  s = s_;
  a = a_;
  c = c_;
  deg = deg_;
}


//--------------------------------------------------------------------------
Polynomial::
Polynomial(double s_, double a_, double c_, double deg_) {
  Init(s_, a_, c_, deg_);
}


//--------------------------------------------------------------------------
Polynomial::
Polynomial(const Polynomial &G) {
  Init();
  DeepCopy(G);
}


//--------------------------------------------------------------------------
Polynomial& Polynomial::
operator= (const Polynomial &G) {
  if (this != &G) {
    DeepCopy(G);
  }
  return *this;
}


//--------------------------------------------------------------------------
void Polynomial::
DeepCopy(const Polynomial &G) {
  s = G.s;
  a = G.a;
  c = G.c;
  deg = G.deg;
}


//--------------------------------------------------------------------------
Polynomial::
~Polynomial() {
  ReleaseAllMemory();
}


//--------------------------------------------------------------------------
double Polynomial::
GetS(void) const {
  return s;
}


//--------------------------------------------------------------------------
double Polynomial::
GetA(void) const {
  return a;
}


//--------------------------------------------------------------------------
double Polynomial::
GetC(void) const {
  return c;
}


//--------------------------------------------------------------------------
double Polynomial::
GetDeg(void) const {
  return deg;
}
