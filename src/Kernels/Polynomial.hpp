// The Polynomial class implements the polynomial kernel
//
//     phi(x,y) = s * (a * x'y + c)^deg.
// 
// Note that the regularization/nugget term lambda is not part of the
// kernel, because some methods treat it differently. To compensate
// the missing term, we allow an additional argument lambda in the
// kernel evaluation. See Eval() below.
//
// The implementation of this class is NOT parallelized, because the
// major computation---kernel evaluation---is more likely to be used
// in the context of evaluating a kernel matrix. The parallelism
// should be used on iterating the matrix elements rather than on an
// individual evaluation.

#ifndef _POLYNOMIAL_
#define _POLYNOMIAL_

#include "../Misc/Common.hpp"

class Polynomial {

public:

  Polynomial();
  Polynomial(double s_, double a_, double c_, double deg_);
  Polynomial(const Polynomial &G);
  Polynomial& operator= (const Polynomial &G);
  ~Polynomial();

  void Init(void);
  void Init(double s_, double a_, double c_, double deg_);
  void ReleaseAllMemory(void);
  void DeepCopy(const Polynomial &G);

  // Kernel properties
  static std::string const GetKernelName(void) { return "Polynomial"; }
  bool IsSymmetric(void) const { return true; }

  // Get kernel parameters
  double GetS(void) const;
  double GetA(void) const;
  double GetC(void) const;
  double GetDeg(void) const;

  // Evaluate the kernel for a pair of points x and y.
  // The class Point must have the following methods:
  //
  //   INTEGER GetD(void) const;
  //   double InProd(const Point &y) const;
  //
  template<class Point>
  double Eval(const Point &x, const Point &y, double lambda = 0.0) const;

protected:

private:

  // Kernel parameters
  double s;
  double a;
  double c;
  double deg;

};

#include "Polynomial.tpp"

#endif
