// The Chi2 class implements the Chi-squared kernel
//
//     phi(x,y) = s * sum(2*x.*y./(x+y)),
// 
// where x,y >= 0. Note that the regularization/nugget term lambda is
// not part of the kernel, because some methods treat it
// differently. To compensate the missing term, we allow an additional
// argument lambda in the kernel evaluation. See Eval() below.
//
// The implementation of this class is NOT parallelized, because the
// major computation---kernel evaluation---is more likely to be used
// in the context of evaluating a kernel matrix. The parallelism
// should be used on iterating the matrix elements rather than on an
// individual evaluation.

#ifndef _CHI2_
#define _CHI2_

#include "../Misc/Common.hpp"

class Chi2 {

public:

  Chi2();
  Chi2(double s_);
  Chi2(const Chi2 &G);
  Chi2& operator= (const Chi2 &G);
  ~Chi2();

  void Init(void);
  void Init(double s_);
  void ReleaseAllMemory(void);
  void DeepCopy(const Chi2 &G);

  // Kernel properties
  static std::string const GetKernelName(void) { return "Chi2"; }
  bool IsSymmetric(void) const { return true; }

  // Get kernel parameters
  double GetS(void) const;

  // Evaluate the kernel for a pair of points x and y.
  // The class Point must have the following methods:
  //
  //   INTEGER GetD(void) const;
  //   double KernelFuncChi2(const Point &y) const;
  //
  template<class Point>
  double Eval(const Point &x, const Point &y, double lambda = 0.0) const;

protected:

private:

  // Kernel parameters
  double s;

};

#include "Chi2.tpp"

#endif
