// The ProductLaplace class implements the product Laplace kernel
//
//     phi(x,y) = s * prod_i exp(-r_i/sigma)
// 
// where r_i = abs(x_i-y_i). Note that the regularization/nugget term
// lambda is not part of the kernel, because some methods treat it
// differently. To compensate the missing term, we allow an additional
// argument lambda in the kernel evaluation. See Eval() below.
//
// The implementation of this class is NOT parallelized, because the
// major computation---kernel evaluation---is more likely to be used
// in the context of evaluating a kernel matrix. The parallelism
// should be used on iterating the matrix elements rather than on an
// individual evaluation.

#ifndef _PRODUCT_LAPLACE_
#define _PRODUCT_LAPLACE_

#include "../Misc/Common.hpp"

class ProductLaplace {

public:

  ProductLaplace();
  ProductLaplace(double s_, double sigma_);
  ProductLaplace(const ProductLaplace &G);
  ProductLaplace& operator= (const ProductLaplace &G);
  ~ProductLaplace();

  void Init(void);
  void Init(double s_, double sigma_);
  void ReleaseAllMemory(void);
  void DeepCopy(const ProductLaplace &G);

  // Kernel properties
  static std::string const GetKernelName(void) { return "ProductLaplace"; }
  bool IsSymmetric(void) const { return true; }

  // Get kernel parameters
  double GetS(void) const;
  double GetSigma(void) const;

  // Evaluate the kernel for a pair of points x and y.
  // The class Point must have the following methods:
  //
  //   INTEGER GetD(void) const;
  //   double Dist1(const Point &y) const;
  //
  template<class Point>
  double Eval(const Point &x, const Point &y, double lambda = 0.0) const;

protected:

private:

  // Kernel parameters
  double s;
  double sigma;

};

#include "ProductLaplace.tpp"

#endif
