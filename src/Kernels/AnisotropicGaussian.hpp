// The AnisotropicGaussian class implements the anisotropic Gaussian
// kernel
//
//     phi(x,y) = s * exp(-r^2/2)
// 
// where r = sum_i [(x_i-y_i)/sigma_i]^2. Note that the
// regularization/nugget term lambda is not part of the kernel,
// because some methods treat it differently. To compensate this
// missing term, we allow an additional argument lambda in the kernel
// evaluation. See Eval() below.
//
// The implementation of this class is NOT parallelized, because the
// major computation---kernel evaluation---is more likely to be used
// in the context of evaluating a kernel matrix. The parallelism
// should be used on iterating the matrix elements rather than on an
// individual evaluation.

#ifndef _ANISOTROPIC_GAUSSIAN_
#define _ANISOTROPIC_GAUSSIAN_

#include "../Misc/Common.hpp"

class AnisotropicGaussian {

public:

  AnisotropicGaussian();
  AnisotropicGaussian(double s_, INTEGER d_, const double *sigma_);
  AnisotropicGaussian(const AnisotropicGaussian &G);
  AnisotropicGaussian& operator= (const AnisotropicGaussian &G);
  ~AnisotropicGaussian();

  void Init(void);
  void Init(double s_, INTEGER d_, const double *sigma_);
  void ReleaseAllMemory(void);
  void DeepCopy(const AnisotropicGaussian &G);

  // Kernel properties
  static std::string const GetKernelName(void) { return "AnisotropicGaussian"; }
  bool IsSymmetric(void) const { return true; }

  // Get kernel parameters
  double GetS(void) const;
  INTEGER GetD(void) const;
  double* GetSigma(void) const;

  // Evaluate the kernel for a pair of points x and y.
  // The class Point must have the following methods:
  //
  //   INTEGER GetD(void) const;
  //   double Dist2(const Point &y) const;
  //
  template<class Point>
  double Eval(const Point &x, const Point &y, double lambda = 0.0) const;

protected:

private:

  // Kernel parameters
  double s;
  INTEGER d;
  double *sigma;

};

#include "AnisotropicGaussian.tpp"

#endif
