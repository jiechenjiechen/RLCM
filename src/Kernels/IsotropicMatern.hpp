// The IsotropicMatern class implements the isotropic Matern kernel
//
//                     r^nu * besselk(nu,r)
//     phi(x,y) = s * ----------------------
//                     2^(nu-1) * gamma(nu)
// 
// where r is the elliptical radius
// 
//     r = sqrt[2*nu] * sqrt[ sum [ (x_i-y_i)^2 / ell^2 ]]
//
// and nu > 0. Note that the regularization/nugget term lambda is not
// part of the kernel, because some methods treat it differently. To
// compensate this missing term, we allow an additional argument
// lambda in the kernel evaluation. See Eval() below.
//
// The implementation of this class is NOT parallelized, because the
// major computation---kernel evaluation---is more likely to be used
// in the context of evaluating a kernel matrix. The parallelism
// should be used on iterating the matrix elements rather than on an
// individual evaluation.

#ifndef _ISOTROPIC_MATERN_
#define _ISOTROPIC_MATERN_

#include "../Misc/Common.hpp"

class IsotropicMatern {

public:

  IsotropicMatern();
  IsotropicMatern(double s_, double nu_, double ell_);
  IsotropicMatern(const IsotropicMatern &G);
  IsotropicMatern& operator= (const IsotropicMatern &G);
  ~IsotropicMatern();

  void Init(void);
  void Init(double s_, double nu_, double ell_);
  void ReleaseAllMemory(void);
  void DeepCopy(const IsotropicMatern &G);

  // Kernel properties
  static std::string const GetKernelName(void) { return "IsotropicMatern"; }
  bool IsSymmetric(void) const { return true; }

  // Get kernel parameters
  double GetS(void) const;
  double GetNu(void) const;
  double GetEll(void) const;

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
  double nu;
  double ell;

  // Internally computed variables
  double k0;      // = pow(2.0, nu-1) * tgamma(nu);
  double sqrt2nu; // = sqrt(2.0*nu);

};

#include "IsotropicMatern.tpp"

#endif
