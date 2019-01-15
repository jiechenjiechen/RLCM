// The ChordalMatern class implements the isotropic Matern kernel
//
//                     r^nu * besselk(nu,r)
//     phi(x,x') = s * ----------------------
//                     2^(nu-1) * gamma(nu)
// 
// where r is the chordal distance
// 
//     r = 2 * sqrt[ sin((lat - lat')/2)^2
//                  + cos(lat) * cos(lat') * sin((lon - lon')/2)^2 ]
//           * sqrt[2*nu] / ell
//
// and nu > 0. Note that the regularization/nugget term lambda is not
// part of the kernel, because some methods treat it differently. To
// compensate this missing term, we allow an additional argument
// lambda in the kernel evaluation. See Eval() below.
//
// Each data point x is a 2-dimensional vector (lat, lon); that is,
// lat is the first dimension and it varies the fastest. Because lat
// and lon are spherical coordinates, the range of lat is [-pi/2,
// pi/2] and that of lon is [-pi, pi], although the kernel still
// returns a value even when they are outside the range.
//
// The implementation of this class is NOT parallelized, because the
// major computation---kernel evaluation---is more likely to be used
// in the context of evaluating a kernel matrix. The parallelism
// should be used on iterating the matrix elements rather than on an
// individual evaluation.

#ifndef _CHORDAL_MATERN_
#define _CHORDAL_MATERN_

#include "../Misc/Common.hpp"

class ChordalMatern {

public:

  ChordalMatern();
  ChordalMatern(double s_, double nu_, double ell_);
  ChordalMatern(const ChordalMatern &G);
  ChordalMatern& operator= (const ChordalMatern &G);
  ~ChordalMatern();

  void Init(void);
  void Init(double s_, double nu_, double ell_);
  void ReleaseAllMemory(void);
  void DeepCopy(const ChordalMatern &G);

  // Kernel properties
  static std::string const GetKernelName(void) { return "ChordalMatern"; }
  bool IsSymmetric(void) const { return true; }

  // Get kernel parameters
  double GetS(void) const;
  double GetNu(void) const;
  double GetEll(void) const;

  // Evaluate the kernel for a pair of points x and y. They must use
  // the DPoint class. This class must have the following methods:
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
  double nu;
  double ell;

  // Internally computed variables
  double k0;      // = pow(2.0, nu-1) * tgamma(nu);
  double sqrt2nu; // = sqrt(2.0*nu);

};

#include "ChordalMatern.tpp"

#endif
