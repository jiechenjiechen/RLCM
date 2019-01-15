// The DPoint class implements a dense data point. The point is
// treated as a dense vector.
//
// The implementation of this class is NOT parallelized, because a
// point is more likely to be used in the context of a point
// array. The parallelism should be used on iterating the points
// rather than on the calculation of an individual point.

#ifndef _DPOINT_
#define _DPOINT_

#include "../Misc/Common.hpp"

class DPoint {

public:

  DPoint();
  DPoint(INTEGER d_);
  DPoint(const DPoint& G);
  DPoint& operator= (const DPoint& G);
  ~DPoint();

  void Init(void);
  void Init(INTEGER d_); // Initialized as a zero vector
  void ReleaseAllMemory(void);
  void DeepCopy(const DPoint &G);

  //-------------------- Utilities --------------------
  //
  // x is the self data point.

  // Get dimension
  INTEGER GetD(void) const;

  // Get the pointer
  double* GetPointer(void) const;

  // Set data
  void SetPoint(const double *x_, INTEGER d_);

  // x = randn()
  void SetStandardNormal(void);

  // Print the data point
  void PrintPoint(const char *name) const;

  //-------------------- Computations --------------------
  //
  // x is the self data point.

  // x = x / norm(x,2)
  void Normalize(void);

  // x'*y
  double InProd(const DPoint &y) const;

  // sum_i |x_i-y_i|
  double Dist1(const DPoint &y) const;

  // sum_i |x_i-y_i|/|sigma_i|
  double Dist1(const DPoint &y, const double *sigma) const;

  // sum_i |x_i-y_i|^2
  double Dist2(const DPoint &y) const;
  
  // sum_i (|x_i-y_i|/|sigma_i|)^2
  double Dist2(const DPoint &y, const double *sigma) const;
  
  // x = x-y  or  z = x-y
  void Subtract(const DPoint &y);
  void Subtract(const DPoint &y, DPoint &z) const;

  // x = (x+y)/2  or  z = (x+y)/2
  void AverageWith(const DPoint &y);
  void AverageWith(const DPoint &y, DPoint &z) const;

  // sum_i 2 * x_i * y_i / (x_i + y_i)
  double KernelFuncChi2(const DPoint &y) const;

protected:

private:

  INTEGER d; // Dimension
  double *x; // Point data

};

#endif
