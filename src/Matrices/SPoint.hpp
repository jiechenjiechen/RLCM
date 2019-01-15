// The SPoint class implements a sparse data point. The point is
// treated as a sparse vector, where d is the length, nnz is the
// number of nonzeros, x is an array storing the nonzero elements, and
// idx is an array storing the indices of these nonzeros. For more
// details, see the class definition below.
//
// The implementation of this class is NOT parallelized, because a
// point is more likely to be used in the context of a point
// array. The parallelism should be used on iterating the points
// rather than on the calculation of an individual point.

#ifndef _SPOINT_
#define _SPOINT_

#include "DPoint.hpp"
#include "../Misc/spblas.hpp"

class SPoint {

public:

  SPoint();
  SPoint(INTEGER d_, INTEGER nnz_);
  SPoint(const SPoint& G);
  SPoint& operator= (const SPoint& G);
  ~SPoint();

  void Init(void);
  void Init(INTEGER d_, INTEGER nnz_); // This is NOT a zero vector
  void ReleaseAllMemory(void);
  void DeepCopy(const SPoint &G);

  //-------------------- Utilities --------------------
  //
  // x is the self data point.

  // Get dimension
  INTEGER GetD(void) const;

  // Get nnz
  INTEGER GetNNZ(void) const;

  // Get the pointer to idx
  INTEGER* GetPointerIdx(void) const;

  // Get the pointer to x
  double* GetPointerX(void) const;

  // Print the data point
  void PrintPoint(const char *name) const;

  //-------------------- Computations --------------------
  //
  // x is the self data point.

  // x'*y
  double InProd(const DPoint &y) const;
  double InProd(const SPoint &y) const;

  // sum_i |x_i-y_i|
  double Dist1(const DPoint &y) const;
  double Dist1(const SPoint &y) const;

  // sum_i |x_i-y_i|/|sigma_i|
  double Dist1(const DPoint &y, double *sigma) const;
  double Dist1(const SPoint &y, double *sigma) const;

  // sum_i |x_i-y_i|^2
  double Dist2(const DPoint &y) const;
  double Dist2(const SPoint &y) const;

  // sum_i (|x_i-y_i|/|sigma_i|)^2
  double Dist2(const DPoint &y, double *sigma) const;
  double Dist2(const SPoint &y, double *sigma) const;

  // sum_i 2 * x_i * y_i / (x_i + y_i)
  double KernelFuncChi2(const DPoint &y) const;
  double KernelFuncChi2(const SPoint &y) const;

protected:

private:

  INTEGER d;    // Dimension
  INTEGER nnz;  // Number of nonzeros
  INTEGER *idx; // Indices (Must be sorted increasingly)
  double *x;    // Values

  // The points are organized in a sparse format. For example, if d =
  // 10, idx = {0, 3, 9}, and x = {1.0, 2.1, 6.3}; then the point has
  // a dimension 10, nnz = 3, x[0] = 1.0, x[3] = 2.1, x[9] = 6.3, and
  // the rest of x is 0.

};

#endif
