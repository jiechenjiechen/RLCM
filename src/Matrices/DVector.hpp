// The DVector class implements a dense vector in double precision. A
// substantial portion of the methods in this class are wrappers of
// BLAS1 routines. The interface tends to be coherent with Matlab
// constructs.
//
// The implementation of this class is parallelized.

#ifndef _DVECTOR_
#define _DVECTOR_

#include "../Misc/Common.hpp"

class DVector {

public:

  DVector();
  DVector(INTEGER N_);
  DVector(const DVector &G);
  DVector& operator= (const DVector &G);
  ~DVector();

  void Init(void);
  void Init(INTEGER N_); // Initialized as a zero vector
  void ReleaseAllMemory(void);
  void DeepCopy(const DVector &G);

  //-------------------- Utilities --------------------
  //
  // a is the self vector.

  // Get length
  INTEGER GetN(void) const;

  // Get a(i)
  double GetEntry(INTEGER i) const;
  // b = a(RowStart:RowStart+nRow-1)
  void GetBlock(INTEGER RowStart, INTEGER nRow, DVector &b) const;
  // b = a(RowStart:RowStart+nRow-1). Memory of b must be preallocated
  void GetBlock(INTEGER RowStart, INTEGER nRow, double *b) const;
  // b = a(idx). n is the length of idx
  void GetBlock(INTEGER *idx, INTEGER n, DVector &b) const;
  // Get the double* pointer
  double* GetPointer(void) const;

  // Set a(i) = b
  void SetEntry(INTEGER i, double b);
  // a(RowStart:RowStart+nRow-1) = b
  void SetBlock(INTEGER RowStart, INTEGER nRow, const DVector &b);
  void SetBlock(INTEGER RowStart, INTEGER nRow, const double *b);

  // a = c
  void SetConstVal(double c);
  // a = rand()
  void SetUniformRandom01(void);
  // a = randn()
  void SetStandardNormal(void);
  // a = trnd(1); each element is a student-t of degree 1
  void SetStudentT1(void);
  // The whole vector is a multivariate student-t of degree 1
  void SetMultivariateStudentT1(void);
  // a = random('sech')
  void SetRandomSech(void);

  // Permutation
  // a_new(i) = a_old(Perm(i))
  void Permute(const INTEGER *Perm, INTEGER N_);
  // a_new(iPerm(i)) = a_old(i)
  void iPermute(const INTEGER *iPerm, INTEGER N_);
  // Note: If Perm is the permutation array and iPerm is the inverse
  // permutation array, then
  //    Permute(Perm, N) is equivalent to iPermute(iPerm, N).    (1)
  // Moreover,
  //    Permute(iPerm, N) is equivalent to iPermute(Perm, N),    (2)
  // and (1) has the opposite effect to that of (2).

  // Sorting
  // a = sort(a, type)
  void Sort(SortType type);
  // a = sort(a, type). Sort accoding to |a|
  void SortByMagnitude(SortType type);

  // Find elements strictly larger than tol. The array idx contains
  // the indices of the elements that satisfy the requirement. The
  // return value is the number of such elements. Memory of idx must
  // be preallocated. For safety, it is recommended that idx has the
  // same length as the vector itself. The returned indices are
  // ordered increasingly.
  INTEGER FindLargerThan(double tol, INTEGER *idx) const;

  // Print the vector in the Matlab form
  void PrintVectorMatlabForm(const char *name) const;

  //-------------------- Vector computations --------------------
  //
  // a is the self vector.

  // a = -a  or  b = -a
  void Negate(void);
  void Negate(DVector &b) const;

  // a = a + b  or  c = a + b
  void Add(double b);
  void Add(const DVector &b);
  void Add(double b, DVector &c) const;
  void Add(const DVector &b, DVector &c) const;

  // a = a - b  or  c = a - b
  void Subtract(double b);
  void Subtract(const DVector &b);
  void Subtract(double b, DVector &c) const;
  void Subtract(const DVector &b, DVector &c) const;

  // a = a * b  or  c = a * b
  void Multiply(double b);
  void Multiply(double b, DVector &c) const;
  // a = a .* b  or  c = a .* b
  void Multiply(const DVector &b);
  void Multiply(const DVector &b, DVector &c) const;

  // a = a / b  or  c = a / b
  void Divide(double b);
  void Divide(double b, DVector &c) const;
  // a = a ./ b  or  c = a ./ b
  void Divide(const DVector &b);
  void Divide(const DVector &b, DVector &c) const;

  // a' * b
  double InProd(const DVector &b) const;

  // a = abs(a)  or  b = abs(a)
  void Abs(void);
  void Abs(DVector &b) const;

  // a = sqrt(a)  or  b = sqrt(a)
  void Sqrt(void);
  void Sqrt(DVector &b) const;

  // a = a.^2  or  b = a.^2
  void Square(void);
  void Square(DVector &b) const;

  // a = 1./a  or  b = 1./a
  void Inv(void);
  void Inv(DVector &b) const;

  // norm(a,2)
  double Norm2(void) const;

  // norm(a,1)
  double Norm1(void) const;

  // min(a). idx is the first location of the min
  double Min(void) const;
  double Min(INTEGER &idx) const;

  // max(a). idx is the first location of the max
  double Max(void) const;
  double Max(INTEGER &idx) const;

  // sum(a)
  double Sum(void) const;

  // mean(a)
  double Mean(void) const;

  //-------------------- Build Response Vector --------------------
  //
  // Build a simulated response vector from test function. The class
  // TestFunction must have the following method:
  //
  //   double Eval(const Point &x) const;
  //
  // The class PointArray must have the following methods:
  //
  //   INTEGER GetN(void) const;
  //   void GetPoint(INTEGER i, Point &x) const;
  template<class TestFunction, class Point, class PointArray>
  void BuildResponseVector(const TestFunction &mTestFunction,
                           const PointArray &X);

protected:

private:

  INTEGER N; // Vector length
  double *a; // Vector data

};

#include "DVector.tpp"

#endif
