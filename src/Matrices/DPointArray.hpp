// The DPointArray class implements a set of dense data points; see
// the DPoint class for an individual dense data point. The point set
// is treated as a dense matrix of size N*d, where N is the number of
// points and d is the dimension of the points.
//
// Note that the dense matrix is in column-major order. This is
// different from the storage of a sparse point array (see
// SPointArray).
//
// The implementation of some methods in this class is parallelized
// whereas others not. The general rule of thumb is that if the
// computation is w.r.t. one point (e.g., GetPoint), then it is not
// parallelized, because such computation can possibly be iterated
// through all points (where parallelism better takes place). If the
// computation is w.r.t. the whole point array (e.g.,
// RandomBipartition), then it is parallelized. The computation
// w.r.t. a subset of points (e.g., GetSubset) lies between the above
// two cases; for the moment this computation is parallelized.

#ifndef _DPOINT_ARRAY_
#define _DPOINT_ARRAY_

#include "DPoint.hpp"
#include "DMatrix.hpp"

class DPointArray {

public:

  DPointArray();
  DPointArray(INTEGER N_, INTEGER d_);
  DPointArray(const DPointArray& G);
  DPointArray& operator= (const DPointArray& G);
  ~DPointArray();

  void Init(void);
  void Init(INTEGER N_, INTEGER d_); // Initialized as a zero matrix
  void ReleaseAllMemory(void);
  void DeepCopy(const DPointArray &G);

  //-------------------- Utilities --------------------
  //
  // X is the self data matrix.

  // Get dimension
  INTEGER GetD(void) const;

  // Get number of points
  INTEGER GetN(void) const;

  // Get the i-th point x
  void GetPoint(INTEGER i, DPoint &x) const;

  // Get a consecutive chunk Y
  void GetSubset(INTEGER start, INTEGER n, DPointArray &Y) const;

  // Get a subset Y
  void GetSubset(INTEGER *idx, INTEGER n, DPointArray &Y) const;

  // Get the pointer
  double* GetPointer(void) const;

  // Get the data matrix
  void AsDMatrix(DMatrix &A) const;

  // Generate points with uniformly random coordinates inside the unit box
  void SetUniformRandom01(void);

  // Generate points uniformly random on the unit sphere
  void SetUniformSphere(void);

  // Similar to Matlab ndgrid. The first dimension varies the fastest
  void SetRegularGrid(INTEGER d_, const INTEGER *Dim,
                      const double *Lower, const double *Upper);

  // Print the point set
  void PrintPointArray(const char *name) const;

  //-------------------- Computations --------------------
  //
  // X is the self data matrix.

  // Center of the point set
  void Center(DPoint &c) const;

  // Root mean squared distances between the center and the points
  double StdvDist(void) const;

  // max_i x_i'*x_i
  double MaxPointNorm2(void) const;

  // y = mode(X) * b (size of X is N*d)
  void MatVec(const DVector &b, DVector &y, MatrixMode ModeX) const;
  void MatVec(const DPoint &b, DVector &y, MatrixMode ModeX) const;
  // Y = mode(X) * mode(B) (size of X is N*d)
  void MatMat(const DMatrix &B, DMatrix &Y,
              MatrixMode ModeX, MatrixMode ModeB) const;
  void MatMat(const DPointArray &B, DMatrix &Y,
              MatrixMode ModeX, MatrixMode ModeB) const;

  // Bipartition and permute the points by using a hyperplane with
  // random orientation.
  //
  // This routine randomly partitions the data points (indexed from
  // start to start+n-1) into two equal halves, if n >= 2*N0. Let the
  // dividing hyperplane have a normal direction 'normal'. All the
  // points are projected along this direction and a median (named
  // 'offset') of the projected values is computed. The hyperplane
  // equation is thus h(x) = normal'*x - offset = 0. In the
  // nondegenerate case, there are m1 = floor(n/2) points x satisfying
  // h(x) < 0 and the rest satisfying h(x) > 0. The routine permutes
  // the points so that those indexed from start to start+m1-1
  // correspond to h(x) < 0 and the rest corresponds to h(x) > 0. The
  // return value of the routine is m1. Additionally, if the
  // permutation array perm is not NULL and is filled with ordering
  // information, then on exit it records the new order of the points.
  INTEGER RandomBipartition(INTEGER start, INTEGER n, INTEGER N0,
                            INTEGER *perm, DPoint &normal, double &offset);

  // Similar to RandomBipartition, except that the orientation of the
  // hyperplane is the principal direction of the points.
  INTEGER PCABipartition(INTEGER start, INTEGER n, INTEGER N0,
                         INTEGER *perm, DPoint &normal, double &offset);

  // Partitioning the points along the longest dimension of the
  // bounding box.
  void ComputeBBox(double **bbox, INTEGER &dim);
  INTEGER BBoxBipartition(INTEGER start, INTEGER n, INTEGER N0,
                          INTEGER *perm, DPoint &normal, double &offset,
                          const double *bbox, INTEGER &which_dim);

protected:

private:

  INTEGER N; // Number of points
  INTEGER d; // Dimension
  double *X; // Points

  // X is an array of length N*d, where X[0] to X[N-1] give the first
  // coordinate value of the points, X[N] to X[2N-1] give the second
  // coordinate values, etc.

  // Called by RandomBipartition() and PCABipartition()
  INTEGER Bipartition(INTEGER start, INTEGER n,
                      INTEGER *perm, const DPoint &normal, double &offset);

};

#endif
