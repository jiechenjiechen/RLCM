// The SPointArray class implements a set of sparse data points; see
// the SPoint class for an individual sparse data point. The point set
// is treated as a sparse matrix, which is represented by three
// numbers and three arrays. Let the matrix have a size N*d with nnz
// nonzero elements, where N is the number of points and d is the
// dimension of the points. The three arrays are start, idx, and X,
// which means respectively, the starting location of a new point, the
// index of the nonzero elements, and the values of the nonzero
// elements. Such a storage format is consistent with the compressed
// sparse row (CSR) matrix format. For more details, see the class
// definition below.
//
// Note that the sparse matrix is effectively in row-major order. This
// is different from the storage of a dense point array (see
// DPointArray).
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

#ifndef _SPOINT_ARRAY_
#define _SPOINT_ARRAY_

#include "SPoint.hpp"
#include "DPointArray.hpp"

class SPointArray {

public:

  SPointArray();
  SPointArray(INTEGER N_, INTEGER d_, INTEGER nnz_);
  SPointArray(const SPointArray& G);
  SPointArray& operator= (const SPointArray& G);
  SPointArray& operator= (const DPointArray& G); // Convert from dense to sparse
  ~SPointArray();

  void Init(void);
  void Init(INTEGER N_, INTEGER d_, INTEGER nnz_); // This is NOT a zero matrix
  void ReleaseAllMemory(void);
  void DeepCopy(const SPointArray &G);

  //-------------------- Utilities --------------------
  //
  // X is the self data matrix.

  // Get dimension
  INTEGER GetD(void) const;

  // Get number of points
  INTEGER GetN(void) const;

  // Get nnz
  INTEGER GetNNZ(void) const;

  // Get the i-th point x
  void GetPoint(INTEGER i, SPoint &x) const;
  void GetPoint(INTEGER i, DPoint &x) const;

  // Get a consecutive chunk Y. Note that SPointArray is in row-major
  // order whereas DPointArray is in column-major order.
  void GetSubset(INTEGER istart, INTEGER n, SPointArray &Y) const;
  void GetSubset(INTEGER istart, INTEGER n, DPointArray &Y) const;

  // Get a subset Y. Note that SPointArray is in row-major order
  // whereas DPointArray is in column-major order.
  void GetSubset(INTEGER *iidx, INTEGER n, SPointArray &Y) const;
  void GetSubset(INTEGER *iidx, INTEGER n, DPointArray &Y) const;

  // Get the pointer to start
  INTEGER* GetPointerStart(void) const;

  // Get the pointer to idx
  INTEGER* GetPointerIdx(void) const;

  // Get the pointer to X
  double* GetPointerX(void) const;

  // Get the data matrix (stored in dense format)
  void AsDMatrix(DMatrix &A) const;

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
  void MatMat(const SPointArray &B, DMatrix &Y,
              MatrixMode ModeX, MatrixMode ModeB) const;

  // Bipartition and permute the points by using a hyperplane with
  // random orientation.
  //
  // This routine randomly partitions the data points (indexed from
  // istart to istart+n-1) into two equal halves, if n >= 2*N0. Let
  // the dividing hyperplane have a normal direction 'normal'. All the
  // points are projected along this direction and a median (named
  // 'offset') of the projected values is computed. The hyperplane
  // equation is thus h(x) = normal'*x - offset = 0. In the
  // nondegenerate case, there are m1 = floor(n/2) points x satisfying
  // h(x) < 0 and the rest satisfying h(x) > 0. The routine permutes
  // the points so that those indexed from istart to istart+m1-1
  // correspond to h(x) < 0 and the rest corresponds to h(x) > 0. The
  // return value of the routine is m1. Additionally, if the
  // permutation array perm is not NULL and is filled with ordering
  // information, then on exit it records the new order of the points.
  INTEGER RandomBipartition(INTEGER istart, INTEGER n, INTEGER N0,
                            INTEGER *perm, DPoint &normal, double &offset);

  // Similar to RandomBipartition, except that the orientation of the
  // hyperplane is the principal direction of the points.
  INTEGER PCABipartition(INTEGER istart, INTEGER n, INTEGER N0,
                         INTEGER *perm, DPoint &normal, double &offset);

  // Partitioning the points along the longest dimension of the
  // bounding box.
  //
  // This functionality is not implemented, because for sparse points,
  // such a partitioning may be highly imbalanced due to the excessive
  // number of zeros. One should not use this a partitioning
  // method. These functions are declared only as a placeholder.
  void ComputeBBox(double **bbox, INTEGER &dim);
  INTEGER BBoxBipartition(INTEGER start, INTEGER n, INTEGER N0,
                          INTEGER *perm, DPoint &normal, double &offset,
                          const double *bbox, INTEGER &which_dim);

protected:

private:

  INTEGER N;      // Number of points
  INTEGER d;      // Dimension
  INTEGER nnz;    // Total number of nonzeros
  INTEGER *start; //
  INTEGER *idx;   //
  double *X;      //

  // The three arrays 'start', 'idx', and 'X', altogether, are used to
  // represent points in a manner consistent with the CSR format of a
  // sparse matrix.
  //
  // 'start' has a length N+1, where start[i] means the location of
  // idx (and of X) where the i-th point begins, for i = 0:N-1. For
  // convenience, start[N] = nnz.
  //
  // 'idx' has a length nnz, where the segment idx[start[i]] to
  // idx[start[i+1]-1] stores the coordinates of the nonzero elements
  // of the i-th point. All 'idx' entries must be < d.
  //
  // 'X' has a length nnz, where the segment X[start[i]] to
  // X[start[i+1]-1] stores the values of the nonzero elements of the
  // i-th point.

  // Called by RandomBipartition() and PCABipartition()
  INTEGER Bipartition(INTEGER istart, INTEGER n,
                      INTEGER *perm, DPoint &normal, double &offset);

};

#endif
