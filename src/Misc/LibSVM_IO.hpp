// This namespace handles IO for ascii files that store a data set in
// the LibSVM format.
//
// The implementation is NOT parallelized.

#ifndef _LIBSVM_IO_
#define _LIBSVM_IO_

#include "../Matrices/DPointArray.hpp"
#include "../Matrices/SPointArray.hpp"

namespace LibSVM_IO {

// Read a data file in the LibSVM format. The points are stored in X
// and the labels are stored in y. The size of X is N*d, where N
// (number of points) is inferred from the data file, and d (dimension
// of points) is a mandatoray input. One may choose to store the data
// points in a dense format (DPointArray) or in a sparse format
// (SPointArray).
//
// Note that LibSVM index starts from 1, whereas label index starts
// from 0.

bool ReadData(const char *path_and_filename, DPointArray &X, DVector &y, INTEGER d);
bool ReadData(const char *path_and_filename, SPointArray &X, DVector &y, INTEGER d);

// Subroutine called by ReadData(). Verify the file format and count N
// and nnz.
bool VerifyFileFormat(const char *path_and_filename, INTEGER d, INTEGER &N, INTEGER &nnz);

} // namespace LibSVM_IO

#endif
