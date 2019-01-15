// This namespace handles IO for ascii files that store a data set in
// a ad hoc LibSVM-like format.
//
// The implementation is parallelized.

#ifndef _LIBSVM_IO_BINARY_
#define _LIBSVM_IO_BINARY_

#include "../Matrices/DPointArray.hpp"

namespace LibSVM_IO_binary {

// Read a data file in a LibSVM-like format. Let the data set X have
// size N*d and the label vector y have size N*1. The file stores the
// matrix [y X] in binary format and in row-major order, where each
// number is a 4-byte, big endian, floating point number.

// This format is an ad hoc imitation of the LibSVM format, serving
// the main purpose of efficient file IO, while also brings in the
// benefit of a smaller file size. The matrix is treated fully dense
// with explicit zeros, if any. There is minimal error check so the
// user is responsible for ensuring that the data file is precisely
// formatted.

bool ReadData(const char *path_and_filename, DPointArray &X, DVector &y, INTEGER d);

} // namespace LibSVM_IO_binary

#endif
