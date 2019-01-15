// This namespace handles IO for binary files that store a dense
// vector.
//
// The implementation is parallelized.

#ifndef _RAW_BINARY_
#define _RAW_BINARY_

#include "../Matrices/DVector.hpp"

namespace Raw_binary {

// Read in an array of real numbers in double precision. Each number
// is a 8-byte and big endian.

bool ReadData(const char *path_and_filename, DVector &y);

} // namespace Raw_binary

#endif
