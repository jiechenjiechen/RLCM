// This namespace handles IO for binary files that store a sparse
// matrix.
//
// The implementation is parallelized.

#ifndef _CSR_IO_BINARY_
#define _CSR_IO_BINARY_

#include "../Matrices/SMatrix.hpp"

namespace CSR_IO_binary {

// Read sparse matrix in the CSR format (data file is binary and big
// endian). The function can read in any square sparse matrix, not
// necessarily a symmetric positive definite one.
//
// File format is the following---
//
// 1. The first four bytes is an integer, either 4 or 8, indicating
//    the length (in bytes) of an integer.
//
// 2. What follow are three integers M, N, and nnz (see SMatrix.hpp).
//
// 3. Then, the integer array start, integer array index, and double
//    precision real array A (see SMatrix.hpp).
//
// This format is an ad hoc replacement of the METIS graph format,
// serving the main purpose of reducing the file size, while also
// brings in the benefit of more efficient file IO. There is minimal
// error check so the user is responsible for ensuring that the data
// file is precisely formatted.

bool ReadData(const char *path_and_filename, SMatrix &A);

} // namespace CSR_IO_binary

#endif
