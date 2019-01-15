// This namespace handles IO for ascii files that store a graph in the
// METIS format.
//
// The implementation is NOT parallelized.

#ifndef _METIS_IO_
#define _METIS_IO_

#include "../Matrices/SMatrix.hpp"

namespace METIS_IO {

// Read a data file in the METIS graph format. The graph must be an
// undirected, weighted graph and it is stored as a symmetric sparse
// matrix A.
//
// The matrix is symmetric positive definite, which means that the
// graph contains self-loops. For example, vertex 1 must have itself
// as a neighbor, with a certain positive weight. Let n be the number
// of vertices (also the number of rows/columns of A) and let A has
// nnz nonzeros. The number of edges is hence (nnz+n)/2.
//
// Note that something funny about METIS is that in the graph file,
// index starts from 1 (see e.g. page 11, Figure 2 of METIS mannual),
// whereas in the API, index starts from 0 (see e.g., page 24, Figure
// 3 of METIS mannual).
//
// Also note that the METIS graph format does not require the
// neighbors of a vertex to be sorted; but the sparse matrix class
// SMatrix requires the columns indices to be sorted.

bool ReadData(const char *path_and_filename, SMatrix &A);

} // namespace METIS_IO

#endif
