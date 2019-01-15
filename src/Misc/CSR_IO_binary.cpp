#include "CSR_IO_binary.hpp"
namespace CSR_IO_binary {


//--------------------------------------------------------------------------
bool ReadData(const char *path_and_filename, SMatrix &A) {

  FILE *fp = NULL;
  fp = fopen(path_and_filename, "r");
  if (!fp) {
    printf("ReadData. Error: Cannot open file %s. Function call takes no effect.\n", path_and_filename);
    return false;
  }

  // Endianness
  bool big_endian = IsBigEndian();

  // First 4 bytes (size_of_int)
  char *buf = NULL;
  INTEGER read_size = 4;
  int size_of_int = 0;
  New_1D_Array<char, INTEGER>(&buf, read_size);
  INTEGER nread = fread(buf, 1, read_size, fp);
  if (nread != read_size) {
    printf("ReadData. Error occurs when calling fread. Function call takes no effect.\n");
    fclose(fp);
    return false;
  }
  if (big_endian) {
    memcpy(&size_of_int, buf, read_size);
  }
  else {
    Swap4Bytes((char *)(&size_of_int), buf);
  }
  if ((INTEGER)sizeof(INTEGER) < size_of_int) {
    printf("ReadData. Error: Input matrix is too large. The program must be compiled with 64bit integer. Function call takes no effect.\n");
    fclose(fp);
    return false;
  }
  Delete_1D_Array<char>(&buf);

  // Next 3 integers (M, N, nnz)
  INTEGER M = 0;
  INTEGER N = 0;
  INTEGER nnz = 0;
  New_1D_Array<char, INTEGER>(&buf, size_of_int);
  for (INTEGER i = 0; i < 3; i++) {
    nread = fread(buf, size_of_int, 1, fp);
    if (nread != 1) {
      printf("ReadData. Error occurs when calling fread. Function call takes no effect.\n");
      fclose(fp);
      return false;
    }
    char *ptr = NULL;
    if (i == 0) {
      ptr = (char *)(&M);
    }
    else if (i == 1) {
      ptr = (char *)(&N);
    }
    else {
      ptr = (char *)(&nnz);
    }
    if (big_endian) {
      memcpy(ptr, buf, size_of_int);
    }
    else {
      SwapNBytes(ptr, buf, size_of_int);
    }
  }
  Delete_1D_Array<char>(&buf);

  // Init
  A.Init(M, N, nnz);
  INTEGER *start = A.GetPointerStart();
  INTEGER *idx = A.GetPointerIdx();
  double *SA = A.GetPointerA();

  // Array start
  New_1D_Array<char, INTEGER>(&buf, (M+1)*size_of_int);
  nread = fread(buf, size_of_int, M+1, fp);
  if (nread != M+1) {
    printf("ReadData. Error occurs when calling fread. Function call takes no effect.\n");
    fclose(fp);
    return false;
  }
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < M+1; i++) {
    if (big_endian) {
      memcpy(start+i, buf+size_of_int*i, size_of_int);
    }
    else {
      SwapNBytes((char *)(start+i), buf+size_of_int*i, size_of_int);
    }
  }
  Delete_1D_Array<char>(&buf);

  // Array idx
  New_1D_Array<char, INTEGER>(&buf, nnz*size_of_int);
  nread = fread(buf, size_of_int, nnz, fp);
  if (nread != nnz) {
    printf("ReadData. Error occurs when calling fread. Function call takes no effect.\n");
    fclose(fp);
    return false;
  }
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < nnz; i++) {
    if (big_endian) {
      memcpy(idx+i, buf+size_of_int*i, size_of_int);
    }
    else {
      SwapNBytes((char *)(idx+i), buf+size_of_int*i, size_of_int);
    }
  }
  Delete_1D_Array<char>(&buf);

  // Array SA
  INTEGER size_of_double = 8;
  New_1D_Array<char, INTEGER>(&buf, nnz*size_of_double);
  nread = fread(buf, size_of_double, nnz, fp);
  if (nread != nnz) {
    printf("ReadData. Error occurs when calling fread. Function call takes no effect.\n");
    fclose(fp);
    return false;
  }
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (INTEGER i = 0; i < nnz; i++) {
    if (big_endian) {
      memcpy(SA+i, buf+size_of_double*i, size_of_double);
    }
    else {
      Swap8Bytes((char *)(SA+i), buf+size_of_double*i);
    }
  }
  Delete_1D_Array<char>(&buf);

  // Clean up
  fclose(fp);
  return true;

}


} // namespace CSR_IO_binary
