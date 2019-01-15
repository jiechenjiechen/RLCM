#include "LibSVM_IO_binary.hpp"
namespace LibSVM_IO_binary {


//--------------------------------------------------------------------------
bool ReadData(const char *path_and_filename, DPointArray &X, DVector &y, INTEGER d){

  // Open file
  FILE *fp = NULL;
  fp = fopen(path_and_filename, "rb");
  if (!fp) {
    printf("ReadData. Error: Cannot open file %s. Function call takes no effect.\n", path_and_filename);
    return false;
  }

  // Floating point number in single precision
  INTEGER size_of_real = 4;

  // Obtain file size
  fseek(fp , 0 , SEEK_END);
  long file_size = ftell(fp);
  rewind(fp);
  INTEGER N = file_size/(d+1)/size_of_real;
  if (file_size != N*(d+1)*size_of_real) {
    printf("ReadData. Error: File size is not a multiple of d+1. Function call takes no effect.\n");
    fclose(fp);
    return false;
  }

  // Initialize X and y
  X.Init(N,d);
  y.Init(N);
  double *mX = X.GetPointer();
  double *my = y.GetPointer();

  // Create buffer
  char *buf = NULL;
  INTEGER buf_size = (d+1)*size_of_real;
  New_1D_Array<char, INTEGER>(&buf, buf_size);

  // Endianness
  bool big_endian = IsBigEndian();

  // Read in file row by row
  for (INTEGER i = 0; i < N; i++) {

    INTEGER nread = fread(buf, 1, buf_size, fp);
    if (nread != buf_size) {
      printf("ReadData. Error occurs when calling fread. Function call takes no effect.\n");
      break;
    }

    float val;
    if (big_endian) {
      memcpy(&val, buf, 4);
    }
    else {
      Swap4Bytes((char *)(&val), buf);
    }
    my[i] = val;
#ifdef USE_OPENMP
#pragma omp parallel for private(val)
#endif
    for (INTEGER j = 0; j < d; j++) {
      if (big_endian) {
        memcpy(&val, buf+(j+1)*size_of_real, 4);
      }
      else {
        Swap4Bytes((char *)(&val), buf+(j+1)*size_of_real);
      }
      mX[i+N*j] = val;
    }

  }

  // Clean up
  fclose(fp);
  Delete_1D_Array<char>(&buf);
  return true;

}


} // namespace LibSVM_IO_binary

