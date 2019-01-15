#include "Raw_binary.hpp"
namespace Raw_binary {


//--------------------------------------------------------------------------
bool ReadData(const char *path_and_filename, DVector &y){

  // Open file
  FILE *fp = NULL;
  fp = fopen(path_and_filename, "rb");
  if (!fp) {
    printf("ReadData. Error: Cannot open file %s. Function call takes no effect.\n", path_and_filename);
    return false;
  }

  // Floating point number in double precision
  INTEGER size_of_real = 8;

  // Obtain file size
  fseek(fp , 0 , SEEK_END);
  long file_size = ftell(fp);
  rewind(fp);
  INTEGER N = file_size/size_of_real;
  if (file_size != N*size_of_real) {
    printf("ReadData. Error: File size is not a multiple of 8. Function call takes no effect.\n");
    fclose(fp);
    return false;
  }

  // Initialize y
  y.Init(N);
  double *my = y.GetPointer();

  // Create buffer
  char *buf = NULL;
  INTEGER batch_size = 1000; // Hardcoded
  INTEGER buf_size = batch_size * size_of_real;
  New_1D_Array<char, INTEGER>(&buf, buf_size);

  // Endianness
  bool big_endian = IsBigEndian();

  // Read in file batch by batch
  INTEGER num_batch = (INTEGER)ceil((double)N/batch_size);
  INTEGER last_batch_size = N - (num_batch-1) * batch_size;
  for (INTEGER i = 0; i < num_batch; i++) {

    INTEGER nread = fread(buf, size_of_real, batch_size, fp);
    if ( (i < num_batch-1 && nread != batch_size) ||
         (i == num_batch-1 && nread != last_batch_size) ) {
      printf("ReadData. Error occurs when calling fread. Function call takes no effect.\n");
      break;
    }

    double val;
#ifdef USE_OPENMP
#pragma omp parallel for private(val)
#endif
    for (INTEGER j = 0; j < nread; j++) {
      if (big_endian) {
        memcpy(&val, buf+j*size_of_real, size_of_real);
      }
      else {
        Swap8Bytes((char *)(&val), buf+j*size_of_real);
      }
      my[i*batch_size+j] = val;
    }

  }

  // Clean up
  fclose(fp);
  Delete_1D_Array<char>(&buf);
  return true;

}


} // namespace Raw_binary

