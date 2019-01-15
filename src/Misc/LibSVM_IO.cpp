#include "LibSVM_IO.hpp"
namespace LibSVM_IO {


//--------------------------------------------------------------------------
bool ReadData(const char *path_and_filename, DPointArray &X, DVector &y, INTEGER d){

  INTEGER N, nnz; // nnz is not used
  bool correct = VerifyFileFormat(path_and_filename, d, N, nnz);
  if (!correct) {
    return false;
  }

  FILE *fp = NULL;
  fp = fopen(path_and_filename, "r");

  // Initilize X and y
  X.Init(N,d);
  y.Init(N);

  // Read the file again and populate X and y
  INTEGER i = 0;
  INTEGER numi = 0;
  double numf = 0.0;
  double *mX = X.GetPointer();
  double *my = y.GetPointer();
  char *line = NULL, *str = NULL, *saveptr = NULL, *subtoken = NULL;
  ssize_t read = 0;
  size_t len = 0;
  while ((read = getline(&line, &len, fp)) != -1) { // Read in a line
    INTEGER cnt = 0;
    for (str = line; ; str = NULL) {
      subtoken = strtok_r(str, ": ", &saveptr); // Tokenize the line
      if (subtoken == NULL) {
        break;
      }
      else {
        cnt++;
      }
      if (cnt == 1) {
        my[i] = atof(subtoken);
      }
      else if (cnt%2 == 0) {
        numi = String2Integer(subtoken);
      }
      else {
        numf = atof(subtoken);
        mX[i+N*(numi-1)] = numf; // LibSVM index starts from 1
      }
    }
    i++;
  }
  if (line) {
    free(line);
  }

  fclose(fp);
  return true;

}


//--------------------------------------------------------------------------
bool ReadData(const char *path_and_filename, SPointArray &X, DVector &y, INTEGER d){

  INTEGER N, nnz;
  bool correct = VerifyFileFormat(path_and_filename, d, N, nnz);
  if (!correct) {
    return false;
  }

  FILE *fp = NULL;
  fp = fopen(path_and_filename, "r");

  // Initilize X and y
  X.Init(N,d,nnz);
  y.Init(N);

  // Read the file again and popular X and y
  INTEGER i = 0, id = 0;
  INTEGER numi = 0;
  double numf = 0.0;
  INTEGER *mstart = X.GetPointerStart();
  INTEGER *midx = X.GetPointerIdx();
  double *mX = X.GetPointerX();
  double *my = y.GetPointer();
  char *line = NULL, *str = NULL, *saveptr = NULL, *subtoken = NULL;
  ssize_t read = 0;
  size_t len = 0;
  while ((read = getline(&line, &len, fp)) != -1) { // Read in a line
    INTEGER cnt = 0;
    for (str = line; ; str = NULL) {
      subtoken = strtok_r(str, ": ", &saveptr); // Tokenize the line
      if (subtoken == NULL) {
        break;
      }
      else {
        cnt++;
      }
      if (cnt == 1) {
        my[i] = atof(subtoken);
        mstart[i] = id;
      }
      else if (cnt%2 == 0) {
        numi = String2Integer(subtoken);
      }
      else {
        numf = atof(subtoken);
        midx[id] = numi-1; // LibSVM index starts from 1
        mX[id] = numf;
        id++;
      }
    }
    i++;
  }
  mstart[i] = id;
  if (line) {
    free(line);
  }

  fclose(fp);
  return true;

}


//--------------------------------------------------------------------------
bool VerifyFileFormat(const char *path_and_filename, INTEGER d, INTEGER &N, INTEGER &nnz){

  FILE *fp = NULL;
  fp = fopen(path_and_filename, "r");
  if (!fp) {
    printf("VerifyFileFormat. Error: Cannot open file %s. Function call takes no effect.\n", path_and_filename);
    return false;
  }

  N = 0;
  nnz = 0;
  char *line = NULL, *str = NULL, *saveptr = NULL, *subtoken = NULL;
  ssize_t read = 0;
  size_t len = 0;
  while ((read = getline(&line, &len, fp)) != -1) { // Read in a line
    N++;
    INTEGER cnt = 0, maxdim = 0;
    for (str = line; ; str = NULL) {
      subtoken = strtok_r(str, ": ", &saveptr); // Tokenize the line
      if (subtoken == NULL) {
        break;
      }
      else {
        cnt++;
      }
      if (cnt%2 == 0) { // Verify format and increase nnz
        INTEGER num = String2Integer(subtoken);
        if (maxdim > num) {
          printf("VerifyFileFormat. Error: Indices in line %ld are not in the ascending order. Stop reading data. Function call takes no effect.\n", (long)N);
          fclose(fp);
          return false;
        }
        else {
          maxdim = num;
        }
        nnz++;
      }
    }
    if (cnt%2 != 1) { // Verify format
      printf("VerifyFileFormat. Error: Line %ld does not conform with a LibSVM format. Stop reading data. Function call takes no effect.\n", (long)N);
      fclose(fp);
      return false;
    }
    else if (d < maxdim) {
      printf("VerifyFileFormat. Error: Line %ld indicates a point of dimension larger than d. Stop reading data. Function call takes no effect.\n", (long)N);
      fclose(fp);
      return false;
    }
  }
  if (line) {
    free(line);
  }

  fclose(fp);
  if (N == 0) {
    printf("VerifyFileFormat. Error: Empty file!\n");
    return false;
  }
  return true;

}


} // namespace LibSVM_IO

