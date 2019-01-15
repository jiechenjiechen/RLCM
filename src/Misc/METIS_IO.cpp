#include "METIS_IO.hpp"
namespace METIS_IO {


//--------------------------------------------------------------------------
bool ReadData(const char *path_and_filename, SMatrix &A) {

  FILE *fp = NULL;
  fp = fopen(path_and_filename, "r");
  if (!fp) {
    printf("ReadData. Error: Cannot open file %s. Function call takes no effect.\n", path_and_filename);
    return false;
  }

  // Read in nvtxs, nedges, and fmt. fmt must be 001
  ssize_t read = 0;
  char *line = NULL, *str = NULL, *saveptr = NULL, *subtoken = NULL;
  size_t len = 0;
  read = getline(&line, &len, fp);
  if (read == -1) {
    printf("ReadData. Error: Fail to read file header. Stop reading data. Function call takes no effect.\n");
    fclose(fp);
    return false;
  }
  INTEGER cnt = 0;
  INTEGER nvtxs = 0, nedges = 0, fmt = 0;
  for (str = line, cnt = 0; cnt < 3; str = NULL, cnt++) {
    subtoken = strtok_r(str, " ", &saveptr); // Tokenize the line
    if (subtoken == NULL) {
      printf("ReadData. Error: File header does not contain three numbers. Stop reading data. Function call takes no effect.\n");
      fclose(fp);
      return false;
    }
    if (cnt == 0) {
      nvtxs = String2Integer(subtoken);
    }
    else if (cnt == 1) {
      nedges = String2Integer(subtoken);
    }
    else /*if (cnt == 2)*/ {
      fmt = String2Integer(subtoken);
      if (fmt != 1) {
        printf("ReadData. Error: The third number in file header is not 001. Stop reading data. Function call takes no effect.\n");
        fclose(fp);
        return false;
      }
    }
  }

  // Initialize A
  INTEGER nnz = nedges*2-nvtxs;
  A.Init(nvtxs, nnz);
  INTEGER *start = A.GetPointerStart();
  INTEGER *idx = A.GetPointerIdx();
  double *mA = A.GetPointerA();

  // Read the rest of the file
  INTEGER cnt_nnz = 0;
  start[0] = 0;
  for (INTEGER i = 0; i < nvtxs; i++) {
    read = getline(&line, &len, fp);
    if (read == -1) {
      printf("ReadData. Error: Fail to read line %ld. Stop reading data. Function call takes no effect.\n", (long)i);
      fclose(fp);
      return false;
    }
    cnt = 0;
    for (str = line; ; str = NULL) {
      subtoken = strtok_r(str, " ", &saveptr); // Tokenize the line
      if (subtoken == NULL) {
        break;
      }
      else {
        cnt++;
      }
      if (cnt%2 == 1) {
        idx[cnt_nnz] = String2Integer(subtoken)-1;  // METIS index starts from 1
      }
      else {
        mA[cnt_nnz] = atof(subtoken);
        cnt_nnz++;
      }
    }
    start[i+1] = cnt_nnz;
    if (cnt%2 != 0) {
      printf("ReadData. Error: The number of numbers in line %ld is not even. Stop reading data. Function call takes no effect.\n", (long)i);
      fclose(fp);
      return false;
    }
  }
  if (cnt_nnz != nnz) {
    printf("ReadData. Error: The number of nonzeros (%ld) does not match the number of edges (%ld) and the number of vertices (%ld). Stop reading data. Function call takes no effect.\n", (long)cnt_nnz, (long)nedges, (long)nvtxs);
    fclose(fp);
    return false;
  }
  if (line) {
    free(line);
  }

  // Sort the neighbor lists
  Elm1 *elm = NULL;
  New_1D_Array<Elm1, INTEGER>(&elm, nvtxs);
  for (INTEGER j = 0; j < nvtxs; j++) {
    INTEGER istart = start[j];
    INTEGER n = start[j+1] - istart;
    for (INTEGER k = 0; k < n; k++) {
      elm[k] = { idx[istart+k], mA[istart+k] };
    }
    qsort(elm, n, sizeof(Elm1), CompareElm1Idx);
    for (INTEGER k = 0; k < n; k++) {
      idx[istart+k] = elm[k].idx;
      mA[istart+k] = elm[k].val;
    }
  }
  Delete_1D_Array<Elm1>(&elm);

  // Verify symmetry
  for (INTEGER i = 0; i < nvtxs; i++) {
    for (INTEGER j = start[i]; j < start[i+1]; j++) {
      INTEGER col_idx = idx[j];
      double val = mA[j];
      if (col_idx > i) {
        INTEGER *loc;
        loc = (INTEGER *)bsearch(&i, idx+start[col_idx],
                                 start[col_idx+1]-start[col_idx],
                                 sizeof(INTEGER),
                                 CompareIntegerNaturalOrderLess);
        if (loc == NULL) {
          printf("ReadData. Error: Matrix is not symmetric. A[%ld,%ld] = %g whereas A[%ld,%ld] = 0. Function call takes no effect.\n", (long)i, (long)col_idx, val, (long)col_idx, (long)i);
          fclose(fp);
          return false;
        }
        if (mA[loc-idx] != val) {
          printf("ReadData. Error: Matrix is not symmetric. A[%ld,%ld] = %g whereas A[%ld,%ld] = %g. Function call takes no effect.\n", (long)i, (long)col_idx, val, (long)col_idx, (long)i, mA[loc-idx]);
          fclose(fp);
          return false;
        }
      }
    }
  }

  fclose(fp);
  return true;

}


} // namespace METIS_IO
