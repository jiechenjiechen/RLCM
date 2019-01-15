#include "Common.hpp"


//--------------------------------------------------------------------------
double Square(double x) {
  return x*x;
}


//--------------------------------------------------------------------------
double Diff1(const double *a, const double *b, INTEGER n) {
  double ret = 0.0;
#ifdef USE_OPENMP
#pragma omp parallel for reduction(+:ret)
#endif
  for (INTEGER i = 0; i < n; i++) {
    ret += fabs(a[i]-b[i]);
  }
  return ret;
}


//--------------------------------------------------------------------------
double Diff1(const INTEGER *a, const INTEGER *b, INTEGER n) {
  double ret = 0.0;
#ifdef USE_OPENMP
#pragma omp parallel for reduction(+:ret)
#endif
  for (INTEGER i = 0; i < n; i++) {
    ret += fabs((double)a[i]-b[i]);
  }
  return ret;
}


//--------------------------------------------------------------------------
bool IsBigEndian(void) {
  INTEGER num = 1;
  if(*(char *)&num == 1) {
    return false;
  }
  else {
    return true;
  }
}


//--------------------------------------------------------------------------
void Swap4Bytes(char *desc, char *src) {
  desc[0] = src[3];
  desc[1] = src[2];
  desc[2] = src[1];
  desc[3] = src[0];
}


//--------------------------------------------------------------------------
void Swap8Bytes(char *desc, char *src) {
  desc[0] = src[7];
  desc[1] = src[6];
  desc[2] = src[5];
  desc[3] = src[4];
  desc[4] = src[3];
  desc[5] = src[2];
  desc[6] = src[1];
  desc[7] = src[0];
}


//--------------------------------------------------------------------------
void SwapNBytes(char *desc, char *src, INTEGER N) {
  for (INTEGER i = 0; i < N; i++) {
    desc[i] = src[N-i-1];
  }
}


//--------------------------------------------------------------------------
// 0: even; 1: odd
INTEGER Parity(INTEGER *ipiv, INTEGER n) {
  INTEGER i = 0, par = 0;
#ifdef USE_OPENMP
#pragma omp parallel for reduction(+:par)
#endif
  for (i = 0; i < n; i++) {
    if (ipiv[i]-1 != i) {
      par++;
    }
  }
  par %= 2;
  return par;
}


//--------------------------------------------------------------------------
void UniformRandom01(double *a, INTEGER n) {
  for (INTEGER j = 0; j < n; j++) {
    a[j] = (double)random()/RAND_MAX;
  }
}


//--------------------------------------------------------------------------
// Box-Muller: U, V from uniform[0,1]
// X = sqrt(-2.0*log(U))*cos(2.0*M_PI*V)
// Y = sqrt(-2.0*log(U))*sin(2.0*M_PI*V)
void StandardNormal(double *a, INTEGER n) {
  for (INTEGER i = 0; i < n/2; i++) {
    double U = (double)random()/RAND_MAX;
    double V = (double)random()/RAND_MAX;
    double common1 = sqrt(-2.0*log(U));
    double common2 = PIx2*V;
    a[2*i] = common1 * cos(common2);
    a[2*i+1] = common1 * sin(common2);
  }
  if (n%2 == 1) {
    double U = (double)random()/RAND_MAX;
    double V = (double)random()/RAND_MAX;
    double common1 = sqrt(-2.0*log(U));
    double common2 = PIx2*V;
    a[n-1] = common1 * cos(common2);
  }
}


//--------------------------------------------------------------------------
// Let X and Y be independent standard normal. Then X/fabs(Y) follows
// student t with 1 degree of freedom.
void StudentT1(double *a, INTEGER n) {
  for (INTEGER i = 0; i < n; i++) {
    double V = (double)random()/RAND_MAX;
    a[i] = tan(PIx2*V); // Why is cot not in math.h??? :-(
    if (V > 0.5) {
      a[i] = -a[i];
    }
  }
}


//--------------------------------------------------------------------------
// Let X and Y be independent standard normal. Then X/fabs(Y) follows
// student t with 1 degree of freedom.
void MultivariateStudentT1(double *a, INTEGER n) {
  StandardNormal(a, n);
  double b = 0.0;
  StandardNormal(&b, 1);
  b = fabs(b);
  for (INTEGER i = 0; i < n; i++) {
    a[i] /= b;
  }
}


//--------------------------------------------------------------------------
// Inverse transform sampling: X from uniform[0,1].
// H = 2/pi * ln( tan(pi/2 * X) )
void RandomSech(double *a, INTEGER n) {
  for (INTEGER i = 0; i < n; i++) {
    double X = (double)random()/RAND_MAX;
    a[i] = M_2_PI * log(tan(M_PI_2 * X));
  }
}


//--------------------------------------------------------------------------
// Knuth shuffle
void RandPerm(INTEGER n, INTEGER k, INTEGER *a) {
  INTEGER *b = NULL;
  New_1D_Array<INTEGER, INTEGER>(&b, n);
  INTEGER i = 0;
  for (i = 0; i < n; i++) {
    b[i] = i;
  }
  for (i = 0; i < k; i++) {
    Swap<INTEGER>(b[i], b[random()%(n-i)+i]);
  }
  memcpy(a, b, k*sizeof(INTEGER));
  Delete_1D_Array<INTEGER>(&b);
}


//--------------------------------------------------------------------------
int CompareNaturalOrderLess(const void *x, const void *y) {
  if      ((*(double*)x) <  (*(double*)y))   { return -1; }
  else if ((*(double*)x) == (*(double*)y))   { return 0; }
  else  /*((*(double*)x) >  (*(double*)y))*/ { return 1; }
}


//--------------------------------------------------------------------------
int CompareNaturalOrderGreater(const void *x, const void *y) {
  if      ((*(double*)x) >  (*(double*)y))   { return -1; }
  else if ((*(double*)x) == (*(double*)y))   { return 0; }
  else  /*((*(double*)x) <  (*(double*)y))*/ { return 1; }
}


//--------------------------------------------------------------------------
int CompareIntegerNaturalOrderLess(const void *x, const void *y) {
  if      ((*(INTEGER*)x) <  (*(INTEGER*)y))   { return -1; }
  else if ((*(INTEGER*)x) == (*(INTEGER*)y))   { return 0; }
  else  /*((*(INTEGER*)x) >  (*(INTEGER*)y))*/ { return 1; }
}


//--------------------------------------------------------------------------
int CompareIntegerNaturalOrderGreater(const void *x, const void *y) {
  if      ((*(INTEGER*)x) >  (*(INTEGER*)y))   { return -1; }
  else if ((*(INTEGER*)x) == (*(INTEGER*)y))   { return 0; }
  else  /*((*(INTEGER*)x) <  (*(INTEGER*)y))*/ { return 1; }
}


//--------------------------------------------------------------------------
int CompareByMagnitudeLess(const void *x, const void *y) {
  if      (fabs(*(double*)x) <  fabs(*(double*)y))   { return -1; }
  else if (fabs(*(double*)x) == fabs(*(double*)y))   { return 0; }
  else  /*(fabs(*(double*)x) >  fabs(*(double*)y))*/ { return 1; }
}


//--------------------------------------------------------------------------
int CompareByMagnitudeGreater(const void *x, const void *y) {
  if      (fabs(*(double*)x) >  fabs(*(double*)y))   { return -1; }
  else if (fabs(*(double*)x) == fabs(*(double*)y))   { return 0; }
  else  /*(fabs(*(double*)x) <  fabs(*(double*)y))*/ { return 1; }
}


//--------------------------------------------------------------------------
LOGICAL SelectLeftHalfPlane(double *WR, double *WI) {
  if (*WR < 0) {
    return FTRUE;
  }
  else {
    return FFALSE;
  }
}


//--------------------------------------------------------------------------
LOGICAL SelectRightHalfPlane(double *WR, double *WI) {
  if (*WR > 0) {
    return FTRUE;
  }
  else {
    return FFALSE;
  }
}


//--------------------------------------------------------------------------
int CompareElm1Idx(const void *x, const void *y) {
  if      (((Elm1*)x)->idx <  ((Elm1*)y)->idx)   { return -1; }
  else if (((Elm1*)x)->idx == ((Elm1*)y)->idx)   { return 0; }
  else  /*(((Elm1*)x)->idx >  ((Elm1*)y)->idx)*/ { return 1; }
};


//--------------------------------------------------------------------------
int CompareElm2Idx(const void *x, const void *y) {
  if      (((Elm2*)x)->idx <  ((Elm2*)y)->idx)   { return -1; }
  else if (((Elm2*)x)->idx == ((Elm2*)y)->idx)   { return 0; }
  else  /*(((Elm2*)x)->idx >  ((Elm2*)y)->idx)*/ { return 1; }
};


//--------------------------------------------------------------------------
int CompareElm2Above(const void *x, const void *y) {
  if      (((Elm2*)x)->above <  ((Elm2*)y)->above)   { return -1; }
  else if (((Elm2*)x)->above == ((Elm2*)y)->above)   { return 0; }
  else  /*(((Elm2*)x)->above >  ((Elm2*)y)->above)*/ { return 1; }
};


//--------------------------------------------------------------------------
int CompareElm2Below(const void *x, const void *y) {
  if      (((Elm2*)x)->below <  ((Elm2*)y)->below)   { return -1; }
  else if (((Elm2*)x)->below == ((Elm2*)y)->below)   { return 0; }
  else  /*(((Elm2*)x)->below >  ((Elm2*)y)->below)*/ { return 1; }
};


//--------------------------------------------------------------------------
int CompareElm3Dim(const void *x, const void *y) {
  if      (((Elm3*)x)->dim <  ((Elm3*)y)->dim)   { return -1; }
  else if (((Elm3*)x)->dim == ((Elm3*)y)->dim)   { return 0; }
  else  /*(((Elm3*)x)->dim >  ((Elm3*)y)->dim)*/ { return 1; }
};


//--------------------------------------------------------------------------
#ifdef USE_MATERN
double besselk(double nu, double x) {
  double zr = x;
  double zi = 0;
  double fnu = fabs(nu);
  int kode = 1;
  int n = 1;
  double cyr[1];
  double cyi[1];
  int nz;
  int ierr;
  zbesk_(&zr, &zi, &fnu, &kode, &n, cyr, cyi, &nz, &ierr);
  if (ierr) {
    printf("Besselk error. nu = %g, x = %g, ierr = %d. cyr = %g, cyi = %g. Return cyr\n", nu, x, ierr, cyr[0], cyi[0]);
  }
  return cyr[0];
}
#endif

