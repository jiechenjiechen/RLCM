// This code is used to configure the Fortran code besselk.f so that
// it works properly for the CPU architecture. The original code can
// be found from http://www.netlib.org/blas/d1mach.f and
// http://www.netlib.org/blas/i1mach.f
//
// This code supports only 32-bit integers.

#include <math.h>
#include <stdio.h>
#include <float.h>
#include <limits.h>

int main(int argc, char **argv) {

  printf("      DATA DMACH(1) / %.17g /\n", DBL_MIN);
  printf("      DATA DMACH(2) / %.17g /\n", DBL_MAX);
  printf("      DATA DMACH(3) / %.17g /\n", DBL_EPSILON/FLT_RADIX);
  printf("      DATA DMACH(4) / %.17g /\n", DBL_EPSILON);
  printf("      DATA DMACH(5) / %.16g /\n", log10((double)FLT_RADIX));

  printf("      DATA IMACH( 1) / %d /\n", 5);
  printf("      DATA IMACH( 2) / %d /\n", 6);
  printf("      DATA IMACH( 3) / %d /\n", 7);
  printf("      DATA IMACH( 4) / %d /\n", 0);
  printf("      DATA IMACH( 5) / %d /\n", 32); // should = sizeof(int)*8
  printf("      DATA IMACH( 6) / %ld /\n", sizeof(int));
  printf("      DATA IMACH( 7) / %d /\n", 2);
  printf("      DATA IMACH( 8) / %d /\n", 31); // should = (#5)-1
  //printf("      DATA IMACH( 9) / %ld /\n", LONG_MAX); // should = (#7)^(#8)-1
  printf("      DATA IMACH( 9) / %d /\n", 2147483647); // should = (#7)^(#8)-1
  printf("      DATA IMACH(10) / %d /\n", FLT_RADIX);
  printf("      DATA IMACH(11) / %d /\n", FLT_MANT_DIG);
  printf("      DATA IMACH(12) / %d /\n", FLT_MIN_EXP);
  printf("      DATA IMACH(13) / %d /\n", FLT_MAX_EXP);
  printf("      DATA IMACH(14) / %d /\n", DBL_MANT_DIG);
  printf("      DATA IMACH(15) / %d /\n", DBL_MIN_EXP);
  printf("      DATA IMACH(16) / %d /\n", DBL_MAX_EXP);

  return 0;
}
