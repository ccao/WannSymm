#ifdef CONSTANTS_H
#else
#define CONSTANTS_H

#define MAXLEN    512
#define MEDIUMLEN 256
#define SHORTLEN  128
#define MAXN   10
#define eps4   1E-4
#define eps5   1E-5
#define eps6   1E-6
#define eps7   1E-7
#define eps8  1E-8
#define sqrt2 1.41421356237309504880168872421

#define PI    3.14159265358979323846264338328
#define cmplx_i _Complex_I
#define MAX_NUM_of_SYMM  1536
#define MAX_NUM_of_atoms 1024
#define MAX_L 3

#include <complex.h>
typedef double complex dcomplex;
#define MKL_Complex16 dcomplex  // this tells MKL about user's MKL_Complex16 type

#endif
