#include <stdio.h>
#include <complex.h>
#include <string.h>
#include <math.h>

#include "constants.h"

#include "mkl.h"
#include <mpi.h>


void generate_Lmatrix(dcomplex * Lmat, int l);
void rotate_Ylm(dcomplex * rot, int l, double axis[3], double alpha, int inv);
void generate_C2Ylm(dcomplex * C2Ylm, dcomplex * Ylm2C, int l);
void rotate_cubic(dcomplex * rot, int l, double axis[3], double alpha, int inv);
