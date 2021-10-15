#include <stdio.h>
#include <complex.h>

#include "constants.h"

#include <mkl.h>


void rotate_spinor(dcomplex * rot, double axis[3], double alpha, int inv);
