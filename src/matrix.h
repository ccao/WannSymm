#ifndef WSMATRIX_H
#define WSMATRIX_H

#include <stdio.h>
#include <stdlib.h>

#include "constants.h"
#include "mkl.h"

void matrix3x3_transpose(double out[3][3], double in[3][3]);
void matrix3x3_inverse(double out[3][3], double in[3][3]);
void matrix3x3_dot(double out[3][3], double a[3][3], double b[3][3]);
double matrix3x3_determinant(double m[3][3]);

void matrix3x3_copy(double out[3][3], double in[3][3]);
void matrix3x1_copy(double out[3], double in[3]);

#endif
