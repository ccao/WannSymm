#include <stdio.h>
#include <complex.h>
#include <string.h>
#include <math.h>

#include "constants.h"
#include "wannorb.h"

#include "mkl.h"
#include <mpi.h>


void generate_Lmatrix(dcomplex * Lmat, int l);
void rotate_Ylm(dcomplex * rot, int l, double axis[3], double alpha, int inv);
void generate_C2Ylm(dcomplex * C2Ylm, dcomplex * Ylm2C, int l);
void rotate_cubic(dcomplex * rot, int l, double axis[3], double alpha, int inv);

void generate_CS2JM(dcomplex * CS2JM, dcomplex * JM2CS, int l);
void generate_YlmS2JM(dcomplex * YlmS2JM, dcomplex * JM2YlmS, int l);
void print_orbital_definition_for_JM(char * seed, wannorb_JM * orb_info, int norb);