#ifdef BNDSTRUCT_H
#else
#define BNDSTRUCT_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>

#include "constants.h"
#include "vector.h"
#include "matrix.h"
#include "wanndata.h"
#include "readinput.h"
#include "rotate_basis.h"
#include "rotate_orbital.h"
#include "rotate_spinor.h"
#include "rotate_ham.h"

#include "mkl.h"


typedef struct __smallgroup{
    int   order;    // order of the small group
    int * element;   // index of every group element
} smallgroup;

// calculate eigenvalues and characters
void bnd_eigcha(double * eig_hk, int * ndegen, dcomplex *** p2sym_chas, dcomplex *** p2sym_eigs, smallgroup * sgrp,
                double lattice[3][3], double rotations[][3][3], double translations[][3], int TR[], double rots_kd[][3][3], int nsymm,
                wannorb * orb_info, int flag_soc, 
                wanndata * hr, vector kpt, double degenerate_tolerance);

// calculate eigenvalues (eig_hk, dim: norb x 1) and eigenvectors (vr_hk, dim: norb x norb) of k-space Hamiltonian
void bnd_eig_hk(double * eig_hk, dcomplex * vr_hk, wanndata * hr, vector kpt);

void bnd_write_characters(const char * fnout, int ikpt, vector kpt, smallgroup * sgrp, 
                  double lattice[3][3], double rotations[][3][3], double translations[][3], int TR[],
                  double * eig_hk, int * ndegen, dcomplex ** sym_chas, int norb, int flag_dbgrp,
                  int en_print_prec, int bnd_print_len);

void bnd_write_eigenvalues(const char * fnout, int ikpt, vector kpt, smallgroup * sgrp, 
                  double lattice[3][3], double rotations[][3][3], double translations[][3], int TR[],
                  double * eig_hk, int * ndegen, dcomplex ** sym_eigs, int norb, int flag_dbgrp,
                  int en_print_prec, int bnd_print_len);

void bnd_write_bands(FILE * fbands, int norb, int nkpath, int nk_per_kpath, double lattice[3][3], vector * kvecs, double ** ebands);

#endif
