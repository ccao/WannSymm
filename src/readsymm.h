#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>

#include "spglib.h"
#include "constants.h"
#include "wanndata.h"
#include "vector.h"
#include "readinput.h"
#include "wannorb.h"

void readsymm(char * fn_symm, double rotations[][3][3], double translations[][3], int TR[], int * p2nsymm, int * p2flag_trsymm);
