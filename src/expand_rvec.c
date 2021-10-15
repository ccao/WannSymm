#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>

#include "constants.h"
#include "vector.h"
#include "wanndata.h"

#include "rotate_orbital.h"
#include "rotate_spinor.h"
#include "rotate_ham.h"

#include "mkl.h"
#include <mpi.h>


void rotate_ham(wanndata * he, wanndata * hin, double lattice[3][3], int nsymm, int rot_in[][3][3], double translation[][3], wannorb * orb_info)
{
    
}
