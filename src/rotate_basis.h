#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>

#include "constants.h"
#include "vector.h"
#include "wanndata.h"

#include "rotate_orbital.h"
#include "rotate_spinor.h"
#include "rotate_ham.h"        //use getrvec_and_site()

#include <mkl.h>
#include <mpi.h>

void get_sym_op_reciprocalspace(dcomplex * sym_op, double lattice[3][3], wannorb * orb_info, int norb, int isym_in_doublegp, double rotation[3][3],double translation[3], int TR, double rot_kd[3][3], vector kpt, int flag_soc, int flag_local_axis);