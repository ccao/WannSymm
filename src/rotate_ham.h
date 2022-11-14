#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>

#include "constants.h"
#include "vector.h"
#include "matrix.h"
#include "wanndata.h"
#include "readinput.h"
#include "rotate_orbital.h"
#include "rotate_spinor.h"

#include "mkl.h"
#include <mpi.h>


void rotate_ham(wanndata * hout, wanndata * hin, double lattice[3][3], double rotation[3][3], double translation[3], wannorb * orb_info, int flag_soc, int flag_local_axis, int index_of_sym);
//--ABANDON--//vector convert_site_to_ws(vector v);
//--ABANDON--//vector convert_site_back_from_ws(vector v);
void getrvec_and_site(vector * p2rvec, vector * p2site, vector loc, wannorb * orb_info, int norb, double lattice[3][3]);
void get_axis_angle_of_rotation(double axis[3], double * angle, int * inv, double rin[3][3], double lattice[3][3]);
int type_of_rotation(double rot[3][3]);
void inverse_symm(double rin[3][3], double rout[3][3], double tin[3], double tout[3]);
double sign(double in);
void trsymm_ham(wanndata * hout, wanndata * hin, wannorb * orb_info, int flag_soc);
void combine_rot_with_local_axis(double rot_combined[3][3], double rotation[3][3], double lattice[3][3], wannorb * orb_info, int io_in, int io_out);
