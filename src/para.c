#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>

#include "constants.h"
#include "wanndata.h"
#include "wannorb.h"
#include "vector.h"
#include "para.h"

#include <mpi.h>

void init_para(){

    //define mpi_vector data type
    int nitems=3;
    int blocklengths[3]={1,1,1};
    MPI_Datatype types_vector[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    //MPI_Datatype mpi_vector;
    MPI_Aint     offsets[3];
    offsets[0] = offsetof(vector, x);
    offsets[1] = offsetof(vector, y);
    offsets[2] = offsetof(vector, z);
    MPI_Type_create_struct(nitems, blocklengths, offsets, types_vector, &mpi_vector);
    MPI_Type_commit(&mpi_vector);

    //define mpi_wannorb data type
    nitems=6;
    int blocklengths_wannorb[6]={1,3,1,1,1,1};
    MPI_Datatype types_wannorb[6] = {mpi_vector, mpi_vector, MPI_INT, MPI_INT, MPI_INT, MPI_INT};
    //MPI_Datatype mpi_wannorb;
    MPI_Aint     offsets_wannorb[6];
    offsets_wannorb[0] = offsetof(wannorb, site);
    offsets_wannorb[1] = offsetof(wannorb, axis);
    offsets_wannorb[2] = offsetof(wannorb, r);
    offsets_wannorb[3] = offsetof(wannorb, l);
    offsets_wannorb[4] = offsetof(wannorb, mr);
    offsets_wannorb[5] = offsetof(wannorb, ms);
    MPI_Type_create_struct(nitems, blocklengths_wannorb, offsets_wannorb, types_wannorb, &mpi_wannorb);
    MPI_Type_commit(&mpi_wannorb);

}
