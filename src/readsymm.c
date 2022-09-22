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

#include <mpi.h>

//#define __DEBUG

void readsymm(char * fn_symm, double rotations[][3][3], double translations[][3], int TR[], int * p2nsymm, int * p2flag_global_trsymm)
{
    FILE * fin;
    char line[MAXLEN];
    char tag[MAXLEN];
    char arg[MAXLEN];
    char msg[MAXLEN];
    int i,j,k;
    int flag_tmp;
    double rotangle;

    fin=fopen(fn_symm, "r");

    while( fgets(line, MAXLEN, fin) ){
        parseline(tag, arg, line, 1);
        if(strcmp(tag, "globaltime-reversalsymmetry") == 0 || 
           strcmp(tag, "trsymm") == 0 || strcmp(tag, "global_trsymm") == 0){
            if( (flag_tmp = str2boolean(arg)) != -1 )
                *p2flag_global_trsymm = flag_tmp;
            else{
                sprintf(msg, "ERROR: %s = %s --- undefined value, should be replaced by T or F\n", tag, arg);
                print_error(msg);
                exit(1);
            }
        }
        else if(strcmp(tag, "nsymm") == 0){
            sscanf(arg, "%d", p2nsymm);
            break;
        }
    }
    

    for(i=0;i<*p2nsymm;i++){
        fgets(line, MAXLEN, fin);
        parseline(tag, arg, line, 1);
        if( strcmp(tag, "angle") == 0 ){
            sscanf(arg, "%lf", &rotangle);
            // convert angle into rotation matrix
            // this part has not been finished 
        } else {
            for(j=0;j<3;j++){
                fgets(line, MAXLEN, fin);
                sscanf(line, "%lf%lf%lf", &rotations[i][j][0], &rotations[i][j][1], &rotations[i][j][2]);
            }
        }
        fgets(line, MAXLEN, fin);
        sscanf(line, "%lf%lf%lf", &translations[i][0], &translations[i][1], &translations[i][2]);
        if ( strchr(line, 'T') ==NULL && strchr(line, 't') == NULL ){
            TR[i] = 0;
        } else if(strchr(line, 'F') !=NULL && strchr(line, 'f') != NULL)
        {
            TR[i] = 0;
        } else {
            TR[i] = 1;
            *p2flag_global_trsymm = 0; //if some symmetry contains time reversal, then the global time-reversal symmetry won't exist.
        }
    }
    fclose(fin);
}
