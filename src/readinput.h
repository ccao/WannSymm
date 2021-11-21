#ifdef READINPUT_H
#else
#define READINPUT_H


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>

#include "spglib.h"
#include "constants.h"
#include "wanndata.h"
#include "wannorb.h"
#include "vector.h"
#include "matrix.h"
#include "rotate_ham.h"
#include <mpi.h>

typedef struct __CrystalInfo {
    int natomtot;           // the total number of atoms
    int natomtype;          // the number of types of elements
    int natomeach[120];     // the number of atoms for each type of elements
    int * offsets;          // 1st orb offset of each atom in array of orbital info
    //vector lattice[3];      // 
    //int nsymm;
    //SymmetryOperator * symm;
} CrystalInfo;

typedef struct __SymmetryOperator {
    double rot[3][3];
    double krot[3][3];
    double trans[3];
    vector Rot[3];
    vector kRot[3][3];
    vector Trans;
    int TR;
} SymmetryOperator;


typedef struct __projgroup {
    // one projection group is one line in projection section of input file
    int npj;      // total projections written in this line
    int norb;     // number of orbital related to this line, e.g. (nonsoc) s::1 p:3  s;p::4  s;p;d::9
    char element[32];
    char orbname[64][32];
} projgroup;

void readinput(char * fn_input,
               int * p2flag_soc,
               char * seed,
               double lattice[3][3],
               double rotations[][3][3], 
               double translations[][3], 
               int TR[], 
               vector * magmom,
               int * p2number_of_symm, 
               wannorb ** p2orb_info, 
               int * p2num_wann, 
               vec_llist ** p2kpts,
               int * p2nkpt,
               vec_llist ** p2kpaths,
               char klabels[][SHORTLEN],
               int * p2nkpath,
               int * p2nk_per_kpath,
               int * p2flag_bands,
               int * p2flag_chaeig, 
               int * p2flag_chaeig_in_kpath,
               int * p2flag_restart, 
               int * p2flag_global_trsymm, 
               int * p2flag_expandrvec, 
               double * p2symm_magnetic_tolerance,
               double * p2ham_tolerance,
               double * p2degenerate_tolerance,
               int * p2flag_everysymm,
               int * p2flag_output_mem_usage,
               int * p2flag_symm_from_file, 
               char * fn_symm);

void parsestr_find_remove( char strout[MAXLEN], char strin[MAXLEN], char target);
// remove $target from sin, output is sout, sin and sout can be same

void parseline( char tag[MAXLEN], char arg[MAXLEN], char line[MAXLEN], int ignorecase);

void setup_codetype(int * p2code_type, char dftcode[MAXLEN]);

void read_pos_info(double lattice[3][3], 
                   double atom_positions[MAX_NUM_of_atoms][3],
                   char   name_of_atoms_each[128][4],
                   int    number_of_atoms_each[128],
                   int *  p2number_of_atoms_total,
                   int *  p2number_of_atomtypes,
                   char * fn_pos,
                   FILE * fin,
                   int begin_tag_exist);

void read_projection_info(projgroup pjgroup[MAXLEN],
                          int *  p2num_pjgroup,
                          FILE * fin);

void derive_projection_info(int *  p2num_wann, 
                            wannorb ** p2orb_info,
                            projgroup pjgroup[MAXLEN],
                            int    num_pjgroup,
                            double lattice[3][3],
                            double atom_positions[MAX_NUM_of_atoms][3],
                            char   name_of_atoms_each[128][4], 
                            int    number_of_atoms_each[128], 
                            int    number_of_atoms_total, 
                            int    number_of_atomtypes,
                            int    code_type,
                            int    flag_soc);

static void show_spg_dataset(double lattice[3][3],
                 const double origin_shift[3],
                 double position[][3],
                 const int num_atom,
                 const int types[],
                 FILE * fout);

int isletter(char candidate);
//vector find_pos(double pos[MAX_NUM_of_atoms][3], double lattice[3][3], int natoms, vector site);
void derive_magmom_from_string(double * magmom_array, char * magmom_string, int number_of_atoms_total, int flag_soc);
void derive_magmom_from_array(vector * magmom, double * magmom_array,int natom, int flag_soc, vector SAXIS);
int get_magnetic_type(vector magmom[], int natom , double eps);
void derive_symm_for_magnetic_materials(double rotations[][3][3], double translations[][3], int TR[], int * p2number_of_symm, double lattice[3][3], double atom_positions[][3], int * atom_types, int natom, vector * magmom, double symm_magnetic_tolerance );
int find_roted_atom(vector atom_roted, vector * atpos, int natom, double lattice[3][3]);
void remove_a_symm(double rotations[][3][3], double translations[][3], int TR[], int r);

//IO
int file_exists(const char * fn);
void print_error(char * msg);
void print_msg(char * msg);
void print_flag(FILE * fout, char * tag, int arg);
void print_symmetry(const char * fnout, double lattice[3][3], int isymm, double rot[3][3], double trans[3], int TR, int flag_showtrans, int flag_showmirror, int flag_dbgrp);

//string check and manipulate
int str2boolean(char * arg);

void read_kpath_info(vec_llist ** p2kpaths, char klabels[][SHORTLEN], int * p2nkpath, FILE * fin);
void read_kpts_info(vec_llist ** p2kpts, int * p2nkpt, FILE * fin);

// not implentmented yet
int readatag(char * tag, char * arg, char * fn, int ignorecase);

#endif
