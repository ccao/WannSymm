#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#define MAXLEN 256
#define MAXSPC 20

#define TWOPI  6.283185307179586

using namespace std;


typedef struct __wannorb{
  int atmidx;
  char orb_l;
  int norb;
} wannorb;

void matrix_multiply(double * out, double lhs[3][3], double rhs[3][3]) {
  for(int ii=0; ii<3; ii++) {
    for(int jj=0; jj<3; jj++) {
      out[ii*3+jj]=lhs[ii][0]*rhs[0][jj]+lhs[ii][1]*rhs[1][jj]+lhs[ii][2]*rhs[2][jj];
    }
  }
}

int read_inwf_file(FILE * fp, wannorb ** orblist) {
  char line[MAXLEN];
  int norb;

  fgets(line, MAXLEN, fp);
  fgets(line, MAXLEN, fp);
  fgets(line, MAXLEN, fp);
  sscanf(line, " %*d %d", &norb);

  int current_idx=-1, idx;
  char current_l=0, ll;

  int ngrp=0, nskip;

  for(int ii=0; ii<norb; ii++) {
    char * lp;
    fgets(line, MAXLEN, fp);
    sscanf(line, "%d ", &nskip);

    if (nskip!=0) {

      lp=strstr(line, "->");
      sscanf(lp, "-> %d:%c", &idx, &ll);

      if (current_idx!=idx || ll!=current_l) {
        ngrp++;
        current_idx=idx;
        current_l=ll;
      }

      for(int jj=0; jj<nskip; jj++) {
        fgets(line, MAXLEN, fp);
      }
    }
  }

  (*orblist)=(wannorb *)malloc(sizeof(wannorb)*ngrp);
  rewind(fp);

  ngrp=0; current_idx=-1; current_l=0;

  fgets(line, MAXLEN, fp);
  fgets(line, MAXLEN, fp);
  fgets(line, MAXLEN, fp);

  for(int ii=0; ii<norb; ii++) {
    char * lp;
    fgets(line, MAXLEN, fp);
    sscanf(line, "%d ", &nskip);

    if (nskip!=0) {

      lp=strstr(line, "->");
      sscanf(lp, "-> %d:%c", &idx, &ll);

      if (current_idx!=idx || ll!=current_l) {
        current_idx=idx;
        current_l=ll;
        ((*orblist)+ngrp)->atmidx=idx-1;
        ((*orblist)+ngrp)->orb_l=ll;
        ((*orblist)+ngrp)->norb=1;
        ngrp++;
      }
      else {
        ((*orblist)+ngrp-1)->norb++;
      }

      for(int jj=0; jj<nskip; jj++) {
        fgets(line, MAXLEN, fp);
      }
    }
  }

  return ngrp;
}


int read_wout_file(FILE * fp, double ** xpos) {
  char line[MAXLEN];
  int npos=-1, state=0;

  double bvec[3][3];
  double xx[3];

  while(!feof(fp) && state<2) {
    fgets(line, MAXLEN, fp);
    switch(state) {
      case 0:
        if (strstr(line, "Reciprocal-Space Vectors")) {
          char * lp;
          fgets(line, MAXLEN, fp);
          lp=strstr(line, "b_1");
          sscanf(lp, "b_1 %lf %lf %lf", bvec[0], bvec[0]+1, bvec[0]+2);
          fgets(line, MAXLEN, fp);
          lp=strstr(line, "b_2");
          sscanf(lp, "b_2 %lf %lf %lf", bvec[1], bvec[1]+1, bvec[1]+2);
          fgets(line, MAXLEN, fp);
          lp=strstr(line, "b_3");
          sscanf(lp, "b_3 %lf %lf %lf", bvec[2], bvec[2]+1, bvec[2]+2);

          for(int ii=0; ii<3; ii++) {
            for(int jj=0; jj<3; jj++) bvec[ii][jj]=bvec[ii][jj]/TWOPI;
          }
        }
        else if (strstr(line, "Final State")) {
          state=1;
          npos++;
        }
        break;
      case 1:
        if (strstr(line, "WF centre and spread")) {
          npos++;
        }
        else if (strstr(line, "Sum of centres and spreads")) {
          state=2;
        }
        break;
      default:
        printf("!!! SHOULD NOT HAPPEN!\n");
        exit(0);
    }
  }

  (*xpos)=(double *) malloc(sizeof(double)*3*npos);

  rewind(fp);
  state=0; npos=-1;

  while(!feof(fp) && state<2) {
    fgets(line, MAXLEN, fp);
    switch(state) {
      case 0:
        if (strstr(line, "Final State")) {
          state=1;
          npos++;
        }
        break;
      case 1:
        if (strstr(line, "WF centre and spread")) {
          sscanf(line, " WF centre and spread %*d ( %lf, %lf, %lf)", xx, xx+1, xx+2);
          for (int ii=0; ii<3; ii++) {
            (*xpos)[3*npos+ii]=xx[0]*bvec[ii][0]+xx[1]*bvec[ii][1]+xx[2]*bvec[ii][2];
          }
          npos++;
        }
        else if (strstr(line, "Sum of centres and spreads")) {
          state=2;
        }
        break;
      default:
        printf("!!! SHOULD NOT HAPPEN!\n");
        exit(0);
    }
  }
  return npos;
}

int read_output0(FILE * fin, double ** locaxis) {
  int nat=0;

  double locrot[3][3];
  double atmrot[3][3];

  char line[MAXLEN];

  int status=0;

  while(!feof(fin) && status<3) {
    fgets(line, MAXLEN, fin);
    switch(status) {
      case 0:
        if (strstr(line, "NOT EQUIV ATOM")) {
          status=1;
          for(int ii=0; ii<3; ii++) {	// READ LOCAL ROTATION
            fgets(line, MAXLEN, fin);
          }
        }
        break;
      case 1:
        if (strstr(line, "EQUIV ATOM")) {
          for(int ii=0; ii<3; ii++) {	// READ ROTATION MATRIX
            fgets(line, MAXLEN, fin);
          }
          nat++;
        }
        else if (strcmp(line, "")) {
          status=0;
        }
        else {
          status=3;
        }
        break;
    }
  }

  (*locaxis)=(double *) malloc(sizeof(double)*nat*9);
  status=0;
  nat=0;
  rewind(fin);

  while(!feof(fin) && status<3) {
    fgets(line, MAXLEN, fin);
    switch(status) {
      case 0:
        if (strstr(line, "NOT EQUIV ATOM")) {
          status=1;
          for(int ii=0; ii<3; ii++) {   // READ LOCAL ROTATION
            fgets(line, MAXLEN, fin);
            sscanf(line, " %lf %lf %lf", locrot[0]+ii, locrot[1]+ii, locrot[2]+ii);
          }
        }
        break;
      case 1:
        if (strstr(line, "EQUIV ATOM")) {
          for(int ii=0; ii<3; ii++) {   // READ ROTATION MATRIX
            fgets(line, MAXLEN, fin);
            sscanf(line, " %lf %lf %lf", atmrot[0]+ii, atmrot[1]+ii, atmrot[2]+ii);
          }
          matrix_multiply((*locaxis)+9*nat, locrot, atmrot);
          nat++;
        }
        else if (strcmp(line, "")) {
          status=0;
        }
        else {
          status=3;
        }
        break;
    }
  }

  return nat;
}

int main(int argc, char ** argv) {
  int  norb, nat;
  double * xx;
  double * locaxis;

  int nn;
  wannorb * orblist;

  FILE * fin;

  if (argc<2) {
    printf(" Usage: mkw2kinput.x WOUT_FILE OUTPUT0 [INWF_FILE] \n");
    printf("   WOUT_FILE is the output of wannier90 containing WCC information.\n");
    printf("   OUTPUT0 is the output of lapw0.\n");
    printf("   Optional INWF_FILE is the inwf file when perform W2W.\n");
    exit(0);
  }

  fin=fopen(argv[1], "r");
  norb=read_wout_file(fin, &xx);
  fclose(fin);

  fin=fopen(argv[2], "r");
  nat=read_output0(fin, &locaxis);
  fclose(fin);

  if (argc>=2) {
    fin=fopen(argv[3], "r");
    nn=read_inwf_file(fin, &orblist);
    fclose(fin);
  }
  else {
    printf(" No INWF file specified, extra input is required.\n");
    printf(" Please input Wannier orbital local axis specification:\n");
    printf("   Input total number of specification lines:\n");
    scanf(" %d", &nn);
    orblist=(wannorb *)malloc(sizeof(wannorb)*nn);

    printf("   Now input the specifications in the following form:\n");
    printf("     ATOM_INDEX  ORBITAL_TYPE\n");
    printf("     For example:\n");
    printf("     1  d\n");
    printf("     2  f\n");
    printf("     ... \n");
    printf("   Please be noted that the atom list should be in the same order as inwf file.\n");
    printf("     which could be in different order as struct file.\n");

    for(int ii=0; ii<nn; ii++) {
      scanf(" %d %c", &((orblist+ii)->atmidx), &((orblist+ii)->orb_l));
      (orblist+ii)->atmidx--;

      switch((orblist+ii)->orb_l) {
        case 's':
          (orblist+ii)->norb=1;
          break;
        case 'p':
          (orblist+ii)->norb=3;
          break;
        case 'd':
          (orblist+ii)->norb=5;
          break;
        case 'f':
          (orblist+ii)->norb=7;
          break;
      }
    }
  }

  // Headers....
  printf("######### WannSymm INPUT File Follows ########\n");
  printf("# Modify according to your system\n");
  printf("# template input file of wannsymm.x\n");
  printf("# anything following '#', '!' or '//' in a line will be regard as comments\n");
  printf("# tag names are case insensitive( SeedName and seednAme are equivalent)\n");
  printf("\n");
  printf("DFTcode  = VASP\n");
  printf("\n");
  printf("Spinors  = T\n");
  printf("\n");
  printf("SeedName ='wannier90'\n");
  printf("\n");
  printf("Use_POSCAR = 'POSCAR'\n");
  printf("\n");
  printf("# Projections_in_Format_of_wannier90\n");
  printf("\n");

  // Projection ...
  printf("begin projections\n");
  int head_orb=0;
  for(int ii=0; ii<nn; ii++) {
    int idx=(orblist+ii)->atmidx;
    char ll=(orblist+ii)->orb_l;

    printf("f= %9.5f,%9.5f,%9.5f : %c : z= %9.5f,%9.5f,%9.5f : x= %9.5f,%9.5f,%9.5f : y= %9.5f,%9.5f,%9.5f\n", 
              xx[3*head_orb], xx[3*head_orb+1], xx[3*head_orb+2], ll,
              locaxis[idx*9+6], locaxis[idx*9+7], locaxis[idx*9+8],
              locaxis[idx*9+0], locaxis[idx*9+1], locaxis[idx*9+2],
              locaxis[idx*9+3], locaxis[idx*9+4], locaxis[idx*9+5]);
    head_orb+=(orblist+ii)->norb;
  }
  printf("end projections\n");

  // End ...
  printf("\n");
  printf("restart = F\n");
  printf("# Kpoint for calculating band's character of every symmetry\n");
  printf("#kpt = 0 0 0\n");
  printf("#chaeig_in_kpath = T\n");
  printf("\n");
  printf("nk_per_kpath = 101\n");
  printf("\n");
  printf("beginkpath\n");
  printf("endkpath\n");


  free(orblist);
  free(xx);
  free(locaxis);

  return 0;
}
