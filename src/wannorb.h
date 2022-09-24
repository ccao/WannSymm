#ifndef __WANNORB_H
#define __WANNORB_H

#include "constants.h"
#include "vector.h"

typedef struct __wannorb {
  vector site;	/* Orbital site */
  vector axis[3];
  int r;    /*main quantum number? not used in symmetry*/
  int l;	/* Maximal angular momentum for this orbital */
  int mr;   /*indicator of cubic harmonic orbital function*/
  int ms;    /*quantum number of spin, 0 for spin up, 1 for down*/
//  vector * dir;	/* Orbital definition*/
/*
 * l & dir definition:
 * s : l=0; dir=NULL (s orbital do not have direction
 * p : l=1; px: dir=(1,0,0) etc...
 * sp2 : l=1; sp2-1: (1, 0, 0) etc...
 * d : l=2; dzx: dir=(1,0,0) and (0,0,1) ...
 */
} wannorb;

void init_wannorb(wannorb * orb,vector site, int l, int mr, int ms, int r, vector axisz, vector axisx, vector axisy);
int find_index_of_wannorb(wannorb * wann, int num_wann, vector site, int r, int l, int mr, int ms);//ignored axis in version 0.01
//void init_wannorb(wannorb * wann, vector * v, int l);
//void copy_wannorb(wannorb * target, wannorb source);
//void finalize_wannorb(wannorb wann);
//int reappear(int * list, int length);
//void symmop_wannorb(wannorb * out, wannorb in, vector shift, vector * symm);
//int match_wannorb(vector * vt, wannorb orb1, wannorb orb2);

#endif   /* __WANNORB_H*/
