#ifdef WANNDATA_H
#else
#define WANNDATA_H

#include <complex.h>

#include "vector.h"
#include "wannorb.h"

typedef struct __wanndata {
  int norb;
  int nrpt;
  double complex * ham;
  int * hamflag;    //indicate if wann.ham[i] is enabled(meaningful) or not by checking wann.hamflag[i] equals '1' or '0'
                    //hamflag is used for counting the number of rotated hamiltonian containing a certain R vector
  vector * rvec;
  int * weight;

  vector kvec;      //only for hamk
} wanndata;

void init_wanndata(wanndata * wann);
void read_ham(wanndata * wann, char * seed);
void finalize_wanndata(wanndata wann);
void write_ham(wanndata * wann, char * seed);
void write_reduced_ham(wanndata * wann, char * seed);

int find_index_of_ham(wanndata * wann, wannorb * orb_info, int num_wann, int irpt, vector site1, int r1, int l1, int mr1, int ms1, vector site2, int r2, int l2, int mr2, int ms2);

// below are not fully implentmented yet
typedef struct __hamblock {
    int norb;
    double complex * ham;
} hamblock;

void hamblock_init(hamblock * hblk, int norb);
void hamblock_finalize(hamblock * hblk);

typedef struct __ham_R {
    // same as wanndata in meaning, but the hamiltonian is stored in block of R vector
    int norb;
    int nrpt;
    hamblock * hblk;
    vector  * rvec;
    int * weight;
} ham_R;

void ham_R_init(ham_R * hr, int norb, int nrpt);
void ham_R_finalize(ham_R * hr);
void wanndata2ham_R(ham_R * hr, wanndata * wann);
void ham_R2wanndata(wanndata * wann, ham_R * hr);
void ham_R_read(ham_R * hr, char * seed);
void ham_R_write(ham_R * hr, char * seed);
void ham_R_addrvec(ham_R * hr, vector * rvec, int nrvec_add);
void ham_R_delrvec(ham_R * hr, vector rvec);

int combine_wanndata(wanndata * out, wanndata * wup, wanndata * wdn );

#endif
