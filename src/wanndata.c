#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>

#include "constants.h"
#include "vector.h"
#include "wanndata.h"
#include "wannorb.h"
#include "readinput.h"

void init_wanndata(wanndata * wann) {
  wann->ham=(double complex *) malloc(sizeof(double complex)*wann->norb*wann->norb*wann->nrpt);
  wann->hamflag = (int *) malloc(sizeof(int)*wann->norb*wann->norb*wann->nrpt);
  wann->rvec=(vector *) malloc(sizeof(vector)*wann->nrpt);
  wann->weight=(int *) malloc(sizeof(int)*wann->nrpt);
   
  int i;
  //set default vaule
  for(i=0;i<wann->nrpt*wann->norb*wann->norb;i++){
      wann->ham[i]     = 0;
      wann->hamflag[i] = 0;
  }
  for(i=0;i<wann->nrpt;i++){
      wann->weight[i] = 1;
  }
}

void read_ham(wanndata * wann, char * seed) {
  FILE * fin;
  int irpt, iorb_right, iorb_left;
  int ielement;
  char line[MAXLEN];
  char * pch;
  int t1, t2, t3;
  vector trvec;
  double a, b;
  char msg[MAXLEN];

  strcpy(line, seed);
  strcat(line, "_hr.dat");

  if( ! file_exists(line) ){
      sprintf(msg, "ERROR!!! trying to read the file \"%s\", but it can not be found\n", line);
      print_error(msg);
      exit(1);
  }
  fin=fopen(line, "r");

  fgets(line, MAXLEN, fin);
  fgets(line, MAXLEN, fin);
  sscanf(line, " %d", &(wann->norb));
  fgets(line, MAXLEN, fin);
  sscanf(line, " %d", &(wann->nrpt));

  init_wanndata(wann);

  for(irpt=0; irpt<wann->nrpt; irpt++) {
    if(irpt%15==0) {
      fgets(line, MAXLEN, fin);
      pch=strtok(line, " ");
    }
    else {
      pch=strtok(NULL, " ");
    }
    sscanf(pch, "%d", wann->weight+irpt);
  }

  for(irpt=0; irpt<wann->nrpt; irpt++) {
      for( ielement = 0; ielement < wann->norb*wann->norb; ielement++){
          fgets(line, MAXLEN, fin);
          sscanf(line, " %d %d %d %d %d %lf %lf ", &t1, &t2, &t3, &iorb_left, &iorb_right, &a, &b);
          iorb_left--;
          iorb_right--;
          if (iorb_right==0 && iorb_left==0)
              init_vector(wann->rvec+irpt, t1, t2, t3);

          wann->ham[irpt*wann->norb*wann->norb+iorb_right*wann->norb+iorb_left]=a+_Complex_I*b;
          wann->hamflag[irpt*wann->norb*wann->norb+iorb_right*wann->norb+iorb_left] = 1;
      }
  }

  fclose(fin);
}

void finalize_wanndata(wanndata wann) {
  free(wann.rvec);
  free(wann.weight);
  free(wann.ham);
  free(wann.hamflag);
}

void write_ham(wanndata * wann, char * seed) {
  int irpt, iorb, jorb;
  FILE * fout;
  char fn[MAXLEN];

  //if (fn==NULL) {
  //  fout=stdout;
  //}
  sprintf(fn, "%s_hr.dat", seed);
  fout=fopen(fn, "w");
  //fprintf(fout, "# Wannier Hamiltonian extended to supercell\n");
  fprintf(fout, "# symmetrized Hamiltonian \n");
  fprintf(fout, "%5d\n%5d", wann->norb, wann->nrpt);
  for(irpt=0; irpt<wann->nrpt; irpt++) {
    if(irpt%15==0) fprintf(fout, "\n");
    fprintf(fout, "%5d", wann->weight[irpt]);
  }
  fprintf(fout, "\n");

  for(irpt=0; irpt<wann->nrpt; irpt++) {
    for(iorb=0; iorb<wann->norb; iorb++) {
      for(jorb=0; jorb<wann->norb; jorb++) {
        fprintf(fout, "%5d%5d%5d%5d%5d%22.16f%22.16f\n", 
          (int)(wann->rvec+irpt)->x, (int)(wann->rvec+irpt)->y, (int)(wann->rvec+irpt)->z, 
          jorb+1, iorb+1,
          creal(wann->ham[irpt*wann->norb*wann->norb+iorb*wann->norb+jorb]),
          cimag(wann->ham[irpt*wann->norb*wann->norb+iorb*wann->norb+jorb]));
      }
    }
  }
  //if (fn)  
  fclose(fout);
}

void write_reduced_ham(wanndata * wann, char * seed) {
  int irpt, iorb, jorb;
  int iirpt;
  long long int nwann;
  vector rv;
  char fn[MAXLEN];

  FILE * fout;

  //if (fn==NULL) {
  //  fout=stdout;
  //}
  //else {
  //}
  sprintf(fn, "%s_hr.dat", seed);
  fout=fopen(fn, "w");

  fprintf(fout, "# Reduced Wannier Hamiltonian (only nonzero value)\n");
  fprintf(fout, "%5d\n%5d", wann->norb, wann->nrpt);
  for(irpt=0; irpt<wann->nrpt; irpt++) {
    if(irpt%15==0) fprintf(fout, "\n");
    fprintf(fout, "%5d", wann->weight[irpt]);
  }
  fprintf(fout, "\n");

  nwann=0;

  for(irpt=0; irpt<wann->nrpt; irpt++) {
    fprintf(fout, "%5d%5d%5d\n", (int)(wann->rvec+irpt)->x, (int)(wann->rvec+irpt)->y, (int)(wann->rvec+irpt)->z);
    for(iorb=0; iorb<wann->norb; iorb++)
      for(jorb=0; jorb<=iorb; jorb++)
        if(cabs(wann->ham[irpt*wann->norb*wann->norb+iorb*wann->norb+jorb])>eps8)
          nwann++;
  }

  fprintf(fout, "%lld\n", nwann);

  for(irpt=0; irpt<wann->nrpt; irpt++) {
    rv=wann->rvec[irpt];
    rv.x=-rv.x; rv.y=-rv.y; rv.z=-rv.z;
    iirpt=find_vector(rv, wann->rvec, wann->nrpt);
    for(iorb=0; iorb<wann->norb; iorb++) {
      for(jorb=0; jorb<=iorb; jorb++) {
        if(cabs(wann->ham[irpt*wann->norb*wann->norb+iorb*wann->norb+jorb])>eps8) {
          fprintf(fout, "%5d%5d%5d%5d%22.16f%22.16f\n",
            iorb+1, jorb+1, irpt+1, iirpt+1,
            creal(wann->ham[irpt*wann->norb*wann->norb+iorb*wann->norb+jorb]),
            cimag(wann->ham[irpt*wann->norb*wann->norb+iorb*wann->norb+jorb]));
          if(fabs(creal(wann->ham[irpt*wann->norb*wann->norb+iorb*wann->norb+jorb])-
                  creal(wann->ham[iirpt*wann->norb*wann->norb+jorb*wann->norb+iorb]))>eps8 ||
             fabs(cimag(wann->ham[irpt*wann->norb*wann->norb+iorb*wann->norb+jorb])+
                  cimag(wann->ham[iirpt*wann->norb*wann->norb+jorb*wann->norb+iorb]))>eps8) 
            fprintf(fout, "!!!!WARNING: Unhermitian Hamiltonian:%5d%5d%5d%5d%5d\n", iorb+1, jorb+1, (int)(wann->rvec+irpt)->x, (int)(wann->rvec+irpt)->y, (int)(wann->rvec+irpt)->z);
        }
      }
    }
  }
  //if (fn) 
  fclose(fout);
}

int find_index_of_ham(wanndata * wann, wannorb * orb_info, int num_wann, int irpt, vector site1, int r1, int l1, int mr1, int ms1, vector site2, int r2, int l2, int mr2, int ms2){
    //
    int iorb, jorb;
    int norb=wann->norb;
    //irpt = find_vector(rvec, wann->rvec, wann->nrpt);
    if (irpt<0) return -1;
    iorb = find_index_of_wannorb(orb_info, num_wann, site1, r1, l1, mr1, ms1);
    if (iorb<0) return -2;
    jorb = find_index_of_wannorb(orb_info, num_wann, site2, r2, l2, mr2, ms2);
    if( jorb<0) return -3;
    return irpt*norb*norb+jorb*norb+iorb;
}


void wanndata_add(wanndata * hout, wanndata * hin1, wanndata * hin2){
    //find every  hin2.rvec in the list of  hin1.rvec
    //if found, 
    //if not found, hout.nrpt = hin1.nrpt + num_not_found
    wanndata lhin1, lhin2;  //local virables
    
}


void hamblock_init(hamblock * hblk, int norb){
  if(norb>0) hblk->norb=norb;
  if(hblk->norb<=0){
    print_error("ERROR: In initializing hamblock, norb is not positive.");
    exit(1);
  }
  hblk->ham = (double complex *) malloc(sizeof(double complex) * hblk->norb * hblk->norb);
}

void hamblock_finalize(hamblock * hblk){
  free(hblk->ham);
}


void ham_R_init(ham_R * hr, int norb, int nrpt){
  // allocate memory to ham_R
  int ii;
  if(norb>0) hr->norb=norb;
  if(nrpt>0) hr->nrpt=nrpt;
  if(hr->norb<=0 || hr->nrpt<=0){
    // error
    print_error("ERROR: In initializing hamiltonian(ham_R), norb or nrpt is not positive.");
    exit(1);
  }
  hr->hblk=(hamblock *) malloc(sizeof(hamblock)*hr->nrpt);
  for(ii=0; ii<hr->nrpt; ii++){
    hamblock_init(&(hr->hblk[ii]), hr->norb);
  }
  hr->rvec=(vector *) malloc(sizeof(vector)*hr->nrpt);
  hr->weight=(int *) malloc(sizeof(int) * hr->nrpt);
  for(ii=0; ii<hr->nrpt; ii++){
    hr->weight[ii] = 1;
  }
}

void ham_R_finalize(ham_R * hr){
  int ii;
  for(ii=0; ii<hr->nrpt; ii++){
    hamblock_finalize(&(hr->hblk[ii]));
  }
  free(hr->hblk);
  free(hr->rvec);
  free(hr->weight);
}

void wanndata2ham_R(ham_R * hr, wanndata * wann){
  int irpt;
  hr->norb = wann->norb;
  hr->nrpt = wann->nrpt;
  ham_R_init(hr, -1, -1);
  memcpy(hr->rvec, wann->rvec, sizeof(vector) * hr->nrpt);
  memcpy(hr->weight, wann->weight, sizeof(int) * hr->nrpt);
  for(irpt=0; irpt<hr->nrpt; irpt++){
    memcpy(hr->hblk[irpt].ham, wann->ham+hr->norb*hr->norb*irpt, sizeof(double complex) * hr->norb * hr->norb);
  }
}
void ham_R2wanndata(wanndata * wann, ham_R * hr){
  int irpt;
  wann->norb = hr->norb;
  wann->nrpt = hr->nrpt;
  init_wanndata(wann);
  memcpy(wann->rvec, hr->rvec, sizeof(vector) * hr->nrpt);
  memcpy(wann->weight, hr->weight, sizeof(int) * hr->nrpt);
  for(irpt=0; irpt<hr->nrpt; irpt++){
    memcpy(wann->ham+hr->norb*hr->norb*irpt, hr->hblk[irpt].ham, sizeof(double complex) * hr->norb * hr->norb);
  }
}

void ham_R_addrvec(ham_R * hr, vector * rvec, int nrvec_add){
    hr->nrpt = hr->nrpt + nrvec_add;
}

void ham_R_delrvec(ham_R * hr, vector rvec);
