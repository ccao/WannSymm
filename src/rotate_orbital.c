#include "rotate_orbital.h"

//#define __DEBUG_rot_orb

void generate_Lmatrix(dcomplex * Lmat, int l) {
/*
 * This function generates the L operator in |lm> basis
 *   which combined with definition of real-space harmonics
 *   can be used for spacial rotation generation.
 *   ASSUMES hbar=1 !!!
 * OUTPUT:
 *   Lmat : L operator 3 x (2l+1) x (2l+1) matrix
 *          Should be allocated before calling.
 * INPUT:
 *   l    : angular momentum
 *          l=0 : s
 *          l=1 : p
 *          l=2 : d
 *          l=3 : f
 */

  const int N=2*l+1;
  dcomplex * lp;	// L+ = L_x + i L_y
  dcomplex * lm;	// L- = L_x - i L_y
  dcomplex * lz;	// lz is L_z

  lp=(dcomplex *) malloc (sizeof(dcomplex)*N*N);
  lm=(dcomplex *) malloc (sizeof(dcomplex)*N*N);
  lz=(dcomplex *) malloc (sizeof(dcomplex)*N*N);

  memset(lp, 0, sizeof(dcomplex)*N*N);
  memset(lm, 0, sizeof(dcomplex)*N*N);
  memset(lz, 0, sizeof(dcomplex)*N*N);
  int m;

  for (m=-l; m<l+1; m++) {
    // Lz|lm> = m |lm>
    lz[(m+l)*(N+1)]=m;

    // L+|lm> = sqrt[(l-m)(l+m+1)] |lm+1>
    if ( m<l ) {
      //lp[(m+l)*N+(m+l+1)]=sqrt((l-m)*(l+m+1));    //wrong
      lp[(m+l+1)*N+(m+l)]=sqrt((l-m)*(l+m+1));      //right
    }

    // L-|lm> = sqrt[(l+m)(l-m+1)] |lm-1>
    if ( m>-l ) {
      //lm[(m+l)*N+(m+l-1)]=sqrt((l+m)*(l-m+1));    //wrong
      lm[(m+l-1)*N+(m+l)]=sqrt((l+m)*(l-m+1));      //right
    }
  }

  for (m=0; m<N*N; m++) {
    // L_x = (L+ + L-)/2
    Lmat[m]=(lp[m]+lm[m])/2.0;
    // L_y = (L+ - L-)/2i
    Lmat[N*N+m]=(lp[m]-lm[m])/(2.0*I);
    // L_z = lz
    Lmat[2*N*N+m]=lz[m];
  }


  free(lz);
  free(lm);
  free(lp);
}

void rotate_Ylm(dcomplex * rot, int l, double axis[3], double alpha, int inv) {
/*
 * This function generates the rotation matrix for Ylm with a specific l
 */

  const int N=2*l+1;
  const dcomplex one=1.0;
  const dcomplex zero=0.0;

  dcomplex * Lmat;
  dcomplex * Ldotn;

  dcomplex * eig;
  dcomplex * vl;
  dcomplex * vr;

  int info;

  Lmat=(dcomplex *) malloc (sizeof(dcomplex)*3*N*N);
  Ldotn=(dcomplex *) malloc (sizeof(dcomplex)*N*N);
  eig=(dcomplex *) malloc (sizeof(dcomplex)*N);
  vl=(dcomplex *) malloc (sizeof(dcomplex)*N*N);
  vr=(dcomplex *) malloc (sizeof(dcomplex)*N*N);


  // Generate L matrix
  generate_Lmatrix(Lmat, l);

  // calculates -i*alpha* n\cdot L
  int i;
  for (i=0; i<N*N; i++) {
    Ldotn[i]= -I * alpha * (axis[0]*Lmat[i] + axis[1]*Lmat[N*N+i] + axis[2]*Lmat[2*N*N+i]);
  }



  // Diagonalize and get eigen vectors
  //error.old//info = LAPACKE_zgeev( LAPACK_ROW_MAJOR, 'N', 'V', N, Ldotn, N, eig, vl, N, vr, N);
  info = LAPACKE_zgeev( LAPACK_ROW_MAJOR, 'V', 'V', N, Ldotn, N, eig, vl, N, vr, N);


  // exp ( M ) = X^-1 exp(A) X
  memset(rot, 0, sizeof(dcomplex)*N*N);
  for (i=0; i<N; i++) {
    rot[i*N+i]=cexp(eig[i]);
  }

  dcomplex * tmp;
  tmp = (dcomplex *) malloc (sizeof(dcomplex)*N*N);

  //rot = vr * rot * ConjTrans(vr)
  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, N, N, N, &one, rot, N, vr, N, &zero, Ldotn, N);
  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, &one, vr, N, Ldotn, N, &zero, rot, N);
  //cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, &one, rot, N, vr, N, &zero, tmp, N);
  //cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans, N, N, N, &one, vl, N, tmp, N, &zero, rot, N);
  free(tmp);
 
  // Not sure if it is needed...
  if (inv) {
    if(l%2 == 1)
        for (i=0; i<N*N; i++)
          rot[i] *= -1;
    //rot[i]=-rot[i];
    //rot[i]=(l%2==1?-rot[i]:rot[i]);
  }

  free(vr);
  free(vl);
  free(eig);
  free(Ldotn);
  free(Lmat);
}

void generate_C2Ylm(dcomplex * C2Ylm, dcomplex * Ylm2C, int l) {
/*
 * This function generates the transform matrix between Ylm and cubic harmonics
 *   OUTPUT:
 *     C2Ylm: From Cubic Harmonics to Ylm, generated if is not NULL
 *     Ylm2C: From Ylm to cubic Harmonics, generated if is not NULL
 *   INPUT:
 *     l: angular momentum
 *
 *   Definition (and order) of cubic harmonics
 *   l=0
 *     s  = |0 0>
 *   l=1
 *     pz = |1 0>
 *     px = (|1 -1> - |1 1>)/sqrt2
 *     py = (|1 -1> + |1 1>)i/sqrt2
 *   l=2
 *     dz2 = |2 0>
 *     dzx = (|2 -1> - |2 1>)/sqrt2
 *     dzy = (|2 -1> + |2 1>)i/sqrt2
 *     dx2 = (|2 -2> + |2 2>)/sqrt2
 *     dxy = (|2 -2> - |2 2>)i/sqrt2
 *   l=3
 *     fz3 = |3 0>
 *     fxz2 = (|3 -1> - |3 1>)/sqrt2
 *     fyz2 = (|3 -1> + |3 1>)i/sqrt2
 *     fzx2 = (|3 -2> + |3 2>)/sqrt2
 *     fxyz = (|3 -2> - |3 2>)i/sqrt2
 *     fx3  = (|3 -3> - |3 3>)/sqrt2
 *     fy3  = (|3 -3> + |3 3>)i/sqrt2
 */
  const int N=2*l+1;
  dcomplex * tt;

  if (Ylm2C!=NULL) {
    tt=Ylm2C;
  }
  else {
    tt=(dcomplex *)malloc(sizeof(dcomplex)*N*N);
  }
  memset(tt, 0, sizeof(dcomplex)*N*N);

  switch(l) {
    case 0:
      tt[0]=1;                    // s
      break;
    case 1:
      tt[1]=1;                    // pz
      tt[1*N+0]=1/sqrt2; tt[1*N+2]=-1/sqrt2;        // px
      tt[2*N+0]=I/sqrt2; tt[2*N+2]=I/sqrt2;     // py
      break;
    case 2:
      tt[2]=1;                   // dz2
      tt[1*N+1]=1/sqrt2; tt[1*N+3]=-1/sqrt2;      // dzx
      tt[2*N+1]=I/sqrt2; tt[2*N+3]=I/sqrt2;   // dzy
      tt[3*N+0]=1/sqrt2; tt[3*N+4]=1/sqrt2;     // dx2-y2
      tt[4*N+0]=I/sqrt2; tt[4*N+4]=-I/sqrt2;    // dxy
      break;
    case 3:
      tt[3]=1;                   // fz3
      tt[1*N+2]=1/sqrt2; tt[1*N+4]=-1/sqrt2;      // fxz2
      tt[2*N+2]=I/sqrt2; tt[2*N+4]=I/sqrt2;   // fyz2
      tt[3*N+1]=1/sqrt2; tt[3*N+5]=1/sqrt2;   // fz(x2-y2)
      tt[4*N+1]=I/sqrt2; tt[4*N+5]=-I/sqrt2;  // fxyz
      tt[5*N+0]=1/sqrt2; tt[5*N+6]=-1/sqrt2;    // fx3-3xy2
      tt[6*N+0]=I/sqrt2; tt[6*N+6]=I/sqrt2;     // fy(3x2-y2)
      break;
  }

#ifdef __DEBUG
    printf("Ylm2C:(original)\n");
    for(int i=0; i<N; i++) {
      for(int j=0; j<N; j++) {
        printf("(%14.9f,%14.9f) ", creal(Ylm2C[i*N+j]), cimag(Ylm2C[i*N+j]));
      }
      printf("\n");
    }
#endif

  if (C2Ylm!=NULL) {
    int * ipiv;
    ipiv=(int *)malloc(sizeof(int)*N);

    memcpy(C2Ylm, tt, sizeof(dcomplex)*N*N);
    LAPACKE_zgetrf(LAPACK_ROW_MAJOR, N, N, C2Ylm, N, ipiv);
    LAPACKE_zgetri(LAPACK_ROW_MAJOR, N, C2Ylm, N, ipiv);

    free(ipiv);
#ifdef __DEBUG
    printf("C2Ylm:(original)\n");
    for(int i=0; i<N; i++) {
      for(int j=0; j<N; j++) {
        printf("(%14.9f,%14.9f) ", creal(C2Ylm[i*N+j]), cimag(C2Ylm[i*N+j]));
      }
      printf("\n");
    }
#endif
  }

  if (Ylm2C==NULL) free(tt);
}

void rotate_cubic(dcomplex * rot, int l, double axis[3], double alpha, int inv) {
/*
 * This function generates the rotation matrix for Cubic harmonics with specific l
 */
  const int N=2*l+1;
  const dcomplex one=1.0;
  const dcomplex zero=0.0;

  dcomplex * C2Ylm;
  dcomplex * rot_Ylm;
  dcomplex * tmp;
  dcomplex * Ylm2C;

  C2Ylm=(dcomplex *)malloc(sizeof(dcomplex)*N*N);
  rot_Ylm=(dcomplex *)malloc(sizeof(dcomplex)*N*N);
  tmp=(dcomplex *)malloc(sizeof(dcomplex)*N*N);
  Ylm2C=(dcomplex *)malloc(sizeof(dcomplex)*N*N);

  generate_C2Ylm(C2Ylm, Ylm2C, l);
#ifdef __DEBUG
    printf("Ylm2C:\n");
    for(int i=0; i<N; i++) {
      for(int j=0; j<N; j++) {
        printf("(%14.9f,%14.9f) ", creal(Ylm2C[i*N+j]), cimag(Ylm2C[i*N+j]));
      }
      printf("\n");
    }
#endif

  rotate_Ylm(rot_Ylm, l, axis, alpha, inv);

#ifdef __DEBUG
    printf("rot_Ylm:\n");
    for(int i=0; i<N; i++) {
      for(int j=0; j<N; j++) {
        printf("(%14.9f,%14.9f) ", creal(rot_Ylm[i*N+j]), cimag(rot_Ylm[i*N+j]));
      }
      printf("\n");
    }
#endif

  // D_cubic = C2Ylm^T · D_sphere · Ylm2c^T
  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasTrans, N, N, N, &one, rot_Ylm, N, Ylm2C, N, &zero, tmp, N);
  cblas_zgemm(CblasRowMajor, CblasTrans, CblasNoTrans, N, N, N, &one, C2Ylm, N, tmp, N, &zero, rot, N);

  //test:
  //cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, &one, C2Ylm , N, rot_Ylm, N, &zero, tmp, N);
  //cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, &one, tmp, N, Ylm2C, N, &zero, rot, N);


  free(Ylm2C);
  free(tmp);
  free(rot_Ylm);
  free(C2Ylm);

}
