#include "rotate_orbital.h"
#include <stdlib.h>

//#define __DEBUG_rot_orb_detail
//#define __DEBUG_rot_and_Ylm2C

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
#ifdef __DEBUG_rot_orb_detail
    printf("Ylm2C:\n");
    for(int i=0; i<N; i++) {
      for(int j=0; j<N; j++) {
        printf("(%14.9f,%14.9f) ", creal(Ylm2C[i*N+j]), cimag(Ylm2C[i*N+j]));
      }
      printf("\n");
    }
#endif

  rotate_Ylm(rot_Ylm, l, axis, alpha, inv);

#ifdef __DEBUG_rot_orb_detail
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

  #ifdef __DEBUG_rot_and_Ylm2C
    if(l==1 && rot[3*1+1] != 1){//tmp, debug
        printf("l=1\n");
        printf("rot:\n");
        for(int i=0; i<N; i++) {
          for(int j=0; j<N; j++) {
            printf("(%14.9f,%14.9f) ", creal(rot[i*N+j]), cimag(rot[i*N+j]));
          }
          printf("\n");
        }
        printf("Ylm2C:\n");
        for(int i=0; i<N; i++) {
          for(int j=0; j<N; j++) {
            printf("(%14.9f,%14.9f) ", creal(Ylm2C[i*N+j]), cimag(Ylm2C[i*N+j]));
          }
          printf("\n");
        }
        printf("rot_Ylm:\n");
        for(int i=0; i<N; i++) {
          for(int j=0; j<N; j++) {
            printf("(%14.9f,%14.9f) ", creal(rot_Ylm[i*N+j]), cimag(rot_Ylm[i*N+j]));
          }
          printf("\n");
        }
        exit(1);
    }
  #endif

  free(Ylm2C);
  free(tmp);
  free(rot_Ylm);
  free(C2Ylm);

}

void print_orbital_definition_for_JM(char * seed, wannorb_JM * orb_info, int norb){
    char fn_orb_JM[MAXLEN];
    sprintf(fn_orb_JM, "%s_JM_orbital_def.dat", seed);
    FILE *fout = fopen(fn_orb_JM, "w");
    if (fout == NULL) {
        perror(fn_orb_JM);
        exit(EXIT_FAILURE);
    }
    fprintf(fout,"Orbital definition in |J,M> basis:\n");
    fprintf(fout,"# No.    Orbitals' position(Direct)   l    J    M  r         z-axis               x-axis             y-axis\n");
    for(int i=0;i<norb;i++){
        fprintf(fout,"%5d %10.5lf%10.5lf%10.5lf%3d%5.1lf%5.1lf%3d%7.3lf%7.3lf%7.3lf%7.3lf%7.3lf%7.3lf%7.3lf%7.3lf%7.3lf\n", i+1,
                         (orb_info+i)->site.x, (orb_info+i)->site.y, (orb_info+i)->site.z,
                         (orb_info+i)->l, (orb_info+i)->J, (orb_info+i)->M, (orb_info+i)->r,
                         (orb_info+i)->axis[2].x, (orb_info+i)->axis[2].y, (orb_info+i)->axis[2].z,
                         (orb_info+i)->axis[0].x, (orb_info+i)->axis[0].y, (orb_info+i)->axis[0].z,
                         (orb_info+i)->axis[1].x, (orb_info+i)->axis[1].y, (orb_info+i)->axis[1].z);
    }
    fprintf(fout,"\n");
    fclose(fout);
    return;
}

void generate_CS2JM(dcomplex * CS2JM, dcomplex * JM2CS, int l) {
/*
 * Generate the transform matrix between cubic harmonics + spinor and JM.
 * Definition:
 *   CS2YlmS = [ C2Ylm , 0; 0 C2Ylm ]
 *   CS2JM = YlmS2JM * CS2YlmS
 *
 * index for cubic harmonics + spinor:
 *   mr runs in 0,...,2*l+1
 *   ms runs in 0,1
 *   idx_CS = mr + ms*(2*l+1)
 * index for JM:
 *   M runs in -J, -J+1,...,J
 *   if J=l-1/2
 *     idx_JM = M + J
 *   if J=l+1/2
 *     idx_JM = 2*l + M + J
 *
 */

    const int N=(2*l+1)*2;
    const int norb=2*l+1;
    const dcomplex one=1.0;
    const dcomplex zero=0.0;
    dcomplex * tt;
    dcomplex * YlmS2JM;
    dcomplex * C2Ylm;
    dcomplex * CS2YlmS;
    int ii;

    if (CS2JM != NULL) {
        tt = CS2JM;
    } else {
        tt = (dcomplex *)malloc(sizeof(dcomplex)*N*N);
    }
    memset(tt, 0, sizeof(dcomplex)*N*N);
    YlmS2JM = (dcomplex *)malloc(sizeof(dcomplex)*N*N);
    generate_YlmS2JM(YlmS2JM, NULL, l);
    C2Ylm   = (dcomplex *)malloc(sizeof(dcomplex)*N*N/4);
    generate_C2Ylm(C2Ylm, NULL, l);
    CS2YlmS = (dcomplex *)malloc(sizeof(dcomplex)*N*N);
    memset(CS2YlmS, 0, sizeof(dcomplex)*N*N);

    // CS2YlmS = [ C2Ylm, 0; 0, C2Ylm ]
    for (ii=0; ii<N*N; ii++) {
        int row=ii/N;
        int col=ii%N;
        int row_spin=row/norb;
        int col_spin=col/norb;
        if (row_spin == col_spin) {
            CS2YlmS[ii]=C2Ylm[(row%norb)*norb+(col%norb)];
        }
    }

    // Compose the basis transforms: CS2JM = YlmS2JM * CS2YlmS.
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                N, N, N, &one, YlmS2JM, N, CS2YlmS, N,
                &zero, tt, N);

    if (JM2CS != NULL) {
        int * ipiv;
        ipiv=(int *)malloc(sizeof(int)*N);

        memcpy(JM2CS, tt, sizeof(dcomplex)*N*N);
        LAPACKE_zgetrf(LAPACK_ROW_MAJOR, N, N, JM2CS, N, ipiv);
        LAPACKE_zgetri(LAPACK_ROW_MAJOR, N, JM2CS, N, ipiv);

        free(ipiv);
    }
    
    if (CS2JM == NULL) {
        free(tt);
    }
    free(YlmS2JM);
    free(C2Ylm);
    free(CS2YlmS);
}

void generate_YlmS2JM(dcomplex * YlmS2JM, dcomplex * JM2YlmS, int l) {
/* 
 * Generate the transform matrix between Ylm + spinor and |JM> (Clebsch–Gordan coefficients for L+S)
 * Definition:
 *   < l, m_l=(M-1/2), s=1/2, m_s=+1/2| J=(l±1/2), M> = ±sqrt(1/2(1±(M/(l+1/2))))
 *   < l, m_l=(M+1/2), s=1/2, m_s=-1/2| J=(l±1/2), M> =  sqrt(1/2(1-1*±(M/(l+1/2))))
 * l=0, s=1/2, J=1/2
 *   J=1/2
 *      |J=1/2, M=±1/2> = |l=0, m_l=0, s=1/2, m_s=±1/2> 
 * l=1, s=1/2, J=1/2,3/2
 *   J=1/2
 *     |J=1/2, M=±1/2> = -1*±sqrt(1/3) |l=1, m_l=0, s=1/2, m_s=±1/2> ± sqrt(2/3) | l=1, m_l=±1, s=1/2, m_s=-1*±1/2>
 *   J=3/2
 *     |J=3/2, M=±3/2> = |l=1, m_l=±1, s=1/2, m_s=±1/2> 
 *     |J=3/2, M=±1/2> = sqrt(2/3) |l=1, m_l=0, s=1/2, m_s=±1/2> + sqrt(1/3) |l=1, m_l=±1, s=1/2, m_s=-1*±1/2>
 * l=2, s=1/2, J=3/2,5/2
 *   J=3/2
 *     |J=3/2, M=±3/2> = -1*±sqrt(1/5) |l=2, m_l=±1, s=1/2, m_s=±1/2> ± sqrt(4/5) | l=2, m_l=±2, s=1/2, m_s=-1*±1/2>
 *     |J=3/2, M=±1/2> = -1*±sqrt(2/5) |l=2, m_l=0,  s=1/2, m_s=±1/2> ± sqrt(3/5) | l=2, m_l=±1, s=1/2, m_s=-1*±1/2>
 *   J=5/2
 *     |J=5/2, M=±5/2> = |l=2, m_l=±2, s=1/2, m_s=±1/2> 
 *     |J=5/2, M=±3/2> = sqrt(4/5) |l=2, m_l=±1, s=1/2, m_s=±1/2> + sqrt(1/5) |l=2, m_l=±2, s=1/2, m_s=-1*±1/2>
 *     |J=5/2, M=±1/2> = sqrt(3/5) |l=2, m_l=0,  s=1/2, m_s=±1/2> + sqrt(2/5) |l=2, m_l=±1, s=1/2, m_s=-1*±1/2>
 * l=3, s=1/2, J=5/2,7/2
 *   J=5/2
 *     |J=5/2, M=±5/2> = -1*±sqrt(1/7) |l=3, m_l=±2, s=1/2, m_s=±1/2> ± sqrt(6/7) | l=3, m_l=±3, s=1/2, m_s=-1*±1/2>
 *     |J=5/2, M=±3/2> = -1*±sqrt(2/7) |l=3, m_l=±1, s=1/2, m_s=±1/2> ± sqrt(5/7) | l=3, m_l=±2, s=1/2, m_s=-1*±1/2>
 *     |J=5/2, M=±1/2> = -1*±sqrt(3/7) |l=3, m_l=0,  s=1/2, m_s=±1/2> ± sqrt(4/7) | l=3, m_l=±1, s=1/2, m_s=-1*±1/2>
 *   J=7/2
 *     |J=7/2, M=±7/2> = |l=3, m_l=±3, s=1/2, m_s=±1/2> 
 *     |J=7/2, M=±5/2> = sqrt(6/7) |l=3, m_l=±2, s=1/2, m_s=±1/2> + sqrt(1/7) |l=3, m_l=±3, s=1/2, m_s=-1*±1/2>
 *     |J=7/2, M=±3/2> = sqrt(5/7) |l=3, m_l=±1, s=1/2, m_s=±1/2> + sqrt(2/7) |l=3, m_l=±2, s=1/2, m_s=-1*±1/2>
 *     |J=7/2, M=±1/2> = sqrt(4/7) |l=3, m_l=0,  s=1/2, m_s=±1/2> + sqrt(3/7) |l=3, m_l=±1, s=1/2, m_s=-1*±1/2>
 *     
 * Definition of index:
 * index for YlmS ->
 *   ml = m_l + l, can be 0,...,2*l
 *   ms = m_s + s = ±1/2 + 1/2, can be 0,1
 *   idx_YlmS = ml + ms*(2*l+1) = m_l + l + (m_s+1/2) * (2*l+1)
 * index for J, M ->
 *   if J=l-1/2, M=-J,...,J, 
 *      idx_JM = M+J
 *   if J=l+1/2, M=-J,...,J,
 *      idx_JM = 2*l + M + J
 */
    const int N=(2*l+1)*2;
    dcomplex * tt;
    int idx_JM, idx_YlmS;

    if (YlmS2JM != NULL) {
        tt = YlmS2JM;
    } else {
        tt = (dcomplex *)malloc(sizeof(dcomplex)*N*N);
    }
    memset(tt, 0, sizeof(dcomplex)*N*N);
    
    for (idx_JM=0; idx_JM < N; idx_JM++) {
        double J = idx_JM < 2*l ? l-0.5 : l+0.5;
        double M = idx_JM < 2*l ? idx_JM - J : idx_JM - 2*l -J;
        double p_or_m = J > l ? 1 : -1; // plus or minus
        for (idx_YlmS=0; idx_YlmS<N; idx_YlmS++) {
            int ml = idx_YlmS % (2*l+1);
            int ms = idx_YlmS / (2*l+1);
            int m_l = ml - l;
            if (fabs(m_l - (M-0.5)) < 1e-8 && ms == 1){
                // m_l == M-1/2, m_s == 1/2
                tt[idx_JM*N+idx_YlmS] = p_or_m * sqrt((1.0 + p_or_m * (M/(l+0.5)))/2);
            } 
            else if (fabs(m_l - (M+0.5)) < 1e-8 && ms == 0) {
                // m_l == M+1/2, m_s == -1/2
                tt[idx_JM*N+idx_YlmS] = sqrt((1.0 - p_or_m * (M/(l+0.5)))/2);
            }
        }
    }

    if (JM2YlmS != NULL) {
        int * ipiv;
        ipiv=(int *)malloc(sizeof(int)*N);

        memcpy(JM2YlmS, tt, sizeof(dcomplex)*N*N);
        LAPACKE_zgetrf(LAPACK_ROW_MAJOR, N, N, JM2YlmS, N, ipiv);
        LAPACKE_zgetri(LAPACK_ROW_MAJOR, N, JM2YlmS, N, ipiv);

        free(ipiv);
    }

    if (YlmS2JM == NULL) {
        free(tt);
    }
}
