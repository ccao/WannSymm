#include "rotate_spinor.h"

void rotate_spinor(dcomplex * rot, double axis[3], double alpha, int inv) {
  /* Rotate spinor
   * INPUT:
   *   axis[3]: rotation axis (unit)
   *   alpha  : rotation angle
   *   inv    : inversion flag?
   * OUTPUT:
   *   rot    : 2x2 SU(2) rotation matrix
   */

  const dcomplex one=1.0;
  const dcomplex zero=0.0;

  //dcomplex pauli_sigma[3][4]={{0, 1, 1, 0}, {0, I, -I, 0}, {1, 0, 0, -1}};	// Pauli Matrix...
  dcomplex pauli_sigma[3][4]={{0, 1, 1, 0}, {0, -I, I, 0}, {1, 0, 0, -1}};	// Pauli Matrix...
  /*           [ 0  1 ]           [ 0 -i ]           [ 1  0 ]
   *  sigma_x= [ 1  0 ], sigma_y= [ i  0 ], sigma_z= [ 0 -1 ]
   */

  dcomplex sigma_n[4]={0, 0, 0, 0};
  dcomplex vl[4], vr[4];
  dcomplex eig[2];

  int info;
  int i;

  for(i=0; i<4; i++) {
    sigma_n[i] = -0.5*(pauli_sigma[0][i]*axis[0] + 
                       pauli_sigma[1][i]*axis[1] +
                       pauli_sigma[2][i]*axis[2]) * alpha * I;
  }

  info = LAPACKE_zgeev( LAPACK_ROW_MAJOR, 'N', 'V', 2, sigma_n, 2, eig, vl, 2, vr, 2);

  /*
   * A = X M X^-1
   * exp( M ) = X^-1 exp(A) X
   */

  // rot = exp(eig)
  rot[0]=cexp(eig[0]); rot[1]=0;
  rot[2]=0; rot[3]=cexp(eig[1]);
  // sigma_n=rot * vr^{-1}
  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, 2, 2, 2, &one, rot, 2, vr, 2, &zero, sigma_n, 2);
  // rot = vr * sigma_n
  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 2, 2, 2, &one, vr, 2, sigma_n, 2, &zero, rot, 2);

  //transverse rot, from one that rotating a spin orbital to one that rotating basis of spin orbital
  //dcomplex tmptmp;
  //tmptmp=rot[1];
  //rot[1]=rot[2];
  //rot[2]=tmptmp;


  //if (inv) {
  //  for (i=0; i<4; i++) {
  //    rot[i] = -I*rot[i];
  //  }
  //}
  //for(i=0;i<4;i++){
  //    rot[i] = conj(rot[i]);
  //}
}
