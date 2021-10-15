#include "matrix.h"

void matrix3x3_transpose(double out[3][3], double in[3][3]){
    int i,j;
    double tmp[3][3];
    for(i=0;i<3;i++)
        for(j=0;j<3;j++)
            tmp[i][j]=in[j][i];
    for(i=0;i<3;i++)
        for(j=0;j<3;j++)
            out[i][j] = tmp[i][j];
}

void matrix3x3_inverse(double out[3][3], double in[3][3]){
    double * tmp;
    int * ipiv;
    int i,j;
    int info;
    int dim=3;

    tmp = (double *)malloc(sizeof(double)*dim*dim);
    ipiv = (int *)malloc(sizeof(int)*dim);
    for(i=0;i<3;i++)
        for(j=0;j<3;j++)
            tmp[i*3+j]=in[i][j];

    LAPACKE_dgetrf(LAPACK_ROW_MAJOR, 3, 3, tmp, 3, ipiv);
    LAPACKE_dgetri(LAPACK_ROW_MAJOR, 3,    tmp, 3, ipiv);
    for(i=0;i<3;i++)
        for(j=0;j<3;j++)
            out[i][j] = tmp[i*3+j];
    free(tmp);
    free(ipiv);
}

void matrix3x3_dot(double out[3][3], double a[3][3], double b[3][3]){
    double tmpa[9];
    double tmpb[9];
    double tmpc[9];
    int i,j;
    int ipiv[3];
    int info;
    const double one=1.0;
    const double zero=0.0;
    for(i=0;i<3;i++)
        for(j=0;j<3;j++){
            tmpa[3*i+j] = a[i][j];
            tmpb[3*i+j] = b[i][j];
        }
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3,  one, tmpa, 3, tmpb, 3, zero, tmpc, 3);
    for(i=0;i<3;i++)
        for(j=0;j<3;j++)
            out[i][j] = tmpc[3*i+j];
}

double matrix3x3_determinant(double m[3][3])
{
    int i;
    double determinant=0;
    for(i = 0; i < 3; i++){
        determinant += (m[0][i] * (m[1][(i+1)%3] * m[2][(i+2)%3] - m[1][(i+2)%3] * m[2][(i+1)%3]));
    }
    return determinant;
}

void matrix3x3_copy(double out[3][3], double in[3][3]){
    int i,j;
    for(i=0;i<3;i++)
        for(j=0;j<3;j++)
            out[i][j] = in[i][j];
}

void matrix3x1_copy(double out[3], double in[3]){
    int i;
    for(i=0;i<3;i++)
        out[i] = in[i];
}