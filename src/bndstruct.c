#include "bndstruct.h"

void bnd_eigcha(double * eig_hk, int * ndegen, dcomplex *** p2sym_chas, dcomplex *** p2sym_eigs, smallgroup * sgrp,
                double lattice[3][3], double rotations[][3][3], double translations[][3], double rots_kd[][3][3], int nsymm,
                wannorb * orb_info, int flag_soc, 
                wanndata * hr, vector kpt, double degenerate_tolerance) {
    // determine the small group (sgrp) that keep kpt unchanged.
    // calculate the band structure (eig_hk) from hr.
    // with eig_hk, derive the number of degeneration (ndegen) of every state
    // calculate symmetry operator  (sym_op) of the given k-point,
    // then derive characters (sym_chas) and eigenvalue (sym_eigs) from sym_op.

    dcomplex * vr_hk;
    dcomplex coeff, scal;
    int io, jo, iio, jjo;
    int irpt, isymm, sgrpi;
    int norb;
    dcomplex * sym_op;         //symmetries' operation
    dcomplex * sym_optmp;
    dcomplex * sym_op_bd;      //bd is abbreviation of block diagonalized
    dcomplex one=1, zero=0;
    vector kpt_roted;
    dcomplex * sym_op_eig;
    dcomplex * sym_op_block;
    dcomplex * vl_sym, * vr_sym;

    norb=hr->norb;


    // get eigenvalues (eig_hk) and eigenvectors (vr_hk)
    vr_hk  = (dcomplex *)malloc(sizeof(dcomplex)*norb*norb);
    bnd_eig_hk(eig_hk, vr_hk, hr, kpt);

    int      ndegentmp = 0;
    // determine number of degenerate states
    for(io=0; io<norb; io++){
        ndegentmp++;
        if(io == norb-1 || fabs(eig_hk[io] - eig_hk[io+1]) > degenerate_tolerance ){
            for(jo=0;jo<ndegentmp;jo++){
                ndegen[io-jo]   = ndegentmp;
            }
            ndegentmp= 0;
        }
    }

    // alloc memory for sym_optmp, sym_op, sym_op_bd
    sym_optmp = (dcomplex *)malloc(sizeof(dcomplex)*norb*norb);
    sym_op    = (dcomplex *)malloc(sizeof(dcomplex)*norb*norb);
    sym_op_bd = (dcomplex *)malloc(sizeof(dcomplex)*norb*norb);

    sgrp->order = 0;
    for(isymm=0;isymm<nsymm;isymm++){
        // the operations keep kpt unchanged sum up to be the small group.
        // we need to calculate the character only for small group.
        kpt_roted = vector_rotate(kpt, rots_kd[isymm]);
        if(kpt_equivalent(kpt_roted, kpt, lattice)){
            sgrp->element[sgrp->order] = isymm;
            // sgrp->element range [0, nsymm-1] for single-valued group
            sgrp->order++;
        }
    }
    if(flag_soc == 1){ // for double group
        sgrp->order *= 2;
        for(sgrpi=0; sgrpi < sgrp->order/2; sgrpi++ ){
            sgrp->element[sgrp->order/2 + sgrpi] = -1 * sgrp->element[sgrpi] - 1; 
            // sgrp->element range [-nsymm, -1] for symmetries introduced by double group
        }
    }

    // alloc memory for sym_chas, sym_eigs 
    *p2sym_chas  = (dcomplex **)malloc(sizeof(dcomplex *)*sgrp->order);
    *p2sym_eigs  = (dcomplex **)malloc(sizeof(dcomplex *)*sgrp->order);
    for(sgrpi=0; sgrpi < sgrp->order; sgrpi++){
        (*p2sym_eigs)[sgrpi] = (dcomplex *)malloc(sizeof(dcomplex)*norb);
        (*p2sym_chas)[sgrpi] = (dcomplex *)malloc(sizeof(dcomplex)*norb);
    }

    dcomplex trace     = 0;
    for(sgrpi=0; sgrpi < sgrp->order; sgrpi++){
        isymm = sgrp->element[sgrpi];
        if(isymm < 0)  isymm = -isymm - 1;
        get_sym_op_reciprocalspace(sym_op, lattice, orb_info, norb, sgrp->element[sgrpi], rotations[isymm], translations[isymm], rots_kd[isymm], kpt, flag_soc);
        
        //sym_op will be block diagonalized by vr_hk 
        //sym_op_block_diag = vr_hk.conj.transe * sym_op * vr_hk
        cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
                    norb, norb, norb, 
                    &one, sym_op, norb, vr_hk, norb, &zero,
                    sym_optmp, norb);
        cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans, 
                    norb, norb, norb, 
                    &one, vr_hk, norb, sym_optmp, norb, &zero,
                    sym_op_bd, norb);

        // calculate the characters from sym_op_bd for every group of degenerate 
        // state corresponding to $kpt.
        for(io=0; io<norb; io+=ndegen[io]){
            trace = 0;
            for(jo=0; jo<ndegen[io]; jo++){
                trace += sym_op_bd[(io+jo)*norb + (io+jo)];
            }
            for(jo=0; jo<ndegen[io]; jo++){
                (*p2sym_chas)[sgrpi][io+jo] = trace;
            }
        }

        for(io=0; io<norb; io+=ndegen[io]){
            sym_op_eig = (*p2sym_eigs)[sgrpi]+io;
            if(ndegen[io] == 1){
                *(sym_op_eig) = sym_op_bd[io*norb+io];
            }
            else{
                vl_sym       = (dcomplex *)malloc(sizeof(dcomplex)*ndegen[io]*ndegen[io]);
                vr_sym       = (dcomplex *)malloc(sizeof(dcomplex)*ndegen[io]*ndegen[io]);
                sym_op_block = (dcomplex *)malloc(sizeof(dcomplex)*ndegen[io]*ndegen[io]);
                for(iio=0;iio<ndegen[io];iio++){
                    for(jjo=0;jjo<ndegen[io];jjo++){
                        sym_op_block[iio*ndegen[io]+jjo] = sym_op_bd[(io+iio)*norb+io+jjo];
                    }
                }
                LAPACKE_zgeev( LAPACK_ROW_MAJOR, 'N', 'V', ndegen[io], sym_op_block,
                               ndegen[io], sym_op_eig, vl_sym, ndegen[io], vr_sym, ndegen[io]);
                free(vl_sym);
                free(vr_sym);
                free(sym_op_block);
            }
        }
    }

    //free the memory allocated this step
    free(vr_hk);
    free(sym_optmp);
    free(sym_op);
    free(sym_op_bd);
}


void bnd_eig_hk(double * eig_hk, dcomplex * vr_hk, wanndata * hr, vector kpt){
    //
    // NOTE: eig_hk and vr_hk must be allocateed before use of this function.
    //
    // Input:  hr, kpt
    // Output: eig_hk, vr_hk
    //
    // 1. calculate recciprocal-space Hamiltonian (hamk) corresponding to the given k-point (kpt)
    // from the real-space Hamiltonian (hr).
    // 2. diagonalize hamk and output the eigenvalues (eig_hk) and eigenvectors (vr_hk)

    int irpt, i;
    dcomplex coeff, scal;
    dcomplex * ham_k;

    int norb=hr->norb;

    // get ham_k from hr
    ham_k = (dcomplex *)malloc(sizeof(dcomplex)*norb*norb);
    for(i=0; i<norb*norb; i++){
        ham_k[i] = 0;
    }
    for(irpt=0;irpt<hr->nrpt;irpt++){
        coeff = cexp( cmplx_i * 2 * PI * dot_product( kpt, hr->rvec[irpt] ) );
        scal  = coeff / (dcomplex) hr->weight[irpt];
        // cblas_zaxpy(n, a, x, 1, y, 1) means: y = a*x + y
        cblas_zaxpy(norb*norb, &scal, hr->ham + irpt*norb*norb, 1, ham_k, 1);
    }
    
    // diagonalize hamk
    LAPACKE_zheev( LAPACK_ROW_MAJOR,'V',  'U', 
                   norb, ham_k, norb, eig_hk);
    // zheev changes ham_k to its right eigenvectors.
    if(vr_hk != NULL){
        memcpy(vr_hk, ham_k, sizeof(dcomplex)*norb*norb);
    }
    free(ham_k);
}

void bnd_write_characters(const char * fnout, int ikpt, vector kpt, smallgroup * sgrp, 
                  double lattice[3][3], double rotations[][3][3], double translations[][3], int TR[],
                  double * eig_hk, int * ndegen, dcomplex ** sym_chas, int norb, int flag_soc,
                  int en_print_prec, int bnd_print_len){

    FILE * fbndsymcha;
    int sgrpi, isymm;
    int io;
    char msg[MAXLEN];
    char msg1[MAXLEN];

    fbndsymcha = fopen(fnout, "a");
    //output this kpt
    //fprintf(fbndsymcha, "kpt No. %5d:%10.7lf %10.7lf %10.7lf\n", ikpt+1, kpt.x, kpt.y, kpt.z);
    fprintf(fbndsymcha, "kpt:%10.7lf %10.7lf %10.7lf\n", kpt.x, kpt.y, kpt.z);
    //output symmetries that keep kpt invariant
    fprintf(fbndsymcha, "Related symmetries:\n");
    fclose(fbndsymcha);
    for(sgrpi=0; sgrpi < sgrp->order; sgrpi++){
        isymm = sgrp->element[sgrpi];
        if(isymm < 0)  isymm = -isymm - 1;
        print_symmetry("bnd_sym_characters", lattice, sgrp->element[sgrpi], rotations[isymm], translations[isymm], TR[isymm], 0, -1, flag_soc);
    }
    fbndsymcha = fopen("bnd_sym_characters", "a");
    //output characters of symm each state corresponding to this kpt
    for(io=0; io<norb; io+=ndegen[io]){
        sprintf(msg1, "%*d - %-*d", bnd_print_len , io+1, bnd_print_len, io+ndegen[io]);
        sprintf(msg, "band No. %s  energy =%*.*lf eV  ndegen = %*d", 
                        msg1, en_print_prec+5, en_print_prec, eig_hk[io], bnd_print_len, ndegen[io]);
        fprintf(fbndsymcha, "%s\n    ", msg);
        for(sgrpi=0; sgrpi<sgrp->order; sgrpi++){
            if(creal(sym_chas[sgrpi][io]) < 0 && fabs(creal(sym_chas[sgrpi][io])) < eps6){
                // display -0.000 as 0.000
                sym_chas[sgrpi][io] = 0 + cimag(sym_chas[sgrpi][io]) * cmplx_i;
            }
            if(cimag(sym_chas[sgrpi][io]) < 0 && fabs(cimag(sym_chas[sgrpi][io])) < eps6){
                // display -0.000 as 0.000
                sym_chas[sgrpi][io] = creal(sym_chas[sgrpi][io]) + 0 * cmplx_i;
            }
            fprintf(fbndsymcha, "%6.3lf%+6.3lfi",  creal(sym_chas[sgrpi][io]), cimag(sym_chas[sgrpi][io]) );
            if(sgrpi%4 == 3 || sgrpi == sgrp->order-1) 
                fprintf(fbndsymcha, "\n    ");
            else
                fprintf(fbndsymcha, ", ");
        }
        fprintf(fbndsymcha, "\n");
    }
    fclose(fbndsymcha);
}

void bnd_write_eigenvalues(const char * fnout, int ikpt, vector kpt, smallgroup * sgrp, 
                  double lattice[3][3], double rotations[][3][3], double translations[][3], int TR[],
                  double * eig_hk, int * ndegen, dcomplex ** sym_eigs, int norb, int flag_soc,
                  int en_print_prec, int bnd_print_len){

    FILE * fbndsymeig;
    int sgrpi, isymm;
    int io, jo;
    char msg[MAXLEN];
    char msg1[MAXLEN];
    
    fbndsymeig = fopen(fnout, "a");
    //output this kpt 
    //fprintf(fbndsymeig, "kpt No. %5d:%10.7lf %10.7lf %10.7lf\n", ikpt+1, kpt.x, kpt.y, kpt.z);
    fprintf(fbndsymeig, "kpt:%10.7lf %10.7lf %10.7lf\n", kpt.x, kpt.y, kpt.z);
    //output symmetries that keep kpt invariant
    fprintf(fbndsymeig, "Related symmetries:\n");
    fclose(fbndsymeig);
    for(sgrpi=0; sgrpi<sgrp->order; sgrpi++){
        isymm = sgrp->element[sgrpi];
        if(isymm < 0)  isymm = -isymm - 1;
        print_symmetry("bnd_sym_eig", lattice, sgrp->element[sgrpi], rotations[isymm], translations[isymm], TR[isymm], 0, -1, flag_soc);
    }
    fbndsymeig = fopen("bnd_sym_eig", "a");
    //output characters of symm each state corresponding to this kpt
    for(io=0;io<norb;io+=ndegen[io]){
        for(jo=0; jo<ndegen[io]; jo++){
            sprintf(msg1, "%*d - %-*d", bnd_print_len, io+1, bnd_print_len, io+ndegen[io]);
            sprintf(msg, "No.%*d of %*d degenerage band (band No. %s) energy =%*.*lf eV ", 
                            bnd_print_len, jo+1, bnd_print_len, ndegen[io], msg1, en_print_prec+5, en_print_prec, eig_hk[io]);
            fprintf(fbndsymeig, "%s\n    ", msg);
            for(sgrpi=0; sgrpi<sgrp->order; sgrpi++){
                if (creal(sym_eigs[sgrpi][io + jo]) < 0 && fabs(creal(sym_eigs[sgrpi][io + jo])) < eps6){
                    // display -0.000 as 0.000
                    sym_eigs[sgrpi][io + jo] = 0 + cimag(sym_eigs[sgrpi][io + jo]) * cmplx_i;
                }
                if (cimag(sym_eigs[sgrpi][io + jo]) < 0 && fabs(cimag(sym_eigs[sgrpi][io + jo])) < eps6){
                    // display -0.000 as 0.000
                    sym_eigs[sgrpi][io + jo] = creal(sym_eigs[sgrpi][io + jo]) + 0 * cmplx_i;
                }
                fprintf(fbndsymeig, "%6.3lf%+6.3lfi",  creal(sym_eigs[sgrpi][io+jo]), cimag(sym_eigs[sgrpi][io+jo]) );
                if(sgrpi%4 == 3 || sgrpi == sgrp->order-1) 
                    fprintf(fbndsymeig, "\n    ");
                else
                    fprintf(fbndsymeig, ", ");
            }
            fprintf(fbndsymeig, "\n");
        }
    }
    fclose(fbndsymeig);
}


void bnd_write_bands(FILE * fbands, int norb, int nkpath, int nk_per_kpath, double lattice[3][3], vector * kvecs, double ** ebands){
    int io, ikpt;
    double kpath_len;
    double rec_latt[3][3];

    // reciprocal lattice = lattice.inv.transpose
    matrix3x3_inverse(rec_latt, lattice);
    matrix3x3_transpose(rec_latt, rec_latt);

    fprintf(fbands, "#   k-path len        Energy\n");
    for(io=0; io<norb; io++){
        kpath_len = 0;
        for(ikpt=0; ikpt<nkpath*nk_per_kpath; ikpt++){
            if(ikpt%nk_per_kpath != 0){
                kpath_len += 2 * PI * vector_norm(vector_rotate(vector_sub(kvecs[ikpt], kvecs[ikpt-1]), rec_latt));
            }
            fprintf(fbands, "%15.9lf %*.*lf", kpath_len, 21, 15, ebands[ikpt][io]);
            //fprintf(fbands, "  # %12.7lf %12.7lf %12.7lf", kvecs[ikpt].x, kvecs[ikpt].y, kvecs[ikpt].z);
            fprintf(fbands, "\n");
            if(ikpt%nk_per_kpath == nk_per_kpath - 1 ){
                fprintf(fbands, "\n");
            }
        }
        fprintf(fbands, "\n");
    }
}