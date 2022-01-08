#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <time.h>


#include "spglib.h"
#include "constants.h"
#include "wanndata.h"
#include "wannorb.h"
#include "vector.h"
#include "matrix.h"
#include "readinput.h"
#include "rotate_ham.h"
#include "rotate_basis.h"
#include "readsymm.h"
#include "version.h"
//#include "para.h"
#include "usefulio.h"
#include "bndstruct.h"

#include "mkl.h"
#include <mpi.h>
#define MPI_INCLUDED


//#define __DEBUG_notexpand
//#define __DEBUG
//#define __DEBUG_bndcha
//#define __DEBUG_sym_eig


int main(int argc, char ** argv){

    MPI_Init(NULL, NULL);                       // Initialize the MPI environment
    int mpi_size;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);   // Get the number of processes
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);   // Get the rank of the process
    MPI_Barrier(MPI_COMM_WORLD);
    //---//init_para();

    int flag_soc=0;                            //Switch : spin orbital coupling (default =0 nonsoc)
    //int flag_dbgrp;                            //Switch : consider double group or not, only make sense with flag_soc == 1
    int flag_global_trsymm=1;                  //switch initialization: global time reversal symmetry
    int flag_expandrvec=1;                     //switch initialization: expand rvec of hr.dat, then all weights are '1'
    int flag_symm_from_file=0;                 //Switch initialization: where symmetries come from? 0: use spglib; 1: use input file
    int flag_everysymm=0;                      //Switch initialization: Output the hr' of hr0 operated by each symmetry operation? 1:yes 0:no
    double ham_tolerance=0.1;                  //tolerance used when flag_expandrvec==0, if some element of 
                                               //symmed_hr.dat is different from the original one, and the 
                                               //difference is larger than ham_tolerance, the code will 
                                               //show a warning.
    double symm_magnetic_tolerance=1E-3;

    char fn_symm[MAXLEN];                       //name of file for user set symmetries

    FILE * fdebug;
    FILE * fstdout;
    char * fn_input;                    //input file name
    char seed[MAXLEN];                          //seed name of wannier90
    char seedup[MAXLEN];                        //seed name for spin up of spin collinear calculation
    char seeddn[MAXLEN];                        //seed name for spin dn of spin collinear calculation
    double lattice[3][3];                       //crystal lattice
    SymmetryOperator S[MAX_NUM_of_SYMM];
    double rotations[MAX_NUM_of_SYMM][3][3];    //symmetry operations -> rotations(double)
    double rots_kd[MAX_NUM_of_SYMM][3][3];      //symmetry rotation matrix in reci lattice
    //int irotations[MAX_NUM_of_SYMM][3][3];    //symmetry operations -> rotations(int)
    double translations[MAX_NUM_of_SYMM][3];    //symmetry operations -> translations
    int TR[MAX_NUM_of_SYMM];                    //symmetry operations -> Time Reversal or not, 1: .true. 0: .false. -1: NULL
    int nsymm;                                  //number of symmetry operations (get from spglib)
    wannorb * orb_info;                         //information for each wannier orbital in one w-s unit cell
    int norb, nwann;                            //number of wannier orbitals in one w-s unit cell

    vector kpt;                                 //Kpoint for calculating bands' character of every symmetry
    int ikpt, nkpt=0;
    int flag_chaeig=1;                     //Switch: calculate bands' character of every symmetry for specified kpts
    int flag_chaeig_in_kpath=0;
    int flag_bands=3;           // 0: do not calculate bands; 1: symmed bands; 2: original bands; 3: symmed and original bands
    vec_llist * kpts;           // for backward compatibility with the 'kpt' tag
    int nk_per_kpath=20, ikpath, nkpath=0;
    vec_llist * kpaths;         // every two vectors make up one k-path
    char klabels[MAXLEN*2][SHORTLEN];

    vec_llist_init(&kpts);
    vec_llist_init(&kpaths);

    wanndata ham_in;                            //hamtonian read from seed_hr.dat

    int i,j,k;                                  
    int irpt, jrpt, isymm;
    int io, jo, mo, no;
    int iio, jjo;
    int iblock;
    int iiham, jjham;
    wanndata * ham_out;                         //array of hamtoinans get by the symmetry operation
    wanndata ham_final;
    wanndata  ham_trsymm;                       // tmp variable for time reversal procedure
    int nrpt_final;
    char ham_out_seed[MAXLEN];

    int flag_debug_rothamk;
    FILE * fbndsymcha;
    FILE * fbndsymeig;
    dcomplex * vl_hk, * vr_hk;
    double   * eig_hk;
    int mark_kpt_inv_under_rot[MAX_NUM_of_SYMM];
    const dcomplex one =1.0;
    const dcomplex zero=0.0;
    int * ndegen;
    dcomplex ** sym_chas, ** sym_eigs;
    double degenerate_tolerance=1E-6;                         // to judge two band is degenerate
    //double translation_kd[3];                  // symmetry's translation in reci lattice
    double lattice_trans[3][3];
    double ATA[3][3];
    double ATA_inv[3][3];
    //double llattice_trans[3][3];
    //double llattice_trans_inv[3][3];
    //vector trans_v, translation_kd_v;

    wanndata hamkr;
    wanndata hamkr1;
    dcomplex * eig_hkr, * eig_hkrtmp, * vr_hkr;
    int dia_err_count;

    int flag_restart=0;           //flag_restart = 0 (start from beginning)
                                //               1 (start from wannier90_symmed_hr.dat)
    int flag_restart_character;
    dcomplex dctmp;

    vector magmom[MAX_NUM_of_atoms];          //magnetic moment
    vector qvector;                 //qvector for anti-ferromagnetism

    char msg[MAXLEN];               //buffer for output message
    char msg1[MAXLEN];               //buffer for output message
    char msg2[MAXLEN];               //buffer for output message
    
    time_t  tm_start, tm_now;
    clock_t clk_start, clk_now;
    double  real_time_elapsed, cpu_time_used;

    int     flag_output_mem_usage=0;
    char    mem_usage_str[MAXLEN];

    tm_start  = time(NULL);
    clk_start = clock();


    if(mpi_rank==0){
        //remove error file wannsymm.err
        remove("wannsymm.err");
    }

    int ver_major, ver_minor, ver_patch, ver_pre_release;
    ver_major=VERSION_MAJOR;
    ver_minor=VERSION_MINOR;
    ver_patch=VERSION_PATCH;
    char version_str[MAXLEN];

    sprintf(version_str, "WannSymm version ");
    sprintf(version_str, "%s%d.%d.%d", version_str, ver_major, ver_minor, ver_patch);
    #ifdef PRE_RELEASE_ALPHA
        sprintf(version_str, "%s-alpha", version_str);
    #elif defined PRE_RELEASE_BETA
        sprintf(version_str, "%s-beta", version_str);
    #elif defined PRE_RELEASE_CANDIDATE
        sprintf(version_str, "%s-rc", version_str);
    #endif
    #if defined(PRE_RELEASE_ALPHA) || defined(PRE_RELEASE_BETA) || defined(PRE_RELEASE_CANDIDATE)
        #if defined(VERSION_PRE)
            ver_pre_release = VERSION_PRE;
            sprintf(version_str, "%s%d", version_str, ver_pre_release);
        #endif
    #endif

    if(mpi_rank == 0){
        if(argc>1){
            if( strcmp(argv[1], "--version") == 0 || strcmp(argv[1], "-V") == 0 || strcmp(argv[1], "-v") == 0 ){
                //print version
                printf("%s\n", version_str);
                exit(0);
            }
            fn_input=argv[1];
            if( ! file_exists(fn_input) ){
                sprintf(msg,"ERROR!!! the input file \"%s\" can not be found\n", fn_input);
                print_error(msg);
                exit(1);
            }
        } else {
            if(file_exists("wannsymm.in")) {
                fn_input="wannsymm.in";
            } else if( file_exists("symmham.in") ){
                fn_input="symmham.in";
            } else {
                sprintf(msg,"ERROR!!! no input file found, default input file name is \"symmham.in\" or \"wannsymm.in\"\n" );
                print_error(msg);
                exit(1);
            }
        }
    }


    //STEP 1 readinput
    //1.1. get the DFT code type, SOC flag and seed name of wannier90
    //1.2. get symmetry operations (in real space) with spglib from crystal info in input file.
    //1.3. get wannier orbital information from input file, which is copy from 
    //     wannier90.wout setting iprint=3 in corresponding wannier90.win
    if(mpi_rank == 0){
        fstdout=fopen("wannsymm.out", "w");
        fprintf(fstdout, "%s\n", version_str);
        fclose(fstdout);
        readinput(fn_input,     
                  &flag_soc,
                  seed,
                  lattice,
                  rotations,
                  translations, 
                  TR, 
                  magmom,
                  &nsymm, &orb_info, &nwann, &kpts, &nkpt, &kpaths, klabels, &nkpath, &nk_per_kpath,
                  &flag_bands, &flag_chaeig, &flag_chaeig_in_kpath, &flag_restart, &flag_global_trsymm, 
                  &flag_expandrvec, &symm_magnetic_tolerance, &ham_tolerance, &degenerate_tolerance,
                  &flag_everysymm, &flag_output_mem_usage,
                  &flag_symm_from_file, fn_symm);
        if( flag_symm_from_file==1)    //read symmetry from symminputfile
            readsymm(fn_symm, rotations, translations, TR, &nsymm, &flag_global_trsymm);
    }

#ifdef __DEBUG_notexpand
    flag_expandrvec=0;
    if(mpi_rank==0){
        printf("WARNING!!! expandrvec=0 : not expanding rvector of hr.dat, the symmed_hr.dat is not in highest symmetry state\n");
        fstdout=fopen("wannsymm.out", "a");
        fprintf(fstdout, "WARNING!!! expandrvec = 0: not expanding rvector of hr.dat, the symmed_hr.dat is not in highest symmetry state\n\n\n");
        fclose(fstdout);
    }
#endif 

    if(mpi_rank == 0){
        fstdout=fopen("wannsymm.out", "a");
        fprintf(fstdout,"Orbital definition:\n");
        fprintf(fstdout,"# No.    Orbitals' position(Direct)   l mr ms  r         z-axis               x-axis             y-axis\n");
        for(i=0;i<nwann;i++){
            fprintf(fstdout,"%5d %10.5lf%10.5lf%10.5lf%3d%3d%3d%3d%7.3lf%7.3lf%7.3lf%7.3lf%7.3lf%7.3lf%7.3lf%7.3lf%7.3lf\n", i+1,
                         (orb_info+i)->site.x, (orb_info+i)->site.y, (orb_info+i)->site.z,
                         (orb_info+i)->l, (orb_info+i)->mr, (orb_info+i)->ms, (orb_info+i)->r,
                         (orb_info+i)->axis[2].x, (orb_info+i)->axis[2].y, (orb_info+i)->axis[2].z,
                         (orb_info+i)->axis[0].x, (orb_info+i)->axis[0].y, (orb_info+i)->axis[0].z,
                         (orb_info+i)->axis[1].x, (orb_info+i)->axis[1].y, (orb_info+i)->axis[1].z);
        }
        fprintf(fstdout,"\n");
        fclose(fstdout);
    }

    double rot_angle, rot_axis[3];
    int inv_flag;
    //double drotation[3][3];
    double complex * orb_rot[5];
    double complex * s_rot;
    int l,mr1,mr2,ii,jj;
    if(mpi_rank ==0){
        //output symmetry axis and angles
        fstdout=fopen("wannsymm.out", "a");
        fprintf(fstdout,"Symmetries' info:\n");
        fclose(fstdout);
        for(isymm=0;isymm<nsymm;isymm++){
            print_symmetry("wannsymm.out", lattice, isymm, rotations[isymm], translations[isymm], TR[isymm], 1, 1, 0);
        }
        if(flag_global_trsymm){
            print_msg("and a global time-reversal symmetry.\n");
        }
        print_msg("\n");
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&flag_soc, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&flag_restart, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&flag_chaeig, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&flag_global_trsymm, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&flag_expandrvec, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&flag_output_mem_usage, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&degenerate_tolerance, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(seed, MAXLEN, MPI_CHAR, 0, MPI_COMM_WORLD);
    for(i=0;i<3;i++)
        for(j=0;j<3;j++)
            MPI_Bcast(&lattice[i][j], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    for(i=0;i<MAX_NUM_of_SYMM;i++)
        for(j=0;j<3;j++)
            for(k=0;k<3;k++)
                MPI_Bcast(&rotations[i][j][k], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    for(i=0;i<MAX_NUM_of_SYMM;i++)
        for(j=0;j<3;j++)
                MPI_Bcast(&translations[i][j], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(TR, MAX_NUM_of_SYMM, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nsymm, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nwann, 1, MPI_INT, 0, MPI_COMM_WORLD);
    norb=nwann;

    if(mpi_rank != 0)
        orb_info=(wannorb *) malloc(sizeof(wannorb)*(nwann));
    MPI_Barrier(MPI_COMM_WORLD);

    int nitems=3;
    int blocklengths[3]={1,1,1};
    MPI_Datatype types_vector[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Datatype mpi_vector;
    MPI_Aint     offsets[3];
    offsets[0] = offsetof(vector, x);
    offsets[1] = offsetof(vector, y);
    offsets[2] = offsetof(vector, z);
    MPI_Type_create_struct(nitems, blocklengths, offsets, types_vector, &mpi_vector);
    MPI_Type_commit(&mpi_vector);

    nitems=6;
    int blocklengths_wannorb[6]={1,3,1,1,1,1};
    MPI_Datatype types_wannorb[6] = {mpi_vector, mpi_vector, MPI_INT, MPI_INT, MPI_INT, MPI_INT};
    MPI_Datatype mpi_wannorb;
    MPI_Aint     offsets_wannorb[6];
    offsets_wannorb[0] = offsetof(wannorb, site);
    offsets_wannorb[1] = offsetof(wannorb, axis);
    offsets_wannorb[2] = offsetof(wannorb, r);
    offsets_wannorb[3] = offsetof(wannorb, l);
    offsets_wannorb[4] = offsetof(wannorb, mr);
    offsets_wannorb[5] = offsetof(wannorb, ms);
    MPI_Type_create_struct(nitems, blocklengths_wannorb, offsets_wannorb, types_wannorb, &mpi_wannorb);
    MPI_Type_commit(&mpi_wannorb);

    MPI_Bcast(orb_info, nwann, mpi_wannorb, 0, MPI_COMM_WORLD);




    //STEP 2 read ham
    //2.1. read seed_hr.dat to ham_in
    if(mpi_rank == 0){
        read_ham(&ham_in, seed);
        if( flag_restart == 1){
            sprintf(ham_out_seed, "%s_symmed", seed);
            read_ham(&ham_final, ham_out_seed);
        }
        if (nwann != ham_in.norb){
            sprintf(msg,"!!! ERROR: Number of wannier orbitals mismatch:\n norb in inputfile %d\n norb in hr.dat %d \n", nwann, ham_in.norb);
            print_error(msg);
            exit(1);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&ham_in.norb, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ham_in.nrpt, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if(mpi_rank != 0)
        init_wanndata(&ham_in);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(ham_in.ham, ham_in.nrpt*ham_in.norb*ham_in.norb, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(ham_in.hamflag, ham_in.nrpt*ham_in.norb*ham_in.norb, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(ham_in.rvec, ham_in.nrpt, mpi_vector, 0, MPI_COMM_WORLD);
    MPI_Bcast(ham_in.weight, ham_in.nrpt, MPI_INT, 0, MPI_COMM_WORLD);

    if( mpi_rank == 0){
        real_time_elapsed = difftime( time(NULL), tm_start);
        cpu_time_used     = ( (double) ( clock() - clk_start) ) / CLOCKS_PER_SEC;
        sprintf(msg, "finished reading input files, real time elasped %.0f s\n\n", real_time_elapsed);
        print_msg(msg);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    //STEP 3 Rotating
    //3.1 get many(nsymm) Hamitonians (ham_out[i]) by the operations of symmetry
    if( flag_restart == 1){     // restarting state, no need to symmetrise it
        //ham_final = ham_in;
        if(mpi_rank != 0) finalize_wanndata(ham_in); 
    }
    else if( flag_restart == 0){
        ham_out = (wanndata *)malloc(sizeof(wanndata)*nsymm);

        int iter;
        int round_total;

        round_total = (int)(nsymm / mpi_size);
        if(nsymm % mpi_size > 0)
            round_total += 1;

        //for(i=0;i<nsymm;i++){
        //    ham_out[i].nrpt=ham_in.nrpt;
        //    ham_out[i].norb=ham_in.norb;
        //    //init_wanndata(ham_out+i);
        //}
        MPI_Barrier(MPI_COMM_WORLD);


        for(iter=0;iter<round_total;iter++){
            i=iter*mpi_size + mpi_rank;
            if( i > nsymm-1){
                continue;
            }
            if(TR[i] == 1){ // ham_out need a procedure of Time reversal
                ham_trsymm.nrpt = ham_in.nrpt;
                ham_trsymm.norb = ham_in.norb;
                init_wanndata(&ham_trsymm);
                trsymm_ham( &ham_trsymm, &ham_in, orb_info, flag_soc);
                rotate_ham(ham_out+i, &ham_trsymm, lattice, rotations[i], translations[i], orb_info, flag_soc);
                finalize_wanndata(ham_trsymm);
            } else{
                rotate_ham(ham_out+i, &ham_in, lattice, rotations[i], translations[i], orb_info, flag_soc);
            }
            real_time_elapsed = difftime( time(NULL), tm_start);
            cpu_time_used     = ( (double) ( clock() - clk_start) ) / CLOCKS_PER_SEC;
            sprintf(msg, "Done: symmetry No.%5d, real time elasped %.0f s, cpu time (thread%3d) %.2f s", i+1, real_time_elapsed, mpi_rank+1, cpu_time_used);
            if( flag_output_mem_usage == 1 ){ 
                sprintf(mem_usage_str, "mem useage (thread%3d) ", mpi_rank+1);
                get_memory_usage_str(mem_usage_str);
                sprintf(msg, "%s, %s", msg, mem_usage_str);
            }
            sprintf(msg, "%s\n", msg);
            print_msg(msg);
        }
        if(mpi_rank != 0) finalize_wanndata(ham_in); 
        MPI_Barrier(MPI_COMM_WORLD);


        int isource;
        for(isymm=0;isymm<nsymm;isymm++){
            isource = isymm%mpi_size; //isource= which thread the ham_out[i] locate
            if(isource != 0){
                // the first thread already has the data of itself (isource=0)
                if(mpi_rank == isource){
                    MPI_Send(&(ham_out[isymm].norb), 1, MPI_INT,       0, 0, MPI_COMM_WORLD);
                    MPI_Send(&(ham_out[isymm].nrpt), 1, MPI_INT,       0, 1, MPI_COMM_WORLD);
                }
                else if(mpi_rank == 0){
                    MPI_Recv(&(ham_out[isymm].norb), 1, MPI_INT, isource, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(&(ham_out[isymm].nrpt), 1, MPI_INT, isource, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    init_wanndata(ham_out + isymm);
                }
                if(mpi_rank == isource){
                    MPI_Send(ham_out[isymm].ham,     ham_out[isymm].nrpt*norb*norb, MPI_DOUBLE_COMPLEX, 0, 2, MPI_COMM_WORLD);
                    MPI_Send(ham_out[isymm].hamflag, ham_out[isymm].nrpt*norb*norb, MPI_INT,            0, 3, MPI_COMM_WORLD);
                    MPI_Send(ham_out[isymm].rvec,   ham_out[isymm].nrpt, mpi_vector, 0, 4, MPI_COMM_WORLD);
                    MPI_Send(ham_out[isymm].weight, ham_out[isymm].nrpt, MPI_INT,    0, 5, MPI_COMM_WORLD);
                }
                else if(mpi_rank == 0){
                    MPI_Recv(ham_out[isymm].ham,     ham_out[isymm].nrpt*norb*norb, MPI_DOUBLE_COMPLEX, isource, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(ham_out[isymm].hamflag, ham_out[isymm].nrpt*norb*norb, MPI_INT,            isource, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(ham_out[isymm].rvec,   ham_out[isymm].nrpt, mpi_vector, isource, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(ham_out[isymm].weight, ham_out[isymm].nrpt, MPI_INT,    isource, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }

                if(mpi_rank == isource) 
                    finalize_wanndata(ham_out[isymm]);
            }
            //MPI_Barrier(MPI_COMM_WORLD);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);


    //3.2 calc the average of the Hams got from (3.1)
    vec_llist * rvecs;
    int nrvec=0;
    vec_llist_init(&rvecs);
    int nerr=1;
    if( mpi_rank == 0 && flag_restart == 0 ){
        //nrpt_final = ham_in.nrpt;
        for(irpt=0; irpt < ham_in.nrpt; irpt++){
            vec_llist_add(&rvecs, ham_in.rvec[irpt]);
            nrvec++;
        }
        if( flag_expandrvec == 1){
            for(i=0;i<nsymm;i++){
                for(j=0;j<ham_out[i].nrpt;j++){
                    if( vec_llist_find(&rvecs, ham_out[i].rvec[j]) == -1 ){
                        //vec_llist_add(&rvecs, ham_out[i].rvec[j]);
                        vec_llist_add_inorder(&rvecs, ham_out[i].rvec[j], 0);
                        nrvec++;
                    }
                }
            }
        }
        ham_final.nrpt=nrvec;
        ham_final.norb=ham_in.norb;
        init_wanndata(&ham_final);
        for( irpt=0; irpt < nrvec; irpt++){
            ham_final.rvec[irpt] = vec_llist_pop(&rvecs, &nerr);
            if( nerr == -1){
                sprintf(msg, "Error occured in finding rotated R vectors.");
                print_error(msg); exit(1);
            }
        }
        if( flag_expandrvec == 0){
            memcpy(ham_final.weight, ham_in.weight, sizeof(int)*ham_in.nrpt);
        }
        //ham_final.weight=1 by default
        for(irpt=0;irpt<ham_final.nrpt;irpt++){
            for(isymm=0;isymm<nsymm;isymm++){
                jrpt = find_vector(ham_final.rvec[irpt], ham_out[isymm].rvec, ham_out[isymm].nrpt);
                if(jrpt == -1) continue;
                for(io=0;io<norb;io++){
                    for(jo=0;jo<norb;jo++){
                        iiham = irpt*norb*norb + io*norb + jo;
                        jjham = jrpt*norb*norb + io*norb + jo;
                        ham_final.ham[iiham]     += ham_out[isymm].ham[jjham];
                        ham_final.hamflag[iiham] += ham_out[isymm].hamflag[jjham];
                    }
                }
            }
            for(io=0;io<norb;io++){
                for(jo=0;jo<norb;jo++){
                    iiham = irpt*norb*norb + io*norb + jo;
                    if(ham_final.hamflag[iiham] != 0 ){
                        ham_final.ham[iiham] /= ham_final.hamflag[iiham];
                        if( flag_expandrvec == 0){
                            ham_final.ham[iiham] *= ham_final.weight[irpt];
                        }
                    }
                }
            }
        } 
        real_time_elapsed = difftime( time(NULL), tm_start);
        cpu_time_used     = ( (double) ( clock() - clk_start) ) / CLOCKS_PER_SEC;
        sprintf(msg, "Done: Average the rotated Hamiltonians. real time elasped %.0f s, cpu time (thread%3d) %.2f s", real_time_elapsed, mpi_rank+1, cpu_time_used);
        if( flag_output_mem_usage == 1 ){ 
            sprintf(mem_usage_str, "mem useage (thread%3d) ", mpi_rank+1);
            get_memory_usage_str(mem_usage_str);
            sprintf(msg, "%s, %s", msg, mem_usage_str);
        }
        sprintf(msg, "%s\n", msg);
        print_msg(msg);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    // Output the hr' derived by hr0 operated with every symmetry operation
    if( mpi_rank == 0 && flag_restart == 0 && flag_everysymm == 1){
        for(i=0;i<nsymm;i++){
            if(nsymm<10-1)
                sprintf(ham_out_seed,"symm%d",i+1);
            else if(nsymm<100-1)
                sprintf(ham_out_seed,"symm%02d",i+1);
            else if(nsymm<1000-1)
                sprintf(ham_out_seed,"symm%03d",i+1);
            else
                sprintf(ham_out_seed,"symm%d",i+1);
            fstdout=fopen("wannsymm.out", "a");
            fprintf(fstdout, "writing %s_hr.dat\n", ham_out_seed);
            fclose(fstdout);
            write_ham(ham_out+i, ham_out_seed);
        }
        //MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // free memory of ham_out[i]
    if(flag_restart==0 && mpi_rank == 0){ //else ham_out not used
        for(i=0;i<nsymm;i++) finalize_wanndata(ham_out[i]);
        free(ham_out);  
    }

    MPI_Barrier(MPI_COMM_WORLD);
    //3.3 use time reversal symmetry to symmetrise ham.final
    if( mpi_rank == 0 && flag_restart==0){
        if(flag_global_trsymm == 1){
            // calculate the ham with time reversal symmetry operation
            ham_trsymm.nrpt = ham_final.nrpt;
            ham_trsymm.norb = ham_final.norb;
            init_wanndata(&ham_trsymm);
            trsymm_ham(&ham_trsymm, &ham_final, orb_info, flag_soc);
            // avearge the ham
            for(i=0;i<ham_final.nrpt*ham_final.norb*ham_final.norb;i++){
                ham_final.ham[i] = (ham_final.ham[i] + ham_trsymm.ham[i])/2;
            }
            finalize_wanndata(ham_trsymm);

            real_time_elapsed = difftime( time(NULL), tm_start);
            cpu_time_used     = ( (double) ( clock() - clk_start) ) / CLOCKS_PER_SEC;
            sprintf(msg, "Done: the global time-reversal symmetry, real time elasped %.0f s, cpu time (thread%3d) %.2f\n", real_time_elapsed, mpi_rank+1, cpu_time_used);
            print_msg(msg);
        }
    }
    // check the difference of symmed_hr and the original one
    if(mpi_rank==0 && flag_restart == 0 && flag_expandrvec == 0){ // check difference
        for(i=0;i<ham_final.nrpt*ham_final.norb*ham_final.norb;i++){
            if(cabs(ham_final.ham[i]  -  ham_in.ham[i]) > ham_tolerance){
                fstdout=fopen("wannsymm.out", "a");
                fprintf(fstdout,"WARNING: mismatch at rvec_orb:%5d%5d%5d%5d%5d\tdis=%16.9lf\n",(int)ham_final.rvec[i/(ham_in.norb*ham_in.norb)].x,
                                                          (int)ham_final.rvec[i/(ham_in.norb*ham_in.norb)].y,
                                                          (int)ham_final.rvec[i/(ham_in.norb*ham_in.norb)].z,
                                                          (i%(ham_in.norb*ham_in.norb))%ham_in.norb+1,
                                                          (i%(ham_in.norb*ham_in.norb))/ham_in.norb+1, cabs(ham_final.ham[i]  -  ham_in.ham[i]));
                fclose(fstdout);
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    //3.4 write the final hamitonian
    if( mpi_rank == 0 && flag_restart == 0){
        sprintf(ham_out_seed, "%s_symmed", seed);
        write_ham(&ham_final, ham_out_seed);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // prepare the rotation matrice in reciprocal-space
    for(isymm=0;isymm<nsymm;isymm++){
        //-----// R_KD = A·A^T·R·(A·A^T)^-1 = (R^T)^-1
        matrix3x3_transpose(rots_kd[isymm], rotations[isymm]);
        matrix3x3_inverse(rots_kd[isymm], rots_kd[isymm]);
    }

    int en_print_prec = 6;  // precision of energy for printing
    int bnd_print_len = 3;
    en_print_prec = (int) fabs(log10(degenerate_tolerance)) + 1;
    bnd_print_len = (int) fabs(log10(norb)) + 1;

    char * fn_out_cha="bnd_sym_characters";
    char * fn_out_eig="bnd_sym_eig";
    char * fn_out_bnd_sym="bands_symmed.dat";
    char * fn_out_bnd_ori="bands_ori.dat";
    FILE * fbands;

    vector firstk_in_kpath, lastk_in_kpath, kpt_inc;
    int ik_in_kpath;
    double ** ebands;
    vector * kpt_in_kpath;

    //3.5 calculate band structures along high-symmetry lines
    if( (flag_bands > 0 || flag_chaeig_in_kpath == 1)  && mpi_rank == 0){
        // derive k-poionts in k-paths
        kpt_in_kpath = (vector *) malloc(sizeof(vector)  *nkpath*nk_per_kpath);
        nerr = 1;
        for(ikpath=0; ikpath<nkpath && nerr != -1; ikpath++){
            firstk_in_kpath = vec_llist_pop(&kpaths, &nerr);
            lastk_in_kpath  = vec_llist_pop(&kpaths, &nerr);
            kpt_inc = vector_scale(1.0/(double)(nk_per_kpath-1), vector_sub(lastk_in_kpath, firstk_in_kpath));
            for(ik_in_kpath=0; ik_in_kpath<nk_per_kpath; ik_in_kpath++){
                ikpt = ikpath*nk_per_kpath + ik_in_kpath;
                kpt_in_kpath[ikpt] = vector_add(firstk_in_kpath, vector_scale((double) ik_in_kpath, kpt_inc));
                if(flag_chaeig_in_kpath == 1){
                    vec_llist_add(&kpts, kpt_in_kpath[ikpt]);
                    nkpt++;
                }
            }
        }

        // calculate bands
        if(flag_bands > 0){
            ebands       = (double **)malloc(sizeof(double *)*nkpath*nk_per_kpath);
            for(ikpt=0; ikpt<nkpath*nk_per_kpath; ikpt++){
                ebands[ikpt] = (double *)malloc(sizeof(double)*norb);
            }
            if( flag_bands%2 == 1){
                for(ikpt=0; ikpt<nkpath*nk_per_kpath; ikpt++){
                    bnd_eig_hk(ebands[ikpt], NULL, &ham_final, kpt_in_kpath[ikpt]);
                }
                fbands = fopen(fn_out_bnd_sym, "w");
                fprintf(fbands, "# Bands derived from symmetrized Hamiltonian ( %s_symmed_hr.dat )\n", seed);
                bnd_write_bands(fbands, norb, nkpath, nk_per_kpath, lattice, kpt_in_kpath, ebands);
                fclose(fbands);
            }
            if( flag_bands/2 == 1){
                for(ikpt=0; ikpt<nkpath*nk_per_kpath; ikpt++){
                    bnd_eig_hk(ebands[ikpt], NULL, &ham_in, kpt_in_kpath[ikpt]);
                }
                fbands = fopen(fn_out_bnd_ori, "w");
                fprintf(fbands, "# Bands derived from original Hamiltonian ( %s_hr.dat ) \n", seed);
                bnd_write_bands(fbands, norb, nkpath, nk_per_kpath, lattice, kpt_in_kpath, ebands);
                fclose(fbands);
            }
            for(ikpt=0; ikpt<nkpath*nk_per_kpath; ikpt++){
                free(ebands[ikpt]);
            }
            free(ebands);
        }
        free(kpt_in_kpath);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    //3.6 calculate characters and eigenvalues
    if( flag_chaeig == 1 && mpi_rank == 0){
        sprintf(msg,"Calculating bands' characters...\n");
        print_msg(msg);
        sprintf(msg, "related characters of symmetry operations will be written to the file \"%s\"\n", fn_out_cha);
        print_msg(msg);
        sprintf(msg, "related eigenvalues of symmetry operations will be written to the file \"%s\"\n\n", fn_out_eig);
        print_msg(msg);

        smallgroup sgrp;
        int sgrpi;
        ndegen = (int *)malloc(sizeof(int)*norb);
        eig_hk = (double *)malloc(sizeof(double)*norb);
        sgrp.element = (int *) malloc(sizeof(int)*nsymm*2); // nsymm*2 for double group

        fbndsymcha = fopen(fn_out_cha, "w");
        if(nkpt > 1)
            fprintf(fbndsymcha,  "Characters of symmetry operations    nkpt= %d  nbnd= %d\n\n", nkpt, norb);
        fclose(fbndsymcha);
        fbndsymeig = fopen(fn_out_eig, "w");
        if(nkpt > 1)
            fprintf(fbndsymeig, "Eigenvalues of symmetry operations    nkpt= %d  nbnd= %d\n\n", nkpt, norb);
        fclose(fbndsymeig);
        nerr = 1;
        for(ikpt=0; ikpt<nkpt && nerr!=-1 ; ikpt++){
            kpt = vec_llist_pop(&kpts, &nerr);

            bnd_eigcha(eig_hk, ndegen, &sym_chas, &sym_eigs, &sgrp,
                       lattice, rotations, translations, TR, rots_kd, nsymm,
                       orb_info, flag_soc, 
                       &ham_final, kpt, degenerate_tolerance);

        
            // write characters
            bnd_write_characters(fn_out_cha, ikpt, kpt, &sgrp, lattice, rotations, translations, TR, eig_hk, ndegen, sym_chas, norb, flag_soc, en_print_prec, bnd_print_len);

            // write eigenvalues
            bnd_write_eigenvalues(fn_out_eig, ikpt, kpt, &sgrp, lattice, rotations, translations, TR, eig_hk, ndegen, sym_eigs, norb, flag_soc, en_print_prec, bnd_print_len);

            //free the memory allocated this step
            for(sgrpi=0; sgrpi<sgrp.order; sgrpi++){
                free(sym_chas[sgrpi]);
                free(sym_eigs[sgrpi]);
            }
            free(sym_chas);
            free(sym_eigs);
        }
        free(ndegen);
        free(eig_hk);
        free(sgrp.element);
    }

    MPI_Barrier(MPI_COMM_WORLD);


    if(mpi_rank==0){
        finalize_wanndata(ham_final);
        finalize_wanndata(ham_in);
    }
    free(orb_info);
    MPI_Barrier(MPI_COMM_WORLD);

    if( mpi_rank == 0){
        real_time_elapsed = difftime( time(NULL), tm_start);
        cpu_time_used     = ( (double) ( clock() - clk_start) ) / CLOCKS_PER_SEC;
        sprintf(msg, "\nAll Done.  Real time elasped %.0f s, cpu time (thread%3d) %.2f\n", real_time_elapsed, mpi_rank+1, cpu_time_used);
        print_msg(msg);
    }
    MPI_Finalize();
    return 0;
}

