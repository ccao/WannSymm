#include "rotate_ham.h"

//#define __DEBUG
//#define __DEBUG_transmit
//#define __DEBUG_rot_orb
//#define __DEBUG_partial_d

void rotate_ham(wanndata * hout, wanndata * hin, double lattice[3][3], double rotation[3][3], double translation[3], wannorb * orb_info, int flag_soc, int flag_local_axis, int index_of_sym)
{
    //--This function generate a new Hamiltonian by performing the symmetry operation defined by the input 
    //  rotation and translation.
    //
    //  The Method is: 
    //      firstly, for every location loc_in  = rvec_in + orbital_site_in, find the rotated loccation:
    //      loc_out = rvec_out + orb_site_out
    //      then for every component of hout->ham[loc_out] corresponding to the loc_out, find the corresponding
    //      linear combination of components in hin->ham
    //      which is 
    //          hout->ham[loc_out] = D * hin->ham[loc_in] D^{dagger}
    //      where D is the rotation matrix.
    //
    //  output:
    //      hout        : generated Hamiltonian
    //  input:
    //      hin         : original Hamiltonian
    //      lattice     : crystal lattice
    //      rotation    : 3x3 matrix defining the rotation part of symmetry operation
    //      translation : 1x3 vector defining the translational part of the symmetry operation
    //      orb_info    : contains sites, l, mr, ms, r, axis of every orbital
    //      flag_soc    : consider SOC or not.
    //
    
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    FILE * fstdout;
    int i, j, k;
    int ii, jj, kk;
    int iorb, jorb;
    int io,jo;
    int irpt, jrpt;
    int irpt_in;
    int isite, jsite;
    int il, jl;
    int imr, jmr;
    int ims, jms;
    vector site_ws_in1, site_ws_out1; //site converted to the 0 0 0 w-s unit cell, for which, 
    vector site_ws_in2, site_ws_out2; // all the conponents of the vector should between -0.5 and 0.5
    vector site_out1, site_in1;
    vector site_out2, site_in2;
    vector site_out, site_in;
    vector rvec_in1, rvec_out1;
    vector rvec_in2, rvec_out2;
    vector rvec_in, rvec_out, rvec_in_roted;
    vector loc_in1, loc_out1;
    vector loc_in2, loc_out2;
    vector loc_in, loc_out;
    
    //double rotation[3][3];
    double rot_combined[3][3];
    double inv_rotation[3][3];
    double inv_translation[3];
    int norb;
    double rot_axis[3];
    double rot_angle;

    int N;
    dcomplex * orb_rot[MAX_L+1];
    dcomplex * s_rot;
    int inv_flag;
    int r, r1, r2, l, l1, l2, mr_i, mr_j, mr1, mr2, ms_i, ms_j, ms1, ms2;
    dcomplex Ham_component;

    int ii_out, ii_in;
    int irpt_out, iorb_out, jorb_out;
    int iorb_in;
    vector site_in1_old, site_in2_old;
    wanndata htmp;

    char msg[MAXLEN];

    norb = hin->norb;

    inverse_symm(rotation, inv_rotation, translation, inv_translation);

    get_axis_angle_of_rotation(rot_axis, &rot_angle, &inv_flag, rotation, lattice);

    s_rot = (dcomplex *)malloc(sizeof(dcomplex)*2*2);                       //get s_rot with input (axis,angle,inv)
    rotate_spinor(s_rot, rot_axis, rot_angle, inv_flag);


    for (l=0;l<=3;l++){                                                     //get orb_rot with input (l,axis,angle,inv)
        N = 2*l + 1;
        orb_rot[l] = (dcomplex *)malloc(sizeof(dcomplex)*N*N);
        rotate_cubic( orb_rot[l], l, rot_axis, rot_angle, inv_flag); 
    }
#ifdef __DEBUG_rot_orb
    // print orbital rotation matrix
    char fnrot[MAXLEN];
    FILE * frot;
    sprintf(fnrot, "orb_rot_%d", mpi_rank+1);
    frot=fopen( fnrot, "w");
    fprintf(frot, "orb_rot:\n");
    for(l=0;l<=3;l++){
        fprintf(frot, "l=%d\n", l);
        fprintf(frot, "[\n");
        for(io=0;io<2*l+1;io++){
            for(jo=0;jo<2*l+1;jo++){
                fprintf(frot, "%12.7lf+%12.7lf*i,", creal(orb_rot[l][io*(2*l+1)+jo]), cimag(orb_rot[l][io*(2*l+1)+jo]) );
            }
            fprintf(frot, ";\n");
        }
        fprintf(frot, "]\n");
    }

    if( flag_soc == 1){
        fprintf(frot, "\n\nspin_rot:\n");
        for(io=0;io<2;io++){
            for(jo=0;jo<2;jo++){
                fprintf(frot, "%12.7lf+%12.7lf*i,", creal(s_rot[io*2+jo]), cimag(s_rot[io*2+jo]) );
            }
            fprintf(frot, "\n");
        }
    }
    fclose(frot);
#endif


    //create a table to speed up finding rvec of (inv) symmetry operated orbital
    vector * site_symmed, * rvec_sup_symmed;
    // site_symmed    : site of orbital resulting from (inv) symmetry operation on each orbital.  
    // rvec_sup_symmed: extra R vector gained from (inv) symmetry operation on each orbital.
    vector * site_invsed, * rvec_sup_invsed; 
    site_symmed = (vector *) malloc(sizeof(vector)*norb);
    site_invsed = (vector *) malloc(sizeof(vector)*norb);
    rvec_sup_symmed = (vector *) malloc(sizeof(vector)*norb);
    rvec_sup_invsed = (vector *) malloc(sizeof(vector)*norb);
    for(iorb=0;iorb<norb;iorb++){
        if( iorb > 0 && equal((orb_info+iorb)->site, (orb_info+iorb-1)->site)){
            *(rvec_sup_symmed + iorb) = *(rvec_sup_symmed + iorb -1);
            *(rvec_sup_invsed + iorb) = *(rvec_sup_invsed + iorb -1);
            *(site_symmed + iorb) = *(site_symmed + iorb -1);
            *(site_invsed + iorb) = *(site_invsed + iorb -1);
            continue;   
        }
        loc_in  = (orb_info+iorb)->site;

        // Get rvec_sup_symmed and site_symmed on symmetry operation
        loc_out = vector_add(vector_rotate(loc_in, rotation), array2vector(translation));
        getrvec_and_site(rvec_sup_symmed + iorb, site_symmed + iorb, loc_out, orb_info, norb, lattice);

        // Get rvec_sup_invsed and site_invsed on inverse symmetry operation
        loc_out = vector_add(vector_rotate(loc_in, inv_rotation), array2vector(inv_translation));
        getrvec_and_site(rvec_sup_invsed + iorb, site_invsed + iorb, loc_out, orb_info, norb, lattice);
    }


    // Rotate all site = {Rvec + tau} to site_rot, if site_rot located in a Rvec not listed in rvec_list, add it.
    vec_llist * rvecs;
    int nrvec=0;

    vec_llist_init(&rvecs);
    for( irpt=0; irpt < hin->nrpt; irpt++){
        vec_llist_add(&rvecs, hin->rvec[irpt]);           // use add (to end) for consistency with v1.0.0-rc
        //vec_llist_add_inorder(&rvecs, hin->rvec[irpt]);   // use add_inorder for more convenient output
        nrvec++;
    }
    for( irpt=0; irpt < hin->nrpt; irpt++){
        rvec_in = hin->rvec[irpt];
        rvec_in_roted = vector_rotate(rvec_in, rotation);
        for(jorb=0; jorb<norb; jorb++){
            if( jorb > 0 && equal((orb_info+jorb)->site, (orb_info+jorb-1)->site)) continue;
            rvec_out2 = vector_add(rvec_sup_symmed[jorb], rvec_in_roted);
            for(iorb=0; iorb<norb; iorb++){
                if( iorb > 0 && equal((orb_info+iorb)->site, (orb_info+iorb-1)->site)) continue;
                rvec_out1 = rvec_sup_symmed[iorb];

                rvec_out = vector_sub(rvec_out2, rvec_out1);
                irpt_out = vec_llist_find(&rvecs, rvec_out);
                if(irpt_out == -1){ 
                    vec_llist_add(&rvecs, rvec_out);           // use add (to end) for consistency with v1.0.0-rc
                    //vec_llist_add_inorder(&rvecs, rvec_out);   // use add_inorder for more convenient output
                    nrvec++;
                }
            }
        } 
    }

    hout->norb=hin->norb;
    hout->nrpt = nrvec;
    init_wanndata(hout);
    
    int nerr=1;
    for( irpt=0; irpt < nrvec; irpt++){
        hout->rvec[irpt] = vec_llist_pop(&rvecs, &nerr);
        if( nerr == -1){
            sprintf(msg, "Error occured in finding rotated R vectors.");
            print_error(msg);
            exit(1);
        }
    }
    // End Rotate all site
    
    vector rvec_out_invsed;
    char foutname[MAXLEN];
    FILE * fout;

    sprintf(foutname,".progress-of-thread%d", mpi_rank+1);
    remove(foutname);
    //fout=fopen(foutname, "w");
    for( irpt_out=0; irpt_out < nrvec; irpt_out++){
        fout=fopen(foutname, "a");
        fprintf(fout, "Symm No. %d, progress %5.2lf%% (%d/%d)\n", index_of_sym+1, (double)(irpt_out+1)/(double)nrvec*100, irpt_out+1, nrvec);
        fclose(fout);
        rvec_out = hout->rvec[irpt_out];
        rvec_out_invsed = vector_rotate(rvec_out, inv_rotation);
        for(jorb_out=0; jorb_out<norb; jorb_out++){
            site_out2 = (orb_info + jorb_out)->site;
            site_in2  = site_invsed[jorb_out];
            rvec_in2  = vector_add(rvec_sup_invsed[jorb_out], rvec_out_invsed);
            r2 = (orb_info+jorb_out)->r;
            l2 = (orb_info+jorb_out)->l;
            mr_j = (orb_info+jorb_out)->mr;
            ms_j = (orb_info+jorb_out)->ms;
            for(iorb_out=0; iorb_out<norb; iorb_out++){
                site_out1 = (orb_info + iorb_out)->site;
                site_in1  = site_invsed[iorb_out];
                rvec_in1  = rvec_sup_invsed[iorb_out];
                r1 = (orb_info+iorb_out)->r;
                l1 = (orb_info+iorb_out)->l; 
                mr_i = (orb_info+iorb_out)->mr;
                ms_i = (orb_info+iorb_out)->ms;

                rvec_in = vector_sub(rvec_in2, rvec_in1);
                irpt_in = find_vector(rvec_in, hin->rvec, hin->nrpt);
                if(irpt_in == -1){
                    continue;
                }
                ii_out = irpt_out*norb*norb + jorb_out*norb + iorb_out;
                hout->hamflag[ii_out] = 1;
                if (flag_soc == 1){
                    for(mr1 = 1; mr1 <= 2*l1+1; mr1++){
                        for(mr2 = 1; mr2 <= 2*l2+1; mr2++){
                            for(ms1 = 0; ms1 <=1; ms1++){
                                for(ms2 = 0; ms2 <= 1; ms2++){
                                    ii_in=find_index_of_ham(hin, orb_info, norb, irpt_in,
                                                            site_in1, r1, l1, mr1, ms1,
                                                            site_in2, r2, l2, mr2, ms2);
                                    if(ii_in < 0) continue;
                                    if(flag_local_axis > 0){
                                        //matrix3x3_dot(rot_combined, rotation, );
                                        iorb_in = ii_in % norb;
                                        combine_rot_with_local_axis(rot_combined, rotation, orb_info, iorb_in, iorb_out);
                                        get_axis_angle_of_rotation(rot_axis, &rot_angle, &inv_flag, rot_combined, lattice);
                                        s_rot = (dcomplex *)malloc(sizeof(dcomplex)*2*2);
                                        //get s_rot with input (axis,angle,inv)
                                        rotate_spinor(s_rot, rot_axis, rot_angle, inv_flag);

                                        for (l=0;l<=3;l++){
                                            //get orb_rot with input (l,axis,angle,inv)
                                            N = 2*l + 1;
                                            orb_rot[l] = (dcomplex *)malloc(sizeof(dcomplex)*N*N);
                                            rotate_cubic( orb_rot[l], l, rot_axis, rot_angle, inv_flag); 
                                        }
                                    }
                                    // roted_H(l1,l2) = D(l1) · S · H(l1,l2) · S.conj.transe · D(l2).conj.transe
                                    hout->ham[ii_out] += orb_rot[l1][(2*l1+1)*(mr_i-1)+mr1-1] *
                                                         s_rot[2*ms_i + ms1] * 
                                                         (hin->ham[ii_in] / hin->weight[irpt_in] )*
                                                         conj(s_rot[2*ms_j + ms2])       * 
                                                         conj(orb_rot[l2][(2*l2+1)*(mr_j-1)+mr2-1]);
                                }
                            }
                        }
                    }
                } else if (flag_soc == 0){
                    ms1=ms2=0;
                    for(mr1 = 1; mr1 <= 2*l1+1; mr1++){
                        for(mr2 = 1; mr2 <= 2*l2+1; mr2++){
                            ii_in=find_index_of_ham(hin, orb_info, norb, irpt_in, 
                                                    site_in1, r1, l1, mr1, ms1,
                                                    site_in2, r2, l2, mr2, ms2);
                            if(ii_in < 0) continue;
                            // roted_H(l1,l2) = D(l1) · H(l1,l2) · D(l2).conj.transe
                            hout->ham[ii_out] += orb_rot[l1][(2*l1+1)*(mr_i-1)+mr1-1] *
                                                 (hin->ham[ii_in] / hin->weight[irpt_in]) * 
                                                 conj(orb_rot[l2][(2*l2+1)*(mr_j-1)+mr2-1]);
                        }
                    }
                }
                // may be only in QE convention(or nonsoc for vasp), the block multiply can be used.
                //ii_in  = get_in_block( lattice, rot, trans, orb_info, ii_out);
                // ii_in should be start of the block.
                // use inv_symm to find the block_in
                // roted_H(l1,l2) = D(l1) · S · H(l1,l2) · S.conj.transe · D(l2).conj.transe
                //--S and D can be merge only for QE convention--//merge_orb_s_rot(orb_s_rot, orb_rot, s_rot);
                //cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, N, N, N, &one, hin->ham+ii_in, N, orb_s_rot, N, &zero, tmp, N);
                //cblas_zgemm(CblasRowMajor, CblasTrans, CblasNoTrans, N, N, N, &one, orb_s_rot, N, tmp, N, &zero, hout->ham+ii_out, N);
            }
        }
    }

    free(s_rot);
    for (l=0;l<=3;l++){
        free(orb_rot[l]);
    }
}

void getrvec_and_site(vector * p2rvec, vector * p2site, vector loc, wannorb * orb_info, int norb, double lattice[3][3]){
    //input:  location of a atom : loc
    //output: rvec and site of the atom : *p2rvec *p2site
    int    i,j,k,ii,jj,kk;
    vector dis, dis_int, dis_rem, dis_rem_cartesian;
    init_vector(&dis, 0,0,0);
    init_vector(&dis_int, 0,0,0);
    init_vector(&dis_rem, 0,0,0);
    init_vector(p2rvec, 0,0,0);
    double eps=1E-3;
    for(i=0;i<norb;i++){
        dis.x = loc.x - (orb_info+i)->site.x;
        dis.y = loc.y - (orb_info+i)->site.y;
        dis.z = loc.z - (orb_info+i)->site.z;
        dis_int.x = round(dis.x);
        dis_int.y = round(dis.y);
        dis_int.z = round(dis.z);
        dis_rem.x = dis.x - dis_int.x;
        dis_rem.y = dis.y - dis_int.y;
        dis_rem.z = dis.z - dis_int.z;
        dis_rem_cartesian = vector_rotate(dis_rem, lattice);
        if( sqrt(dot_product( dis_rem_cartesian, dis_rem_cartesian )) < eps){
            *p2site = (orb_info+i)->site;
            *p2rvec = dis_int;
            break;
        }
    }
    if( sqrt(dot_product( dis_rem_cartesian, dis_rem_cartesian )) >= eps){
        fprintf(stderr,"ERROR: can not find rotated rvec and site.\n");
        fprintf(stderr,"loc=(%10.5lf,%10.5lf,%10.5lf) dis=(%10.5lf,%10.5lf,%10.5lf) dis_int=(%5d%5d%5d) dis_rem=(%10.5lf,%10.5lf,%10.5lf)\n",
                loc.x, loc.y, loc.z, dis.x, dis.y, dis.z, (int)dis_int.x, (int)dis_int.y, (int)dis_int.z, dis_rem.x, dis_rem.y, dis_rem.z);
        exit(1);
    }
}

void get_axis_angle_of_rotation(double axis[3], double * angle, int * inv, double rin[3][3], double lattice[3][3]){
    double determinant=0;
    double rot[3][3];
    double trace_rot=0;
    vector axis_v;
    int i,j;
    int ii,jj;
    int tr;             //type of rotation symmetry
    double norm;

    double tmp3x3[3][3];
    double rotation_cartesian[3][3];
    double c_m_lattice[3][3];
    double inv_lattice[3][3];

    for(i=0;i<3;i++)
        for(j=0;j<3;j++)
            c_m_lattice[i][j]=lattice[j][i];                 //row major to cloumn major

    matrix3x3_inverse(inv_lattice, c_m_lattice);             //inv_lattice = inverse matrix of lattice
    matrix3x3_dot(tmp3x3, c_m_lattice, rin);                 //tmp3x3 = lattice' * rotation
    matrix3x3_dot(rotation_cartesian, tmp3x3, inv_lattice);  //rotation_cartesian= lattice*rotation * inv_lattice

    for(i = 0; i < 3; i++){
        determinant += (rotation_cartesian[0][i] * (rotation_cartesian[1][(i+1)%3] * rotation_cartesian[2][(i+2)%3] - rotation_cartesian[1][(i+2)%3] * rotation_cartesian[2][(i+1)%3]));
    }

    if(determinant<0)
            *inv=1;
    else
            *inv=0;

//trans the rotation to proper rotation
    for(i=0;i<3;i++){
        for(j=0;j<3;j++){
            //rot[i][j] = rin[i][j] * sign(determinant);
            rot[i][j] = rotation_cartesian[i][j] * sign(determinant);
        }
    }
    

    tr=type_of_rotation(rot);

    if(tr > 6 || tr < 0){// no use
            printf("! ERROR, symmetry type not compatible to find axis and angle\n");
            //exit(0);
    }
    
    if(tr == 1 || tr==2){
        // Identity rotation, Inversion
        axis[0]=0;
        axis[1]=0;
        axis[2]=1;
        *angle = 0;
    }
    else if(tr==4 || tr==5){
        // 180 proper rotation
        // First the case where the axis is parallel to a coordinate axis
        if ( tr==5){    // mirror
             for(i=0;i<3;i++)
                 for(j=0;j<3;j++){
                     rot[i][j] *= -1.0;
                 }
        }
        *angle = PI;    //tr==4 : a 180 rotation
        for (i=0;i<3;i++){
            axis[i]=0;
            if(fabs(rot[i][i]-1.0)<eps7 ) axis[i] =1.0;
        }
        norm = sqrt(pow(axis[0], 2) + pow(axis[1],2) + pow(axis[2],2));

        if (fabs(norm)<eps7) {
            // then the general case
            for(i=0;i<3;i++){
                    axis[i]=sqrt(fabs(rot[i][i] + 1.0) /2.0);
            }
            for(i=0;i<3;i++)
                for(j=i+1;j<3;j++)
                    if(fabs(axis[i]*axis[j]) > eps7)   axis[i]=0.5 * rot[i][j]/axis[j];
        }
    }
    else{   //tr==3 or tr==6 It is not a 180 rotation: compute the rotation axis
        axis[0] = rot[2][1] - rot[1][2];
        axis[1] = rot[0][2] - rot[2][0];
        axis[2] = rot[1][0] - rot[0][1];
        //normalize
        norm = sqrt(pow(axis[0], 2) + pow(axis[1],2) + pow(axis[2],2));
        for (i=0;i<3;i++) 
            axis[i] /= norm;

        for(i=0;i<3;i++) 
            trace_rot += rot[i][i];
        *angle = acos((trace_rot-1.0)/2.0);

        if( fabs(axis[1]*axis[0]*(1-cos(*angle)) + axis[2]*sin(*angle) - rot[1][0]) > 1e-5 )
            *angle *= -1;
        if( fabs(axis[1]*axis[0]*(1-cos(*angle)) + axis[2]*sin(*angle) - rot[1][0]) > 1e-3 )
            *angle *= -1;
        if( fabs(axis[1]*axis[0]*(1-cos(*angle)) + axis[2]*sin(*angle) - rot[1][0]) > 0.2 )
            *angle *= -1;

        //if(sign(rot[1][0]) != sign(*angle))
        if( fabs(axis[1]*axis[0]*(1-cos(*angle)) + axis[2]*sin(*angle) - rot[1][0]) > 0.1 ){
            fprintf(stderr, "!ERROR: Can not find corresponding axis & angle \n");
            fprintf(stderr, "rot_direct:        rot_Cartesian:\n");
            for(ii=0;ii<3;ii++){
                for(jj=0;jj<3;jj++)
                    fprintf(stderr, "%9.5lf ",rin[ii][jj]);
                fprintf(stderr, "  ||   ");
                for(jj=0;jj<3;jj++)
                    fprintf(stderr, "%9.5lf ",rot[ii][jj]);
                fprintf(stderr,"\n");
            }
            fprintf(stderr,"rot_Cartesian[1][0] =%15.9lf\naxis[1]*axis[0]*(1-cos(*angle)) + axis[2]*sin(*angle)=%15.9lf\n", 
                            rot[1][0], axis[1]*axis[0]*(1-cos(*angle)) + axis[2]*sin(*angle));
            fprintf(stderr,"Generally, rot_Cartesian[1][0] should be equal to axis[1]*axis[0]*(1-cos(*angle)) + axis[2]*sin(*angle)\n");
            fprintf(stderr,"angle=%6.2lf, axis=%7.3lf%7.3lf%7.3lf\n",(*angle)/PI*180,axis[0],axis[1],axis[2]);
            exit(0);
        }

    }
    
    //The direction of the axis is arbitrarily chosen, with positive z. In the
    //xy plane with positive x, and along y with positive y.
    if ( axis[2] < -eps7 ){
        for(i=0;i<3;i++) axis[i] = -axis[i]; 
        *angle *= -1;
    }
    else if ( fabs(axis[2]) < eps7 && axis[0] < -eps7 ){
        for(i=0;i<3;i++) axis[i] = -axis[i];
        *angle *= -1;
    }
    else if ( fabs(axis[2]) < eps7 && fabs(axis[0]) < eps7 && axis[1] < -eps7){
        for(i=0;i<3;i++) axis[i] = -axis[i];
        *angle *= -1;
    }

    norm = sqrt(pow(axis[0], 2) + pow(axis[1],2) + pow(axis[2],2));
    //norm = sqrt(axis[0]**2 + axis[1]**2 + axis[2]**2);
    if ( norm < eps7) {
            fprintf(stderr,"! ERROR \"get_axis_angle_of_rotation\" problem with the matrix\n");
            exit(0);
    }
    for (i=0;i<3;i++) axis[i] /= norm;
    
    if( *angle < -PI + 1E-3){
        *angle += 2*PI; // rot_angle range (-PI, PI]
    }

}


int type_of_rotation(double rot[3][3]){
//     1   Identity
//     2   Inversion
//     3   Proper rotation of an angle <> 180 degrees
//     4   Proper rotation of 180 degrees
//     5   Mirror symmetry
//     6   Improper rotation
    double det, det1;
    int i;
    //
    //  Check for Identity
    //
    if ( fabs(rot[0][0]-1.0) < eps7 &&
         fabs(rot[1][1]-1.0) < eps7 &&
         fabs(rot[2][2]-1.0) < eps7 &&
         fabs(rot[0][1]) < eps7 && fabs(rot[1][2]) < eps7 && fabs(rot[0][2]) < eps7 &&
         fabs(rot[1][0]) < eps7 && fabs(rot[2][1]) < eps7 && fabs(rot[2][0]) < eps7 )
    {   return 1;}

    //
    //   Check for inversion
    //
    if ( fabs(rot[0][0]+1.0) < eps7 &&
         fabs(rot[1][1]+1.0) < eps7 &&
         fabs(rot[2][2]+1.0) < eps7 &&
         fabs(rot[0][1]) < eps7 && fabs(rot[1][2]) < eps7 && fabs(rot[0][2]) < eps7 &&
         fabs(rot[1][0]) < eps7 && fabs(rot[2][1]) < eps7 && fabs(rot[2][0]) < eps7 )
    {   return 2;}

    //
    //  Compute the determinant
    //
    det=0;
    for(i=0;i<3;i++){
        det += (rot[0][i] * (rot[1][(i+1)%3] * rot[2][(i+2)%3] - rot[1][(i+2)%3] * rot[2][(i+1)%3]));
    }
    
    //
    //  Determinant equals to 1: proper rotation
    //
    if(fabs(det-1.0)<eps7){
        //check if an eigenvalue is equal to -1.d0 (180 rotation)
        det1 = (rot[0][0]+1.0)*((rot[1][1]+1.0)*(rot[2][2]+1.0)-rot[2][1]*rot[1][2])-   
                rot[0][1]*      (rot[1][0]*     (rot[2][2]+1.0)-rot[2][0]*rot[1][2])+   
                rot[0][2]*      (rot[1][0]*rot[2][1]       -rot[2][0]*(rot[1][1]+1.0));

        if(fabs(det1) < eps7 )
            return 4;   //180 proper rotation
        else
            return 3;   //proper rotation <> 180
        //endif
    }

    //
    // Determinant equal to -1: mirror symmetry or improper rotation
    //
    if(fabs(det+1.0) < eps7) {
        // check if an eigenvalue is equal to 1.d0 (mirror symmetry)
        det1 = (rot[0][0]-1.0)*((rot[1][1]-1.0)*(rot[2][2]-1.0)-rot[2][1]*rot[1][2])-   
                rot[0][1]*      (rot[1][0]*     (rot[2][2]-1.0)-rot[2][0]*rot[1][2])+   
                rot[0][2]*      (rot[1][0]*rot[2][1]       -rot[2][0]*(rot[1][1]-1.0));

        if(fabs(det1)<eps7)
            return 5;   //mirror symmetry
        else
            return 6;   //improper rotation
    }

    // error
    return -1;
}

void inverse_symm(double rin[3][3], double rout[3][3], double tin[3], double tout[3]){
    double determinant=0;
    int i,j;

    matrix3x3_inverse(rout, rin);

    vector tin_v, tout_v;
    init_vector(&tin_v, tin[0],tin[1],tin[2]);
    tout_v = vector_scale(-1.0, vector_rotate(tin_v, rout));
    tout[0]=tout_v.x;
    tout[1]=tout_v.y;
    tout[2]=tout_v.z;

}

double sign(double in){
    if(in>0)
        return 1.0;
    else if(in<0)
        return -1.0;
    else
        return 0.0;
}


void trsymm_ham(wanndata * hout, wanndata * hin, wannorb * orb_info, int flag_soc){
    int irpt,iorb,jorb,ii;
    int num_orb;
    int r, r1, r2, l, l1, l2, mr_i, mr_j, mr1, mr2, ms_i, ms_j, ms1, ms2;
    vector rvec_in1, rvec_out1;
    vector rvec_in2, rvec_out2;
    vector rvec_in,  rvec_out;
    vector site_out1, site_in1;
    vector site_out2, site_in2;
    //dcomplex sigma_y[4] = {0,-I,I,0};
    //dcomplex tr_factor[4]  = {0,-1,1,0};
    dcomplex tr_factor[4]  = {0, 1,-1,0};
    
    num_orb = hin->norb;

    /* The symmetry operation should not change Rvec and weights */
    memcpy(hout->rvec,   hin->rvec, sizeof(vector)*hin->nrpt);
    memcpy(hout->weight, hin->weight, sizeof(int)*hin->nrpt);

    if(flag_soc == 0){
        for(irpt=0; irpt< hin->nrpt; irpt++){
            for(jorb=0;jorb < num_orb; jorb++){
                for(iorb=0;iorb < num_orb; iorb++){
                    ii= irpt*num_orb*num_orb + jorb*num_orb + iorb;
                    hout->ham[ii] = conj(hin->ham[ii]);
                }
            }
        }
    }
    else if (flag_soc==1){
        for(irpt=0; irpt< hin->nrpt; irpt++){
            rvec_out = rvec_out2  = hin->rvec[irpt];    // rvec added to the second oribital
            init_vector(&rvec_out1, 0, 0, 0);
            rvec_in = rvec_out;
            for(jorb=0;jorb < num_orb; jorb++){
                r2   = (orb_info+jorb)->r;
                l2   = (orb_info+jorb)->l;
                mr_j = (orb_info+jorb)->mr;
                mr2  = mr_j;
                ms_j = (orb_info+jorb)->ms;
                site_out2 = (orb_info+jorb)->site;
                site_in2=site_out2;
                for(iorb=0;iorb < num_orb; iorb++){
                    r1   = (orb_info+iorb)->r;
                    l1   = (orb_info+iorb)->l;
                    mr_i = (orb_info+iorb)->mr;
                    mr1  = mr_i;
                    ms_i = (orb_info+iorb)->ms;
                    site_out1 = (orb_info+iorb)->site;
                    site_in1=site_out1;
                    hout->ham[irpt*num_orb*num_orb + jorb*num_orb + iorb] = 0.0 ;
                    for(ms1 = 0; ms1 <=1; ms1++){
                        for(ms2 = 0; ms2 <= 1; ms2++){
                            ii=find_index_of_ham(hin, orb_info, num_orb, irpt,
                                                 site_in1, r1, l1, mr1, ms1,
                                                 site_in2, r2, l2, mr2, ms2);
                            hout->ham[irpt*num_orb*num_orb
                                      + jorb*num_orb + iorb] += tr_factor[ms_i*2+ms1] *
                                                                tr_factor[ms_j*2+ms2] *
                                                                conj(hin->ham[ii]);
                                                                //error// hin->ham[ii];
                        }
                    }
                }
            }
        }
    }
}

void combine_rot_with_local_axis(double rot_combined[3][3], double rotation[3][3], wannorb * orb_info, int io_in, int io_out){
    // rot_combined =  inv(Axes_out) * rot * Axes_in
    double axes_in[3][3];
    double axes_out[3][3];
    double mtmp[3][3];
    int ii;

    for(ii=0;ii<3;ii++){
        axes_in[ii][0] = orb_info[io_in].axis[ii].x;
        axes_in[ii][1] = orb_info[io_in].axis[ii].y;
        axes_in[ii][2] = orb_info[io_in].axis[ii].z;
        axes_out[ii][0] = orb_info[io_out].axis[ii].x;
        axes_out[ii][1] = orb_info[io_out].axis[ii].y;
        axes_out[ii][2] = orb_info[io_out].axis[ii].z;
    }
    matrix3x3_dot(rot_combined, rotation, axes_in);
    matrix3x3_inverse(mtmp, axes_out);
    matrix3x3_dot(rot_combined, mtmp, rot_combined);
}

