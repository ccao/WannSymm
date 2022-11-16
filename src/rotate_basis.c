#include "rotate_basis.h"

void get_sym_op_reciprocalspace(dcomplex * sym_op, double lattice[3][3], wannorb * orb_info, int norb, int isym_in_doublegp, double rotation[3][3],double translation[3], int TR, double rot_kd[3][3], vector kpt, int flag_soc, int flag_local_axis) {
   //
   //  trans from ccao's fortran code
   //
   //  This function defines the rotation in the Hilbert space
   //    <W'|R|W>
   //  For each Wannier orbital, we indeed have 4 labels:
   //    |W> = | Rv, tau, l, m>
   //    where Rv is the lattice vector
   //          tau is the atomic position
   //          l  is the angular momentum
   //          m  indexes the actual orbital (mostly cubic harmonics)
   //  Therefore
   //    <Rv', tau', l, m' | R | Rv, tau, l, m> is the matrix element
   //    lm part rotation matrix is exacly what we had for Ylm or cubic harmonics
   //    <Rv', tau' | R | Rv, tau> is a delta function
   //  For real-space rotation, above is just fine
   //    However, here we rotate |kv, tau, l, m> = \sum(Rv) exp(-i (kv*Rv)) | Rv, tau, l, m>
   //    under rotation, kv -> kv', Rv+tau -> Rv'+tau'
   //    notice that Rv' and R*Rv can differ by a lattice vector Rv", Rv'=R*Rv+R"
   //    thus under rotation, exp(-i (kv*Rv))=exp(-i (kv'*R*Rv))=exp(-i (kv'*Rv'))*exp(i(kv'*Rv"))
   //  Therefore
   //  \sum(Rv') exp(-i (kv'*Rv')) | Rv', tau', l, m>
   //     =\sum(Rv") exp(-i (kv'*Rv")) *exp(-i (kv'* R*Rv)) |Rv', tau', l, m>
   // 
   //   isym_in_doublegp : index of symm in double group, range [-nsymm, nsymm-1], nsymm is num of symm in single-valued group

    // if( fabs(translation[0]) > 1E-3 || fabs(translation[1]) > 1E-3 || fabs(translation[2]) > 1E-3){
    //     fprintf(stderr, "Warning, non symmorphic symmetry detected\n");
    // }

    int i,j,k;
    int ii,jj,kk;
    int io,jo;

    dcomplex * orb_rot[MAX_L+1];
    dcomplex * s_rot;
    double rot_axis[3];
    double rot_angle;
    int inv_flag;
    double inv_rot[3][3];

    vector kpt_roted;
    vector rvec_supp;       // tau' = rot*tau - rvec_supp
    vector tau1,tau2,tau2_roted, tau2_symed;
    vector translation_v;
    vector inv_trans;
    //vector translation_kd_v;
    int l, N;
    int mr1,mr2;
    int ms1,ms2;

    // vaiable used for deriving rotation matrix for local-axis
    double rot_combined[3][3];
    double la_rot_axis[3]; 
    double la_rot_angle;
    int la_inv_flag;

    init_vector(&translation_v, translation[0], translation[1], translation[2]);
    //init_vector(&translation_kd_v, translation_kd[0], translation_kd[1], translation_kd[2]);
    
    matrix3x3_inverse(inv_rot, rotation);
    inv_trans = vector_rotate(translation_v, inv_rot);
    inv_trans.x *= -1;
    inv_trans.y *= -1;
    inv_trans.z *= -1;

    get_axis_angle_of_rotation(rot_axis, &rot_angle, &inv_flag, rotation, lattice);
    if(flag_soc == 1) {
        if(isym_in_doublegp < 0) {
            rot_angle += 2*PI; // element introduced by double group
        }
    }
    //rot_angle *= -1; // proper rotation or inverse rotation
    //get s_rot with input (axis,angle,inv)
    s_rot = (dcomplex *)malloc(sizeof(dcomplex)*2*2);
    rotate_spinor(s_rot, rot_axis, rot_angle, inv_flag);
    if(flag_soc == 0){
        s_rot[0] = s_rot[3] = 1;
        s_rot[1] = s_rot[2] = 0;
    }

    //get orb_rot with input (l,axis,angle,inv)
    for (l=0;l<=3;l++){
        N = 2*l + 1;
        orb_rot[l] = (dcomplex *)malloc(sizeof(dcomplex)*N*N);
        rotate_cubic( orb_rot[l], l, rot_axis, rot_angle, inv_flag);
    }

    kpt_roted = vector_rotate(kpt, rot_kd);
    if(TR==1){
        kpt_roted = vector_scale(-1.0, kpt_roted);
    }
    for(io=0;io<norb;io++){
        for(jo=0;jo<norb;jo++){
            //if( ! kpt_equivalent(kpt_roted, kpt) ){
            //    sym_op[io*norb+jo] = 0;
            //    continue;
            //}
            if( (orb_info+io)->l != (orb_info+jo)->l ){
                sym_op[io*norb+jo] = 0;
                continue;
            }
            tau1 = (orb_info+io)->site;
            tau2 = (orb_info+jo)->site;
            tau2_roted = vector_rotate(tau2, rotation);
            tau2_roted = vector_add( tau2_roted, translation_v);
            //tau2_roted = vector_rotate(tau2, inv_rot);
            //tau2_roted = vector_add( tau2_roted, inv_trans);

            getrvec_and_site(&rvec_supp, &tau2_symed, tau2_roted, orb_info, norb, lattice);

            if( ! equale(tau1, tau2_symed, 1E-5)){
                sym_op[io*norb+jo] = 0;
                continue;
            }
            l   = (orb_info+jo)->l;
            mr1 = (orb_info+io)->mr;
            mr2 = (orb_info+jo)->mr;
            ms1 = (orb_info+io)->ms;
            ms2 = (orb_info+jo)->ms;
            if(flag_local_axis == 1){
                combine_rot_with_local_axis(rot_combined, rotation, lattice, orb_info, jo, io); // rotation from jo to io
                get_axis_angle_of_rotation(la_rot_axis, &la_rot_angle, &la_inv_flag, rot_combined, lattice);
                rotate_cubic( orb_rot[l], l, la_rot_axis, la_rot_angle, la_inv_flag);
            }
            //sym_op[io*norb+jo] =cexp(2*PI*cmplx_i * (dot_product(kpt_roted, vector_sub(rvec_supp,translation_v))))*
            //sym_op[io*norb+jo] =cexp(-2*PI*cmplx_i * (dot_product(kpt_roted, rvec_supp)))*
            //sym_op[io*norb+jo] =cexp(-2*PI*cmplx_i * (dot_product(kpt_roted, vector_add(rvec_supp,translation_v))))*
            sym_op[io*norb+jo] =cexp(-2*PI*cmplx_i * (dot_product(kpt_roted, vector_sub(rvec_supp,translation_v))))*
                                (orb_rot[l][(2*l+1)*(mr1-1) + mr2-1]) * 
                                (s_rot[2*ms1 + ms2]);
            //sym_op[io*norb+jo] = conj(sym_op[io*norb+jo]);
        }
    }

}

