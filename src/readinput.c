#include "readinput.h"

//#define __DEBUG
//#define __DEBUG1
//#define __DEBUG2
//#define __DEBUG_mpi

//#define __DEBUG3
//#define __DEBUG144_1

void readinput(char * fn_input,
               int * p2flag_soc,
               char * seed,
               double lattice[3][3],
               double rotations[][3][3], 
               double translations[][3], 
               int TR[], 
               vector * magmom,
               int * p2number_of_symm, 
               wannorb ** p2orb_info, 
               int * p2num_wann, 
               vec_llist ** p2kpts,
               int * p2nkpt,
               vec_llist ** p2kpaths,
               char klabels[][SHORTLEN],
               int * p2nkpath,
               int * p2nk_per_kpath,
               int * p2flag_bands,
               int * p2flag_chaeig, 
               int * p2flag_chaeig_in_kpath,
               int * p2flag_restart, 
               int * p2flag_global_trsymm, 
               int * p2flag_expandrvec, 
               double * p2symm_magnetic_tolerance,
               double * p2ham_tolerance,
               double * p2degenerate_tolerance,
               int * p2flag_everysymm,
               int * p2flag_output_mem_usage,
               int * p2flag_symm_from_file, 
               char * fn_symm)
{
    //input: 
    //          InputFile Name: fn_input
    //output:
    //          flag of soc:                    *p2flag_soc
    //          wannier hr.dat seed name:       seed
    //          rotation part of symmetry:      rotations
    //          translation part of symmetry:   translations
    //          number of symmetries:           *p2number_of_symm
    //          orbital infomation:             *p2orb_info
    //          number of wannier orbitals:     *p2num_wann
    
    double lattice4spg[3][3];                      //column major lattice basis vector for spglib to use
    int irotations[MAX_NUM_of_SYMM][3][3];
    int number_of_atomtypes;
    int * atom_types;
    double atom_positions[MAX_NUM_of_atoms][3];
    char name_of_atoms_each[128][4];
    int number_of_atoms_each[128]={0};
    int number_of_atoms_total;
    int max_num_of_symm = MAX_NUM_of_SYMM;
    int i, j, k;
    int ii,jj,kk;
    char msg[MAXLEN];

    FILE * fin;
    FILE * fout;
    FILE * fsymm;
    FILE * fstdout;
    FILE * fstderr;
    char fn[MAXLEN];
    char line[MAXLEN];
    char tag[MAXLEN]="trivial";
    char arg[MAXLEN]="trivial";
    int flag_tmp;

    double magmom_array[MAX_NUM_of_atoms*3];
    char magmom_string[MAX_NUM_of_atoms*30]="trivial";
    int  flag_magnetism=0;  // 0:without magnetism, 1:with magnetism input
    int  magnetic_type=0;   
    vector SAXIS;

    char tag_dftcode[]="dftcode";
    char tag_seedname[]="seedname";
    char tag_spinors[]="spinors";
    char tag_useposcar[]="useposcar";

    int code_type;
    int begin_tag_exist=0; // special tag, for begin structure_in_format_of_POSCAR
 
    projgroup pjgroup[MAXLEN];
    int  num_pjgroup;
    //char projgroup_element[MAXLEN][8];
    //char projgroup_orbname[MAXLEN][8];

    vector tmpkpt;
    int flag_bands_symm=1;
    int flag_bands_ori=1;

    //init SAXIS
    init_vector(&SAXIS, 0, 0, 1);
    
    fin=fopen(fn_input, "r");
    while(fgets(line, MAXLEN, fin) != NULL){
        parseline(tag, arg, line, 1);
        if( tag[0] == '\0' || strcmp(tag, "") == 0 ){
            //empty string do nothing
            continue;
        }
        else if( strcmp(tag, tag_dftcode ) == 0 ){
            setup_codetype(&code_type, arg);
        }
        else if( strcmp(tag, "seedname" ) == 0 ){
            strcpy(seed,arg);
        }
        else if( strcmp(tag, "spinors" ) == 0 ){
            if( (flag_tmp = str2boolean(arg)) != -1 )
                *p2flag_soc = flag_tmp;
            else{
                sprintf(msg, "ERROR: spinors = %s --- undefined spinors type. ", arg);
                print_error(msg);
                exit(1);
            }
        }
        else if( strcmp(tag, "useposcar" ) == 0 || strcmp(tag, "use_poscar")==0 ){
            read_pos_info(lattice, 
                          atom_positions,
                          name_of_atoms_each,
                          number_of_atoms_each,
                          &number_of_atoms_total,
                          &number_of_atomtypes,
                          arg, NULL, 0);
        }
        else if( strcmp(tag, "structure_in_format_of_poscar")==0 || strcmp(tag, "beginstructure_in_format_of_poscar" )==0){
            if(strcmp(tag, "beginstructure_in_format_of_poscar" )==0){
                begin_tag_exist=1;
            }
            read_pos_info(lattice, 
                          atom_positions,
                          name_of_atoms_each,
                          number_of_atoms_each,
                          &number_of_atoms_total,
                          &number_of_atomtypes,
                          NULL, fin, begin_tag_exist);
        }
        else if( strcmp(tag, "beginprojections")==0 || strcmp(tag, "beginprojection" )==0){
            read_projection_info(pjgroup, 
                                 &num_pjgroup, 
                                 fin);
        }
        else if( strcmp(tag, "beginkpath")==0 || strcmp(tag, "beginkpaths" )==0 || strcmp(tag, "beginkpoint_path")==0 ){
            read_kpath_info(p2kpaths, klabels, p2nkpath, fin);
        }
        else if( strcmp(tag, "nk_per_kpath")==0 || strcmp(tag, "nkperkpath")==0 ){
            sscanf(arg, "%d", p2nk_per_kpath);
            if(*p2nk_per_kpath < 2){
                *p2nk_per_kpath = 2;
            }
        }
        else if( strcmp(tag, "beginkpts")==0 || strcmp(tag, "beginkpoints" )==0 || 
                 strcmp(tag, "beginkpt")==0  || strcmp(tag, "beginkpoint" )==0 ){
            read_kpts_info(p2kpts, p2nkpt, fin);
        }
        else if( strcmp(tag, "kpt")==0){
            sscanf(arg, "%lf%lf%lf", &(tmpkpt.x), &(tmpkpt.y), &(tmpkpt.z));
            vec_llist_add(p2kpts, tmpkpt);
            (*p2nkpt)++;
        }
        else if( strcmp(tag, "bands_symmed")==0 || strcmp(tag, "bands_symmed_plot")==0 || 
                 strcmp(tag, "bands_symmetrized")==0 || strcmp(tag, "bands_symmetrized_plot")==0 ||
                 strcmp(tag, "calculatebandssymed")==0 || strcmp(tag, "calculate_bands_symmed")==0 ){
            if( (flag_tmp = str2boolean(arg)) != -1 ){
                flag_bands_symm=flag_tmp;
            } else {
                sprintf(msg, "ERROR: %s = %s --- undefined value, should be replaced by T or F\n", tag, arg);
                print_error(msg);
                exit(1);
            }
        }
        else if( strcmp(tag, "bands_ori")==0 || strcmp(tag, "bands_ori_plot")==0 || 
                 strcmp(tag, "bands_original")==0 || strcmp(tag, "bands_original_plot")==0 || 
                 strcmp(tag, "calculatebandsori")==0 || strcmp(tag, "calculate_bands_ori")==0 ){
            if( (flag_tmp = str2boolean(arg)) != -1 ){
                flag_bands_ori =flag_tmp;
            } else {
                sprintf(msg, "ERROR: %s = %s --- undefined value, should be replaced by T or F\n", tag, arg);
                print_error(msg);
                exit(1);
            }
        }
        else if( strcmp(tag, "bands")==0 || strcmp(tag, "bands_plot")==0 || 
                 strcmp(tag, "calculatebands")==0 || strcmp(tag, "calculate_bands")==0 ){
            if( (flag_tmp = str2boolean(arg)) != -1 ){
                flag_bands_symm=flag_tmp;
                flag_bands_ori =flag_tmp;
            } else {
                sprintf(msg, "ERROR: %s = %s --- undefined value, should be replaced by T or F\n", tag, arg);
                print_error(msg);
                exit(1);
            }
        }
        else if( strcmp(tag, "chaeig")==0 || strcmp(tag, "calculate_chaeig")==0 || strcmp(tag, "calculatechaeig")==0 ||
                 strcmp(tag, "charactersandeigenvalues")==0 || strcmp(tag, "calculatecharactersandeigenvalues")==0 ||
                 strcmp(tag, "characters_and_eigenvalues")==0 || strcmp(tag, "calculate_characters_and_eigenvalues")==0){
            if( (flag_tmp = str2boolean(arg)) != -1 )
                *p2flag_chaeig = flag_tmp;
            else{
                sprintf(msg, "ERROR: %s = %s --- undefined value, should be replaced by T or F\n", tag, arg);
                print_error(msg);
                exit(1);
            }
        }
        else if( strcmp(tag, "characters_in_kpath")==0 || strcmp(tag, "charactersinkpath")==0 ||
                 strcmp(tag, "chaeig_in_kpath")==0 || strcmp(tag, "chaeiginkpath")==0 ||
                 strcmp(tag, "calculate_chaeig_in_kpath")==0 || strcmp(tag, "calculatechaeiginkpath")==0 ||
                 strcmp(tag, "calc_band_symm")==0 || strcmp(tag, "calcbandsymm")==0){
            if( (flag_tmp = str2boolean(arg)) != -1 ){
                *p2flag_chaeig_in_kpath = flag_tmp;
            }
            else{
                sprintf(msg, "ERROR: %s = %s --- undefined value, should be replaced by T or F\n", tag, arg);
                print_error(msg);
                exit(1);
            }
            if(*p2flag_chaeig_in_kpath == 1){
                *p2flag_chaeig = 1;
            }
        }
        else if( strcmp(tag, "restart")==0){
            if( (flag_tmp = str2boolean(arg)) != -1 )
                *p2flag_restart = flag_tmp;
            else{
                sprintf(msg, "ERROR: %s = %s --- undefined value, should be replaced by T or F\n", tag, arg);
                print_error(msg);
                exit(1);
            }
        }
        else if( strcmp(tag, "trsymm") == 0 || strcmp(tag, "global_trsymm") == 0 ){
            if( (flag_tmp = str2boolean(arg)) != -1 )
                *p2flag_global_trsymm = flag_tmp;
            else{
                sprintf(msg, "ERROR: %s = %s --- undefined value, should be replaced by T or F\n", tag, arg);
                print_error(msg);
                exit(1);
            }
        }
        else if( strcmp(tag, "expandrvec") == 0){
            if( (flag_tmp = str2boolean(arg)) != -1 )
                *p2flag_expandrvec = flag_tmp;
            else{
                sprintf(msg, "ERROR: %s = %s --- undefined value, should be replaced by T or F\n", tag, arg);
                print_error(msg);
                exit(1);
            }
        }
        else if( strcmp(tag, "symmmagnetictolerance") == 0 || strcmp(tag, "symm_magnetic_tolerance")==0){
             sscanf( arg, "%lf", p2symm_magnetic_tolerance );
        }
        else if( strcmp(tag, "hamtolerance") == 0 || strcmp(tag, "ham_tolerance") == 0){
             sscanf( arg, "%lf", p2ham_tolerance );
        }
        else if( strcmp(tag, "degeneratetolerance") == 0 || strcmp(tag, "degenerate_tolerance") == 0){
             sscanf( arg, "%lf", p2degenerate_tolerance );
        }
        else if( strcmp(tag, "usesymmetry") == 0 || strcmp(tag, "use_symmetry") == 0 || 
                 strcmp(tag, "symminputfile") == 0 || strcmp(tag, "symm_input_file") == 0){
            *p2flag_symm_from_file = 1;
            strcpy(fn_symm, arg);
        }
        else if( strcmp(tag, "everysymm")==0){
            if( (flag_tmp = str2boolean(arg)) != -1 )
                *p2flag_everysymm = flag_tmp;
            else{
                sprintf(msg, "ERROR: %s = %s --- undefined value, should be replaced by T or F\n", tag, arg);
                print_error(msg);
                exit(1);
            }
        }
        else if( strcmp(tag, "outputmemusage")==0 || strcmp(tag, "output_mem_usage")==0 ){
            if( (flag_tmp = str2boolean(arg)) != -1 )
                *p2flag_output_mem_usage = flag_tmp;
            else{
                sprintf(msg, "ERROR: %s = %s --- undefined value, should be replaced by T or F\n", tag, arg);
                print_error(msg);
                exit(1);
            }
        }
        else if( strcmp(tag, "magmom") == 0 ){
            flag_magnetism = 1;
            if( strcmp(arg,"NULL") == 0){
                sprintf(msg, "WARNING: MAGMOM is set to NULL");
                print_error(msg);
                flag_magnetism = 0;
            }
            else{
                strcpy(magmom_string, arg);
            }
        }
        else if( strcmp(tag, "saxis") == 0 ){
            sscanf(arg, "%lf%lf%lf", &(SAXIS.x), &(SAXIS.y), &(SAXIS.z));
        }
        else if( strcmp(tag, "symmtolerance") == 0 || strcmp(tag, "symm_tolerance") == 0 ){
            sprintf(msg, "WARNING: This TAG (\"%s\") will be abandoned, please use \"ham_tolerance\"\n", tag);
            print_error(msg);
            if( *p2ham_tolerance == 0.1 ){ // only give the new value if it is unchanged
                sscanf( arg, "%lf", p2ham_tolerance );
            }
        }
        else if( strcmp(tag, "NULL") != 0 ){
            sprintf(msg, "WARNING: UNKNOWN TAG: \"%s\"\n", tag);
            print_error(msg);
        }
    }
    fclose(fin);

    // if no k-point or k-path specified, no need to calculated bands, chas and eigs
    if(*p2nkpath == 0) {
        flag_bands_symm = 0;
        flag_bands_ori  = 0;
        *p2flag_chaeig_in_kpath = 0;
    }
    *p2flag_bands = flag_bands_ori * 2 + flag_bands_symm;
    if(*p2nkpt == 0 && (*p2nkpath == 0 || *p2flag_chaeig_in_kpath == 0)){
        *p2flag_chaeig = 0;
    }

    //---- derive orbital infomation from projections
    derive_projection_info(p2num_wann, 
                           p2orb_info, 
                           pjgroup,
                           num_pjgroup, 
                           lattice, 
                           atom_positions, 
                           name_of_atoms_each, 
                           number_of_atoms_each, 
                           number_of_atoms_total, 
                           number_of_atomtypes,
                           code_type,
                           *p2flag_soc);

    //---- derive MAGMOM from magmom_string
    if(flag_magnetism!=0){
        derive_magmom_from_string(magmom_array, magmom_string, number_of_atoms_total, *p2flag_soc);
        derive_magmom_from_array(magmom, magmom_array, number_of_atoms_total, *p2flag_soc, SAXIS);
        magnetic_type = get_magnetic_type(magmom, number_of_atoms_total, *p2symm_magnetic_tolerance);

        if(magnetic_type!=0){
            *p2flag_global_trsymm = 0;
        }
    }

    //---- derive symmetry operation in real space (use lattice vector as basis) with spglib
    for(i=0;i<3;i++)
        for(j=0;j<3;j++)
            lattice4spg[i][j] = lattice[j][i];

    atom_types=(int *) malloc(sizeof(int)*number_of_atoms_total);
    k=0;
    for(i=0;i<number_of_atomtypes;i++){
            for(j=0;j<number_of_atoms_each[i];j++){
                    atom_types[k++]=i+1; // the (i+1)th type of atoms
            }
    }

    *p2number_of_symm = spg_get_symmetry(irotations, translations, max_num_of_symm,
                                      lattice4spg, atom_positions, atom_types,
                                      number_of_atoms_total, 1E-5);

    for(i=0;i<*p2number_of_symm;i++){
        for(j=0;j<3;j++){
            for(k=0;k<3;k++){
                rotations[i][j][k] = (double)irotations[i][j][k];
            }
        }
        TR[i] = -1; // initialize TRS to NULL
    }

    //derive symmetry for magnetic materials from the symmetry found by spglib
    if(flag_magnetism !=0 )
    {
        derive_symm_for_magnetic_materials(rotations,
                                           translations,
                                           TR,
                                           p2number_of_symm,
                                           lattice,
                                           atom_positions,
                                           atom_types,
                                           number_of_atoms_total,
                                           magmom,
                                           *p2symm_magnetic_tolerance);
    }


    //---- output configuration & symmetries
    fstdout=fopen("wannsymm.out", "a");
    fprintf(fstdout, "seedname      = %s\n", seed);
    fprintf(fstdout, "restart state = " );
    if(*p2flag_restart==1){
        fprintf(fstdout,"True. (Directly read %s_symmed_hr.dat and skip the symmetrization procedure)\n",seed);
    }
    else if(*p2flag_restart==0){
        fprintf(fstdout,"False. (Read %s_hr.dat and symmetrize it)\n",seed);
    }
    fprintf(fstdout, "DFT_code_type = %d\n", code_type);
    print_flag(fstdout, "flag_soc                   ", *p2flag_soc);
    print_flag(fstdout, "flag_symm_from_file        ", *p2flag_symm_from_file);
    //if(*p2flag_symm_from_file == 1){
    //    fprintf(fstdout,"symmetry file name:%s\n", fn_symm);
    //}
    print_flag(fstdout, "flag_consider_global_trsymm", *p2flag_global_trsymm);
    print_flag(fstdout, "flag_expandrvec            ", *p2flag_expandrvec);
    if( *p2flag_expandrvec == 0){
        fprintf(fstdout,"    ham_tolerance           = %.5lf\n", *p2ham_tolerance);
    }
    print_flag(fstdout, "flag_calculate_chaeig      ", *p2flag_chaeig);
    print_flag(fstdout, "flag_chaeig_in_kpath       ", *p2flag_chaeig_in_kpath);
    fprintf(fstdout,    "num of k-point              = %d (for cha and eig calculation)\n", *p2nkpt);
    print_flag(fstdout, "flag_bands_symmetrized     ", flag_bands_symm);
    print_flag(fstdout, "flag_bands_original        ", flag_bands_ori);
    fprintf(fstdout,    "num of kpt per kpath        = %d\n", *p2nk_per_kpath);
    fprintf(fstdout,    "num of k-path               = %d (for bands calculation)\n", *p2nkpath);
    fprintf(fstdout,    "num of atomtypes            = %d\n",number_of_atomtypes);
    fprintf(fstdout,    "name_of_atoms               = ");
    for(i=0;i<number_of_atomtypes;i++){
            fprintf(fstdout,"%-6s",name_of_atoms_each[i]);
    }
    fprintf(fstdout,"\n");
    fprintf(fstdout,    "number_of_atoms_each        = ");
    for(i=0;i<number_of_atomtypes;i++){
        fprintf(fstdout,"%-6d",number_of_atoms_each[i]);
    }
    fprintf(fstdout,"\n");
    //fprintf(fstdout,"atom_types(for spglib)=");
    //for(i=0;i<number_of_atoms_total;i++){
    //        fprintf(fstdout,"%3d ", atom_types[i]);
    //        if(i!=0 && i%10==0){
    //                fprintf(fstdout,"\n                       ");
    //        }
    //}
    //fprintf(fstdout,"\n");
    fprintf(fstdout,    "number of symmetries        = %d\n", *p2number_of_symm);
    if(*p2flag_symm_from_file != 1){
        fprintf(fstdout,"The symmetry operations (in real space, using lattice vectors as basis) are in this file: symmetries.dat\n");
    } else {
        fprintf(fstdout,"The symmetry operations are read from this file: %s\n", fn_symm);
    }
    fprintf(fstdout,"\n");
    fprintf(fstdout,"lattice:\n");
    for(i=0;i<3;i++){
        for(j=0;j<3;j++)
            fprintf(fstdout,"%9.5lf", lattice[i][j]);
        fprintf(fstdout,"\n");
    }
    fprintf(fstdout,"\n");

    fprintf(fstdout,"atoms' coordinates (Direct):\n");
    for(i=0;i<number_of_atoms_total;i++){
        for(j=0;j<3;j++)
            fprintf(fstdout,"%9.5lf", atom_positions[i][j]);
        fprintf(fstdout,"\n");
    }
    fprintf(fstdout,"\n");
    fclose(fstdout);
    if(flag_magnetism != 0){
        sprintf(msg, "magnetization:\n");
        print_msg(msg);
        sprintf(msg, "# of ion      x         y         z  \n");
        print_msg(msg);
        sprintf(msg, "-----------------------------------------\n");
        print_msg(msg);
        for(ii=0;ii<number_of_atoms_total; ii++){
            sprintf(msg, "%5d   %10.5lf%10.5lf%10.5lf\n", ii+1, magmom[ii].x, magmom[ii].y, magmom[ii].z);
            print_msg(msg);
        }
        sprintf(msg, "\n");
        print_msg(msg);
        sprintf(msg,"symm_magnetic_tolerance = %E\n", *p2symm_magnetic_tolerance);
        print_msg(msg);
        sprintf(msg, "magnetic order type: ");
        print_msg(msg);
        switch(magnetic_type) {
            case -1: 
                sprintf(msg, "anti-ferromagnetic\n");
                break;
            case 0: 
                sprintf(msg, "non-magnetic\n");
                break;
            case 1: 
                sprintf(msg, "other magnetic type (undefined in the program)\n");
                break;
            case 2: 
                sprintf(msg, "ferromagnetic or ferrimagnetic\n");
                break;
            default:
                sprintf(msg, "magnetic type undefined (in the program)\n");
        }
        print_msg(msg);
        sprintf(msg, "\n");
        print_msg(msg);
    }

    //Use spglib to derive symmetry group number and symbol
    int spg_num_international;    //space group number
    int spg_num_schoenflies;      //space group number
    char spg_symbol_international[11], spg_symbol_schoenflies[7];
    double origin_shift[3]={0,0,0};
    if(*p2flag_symm_from_file != 1){
        fsymm=fopen("symmetries.dat", "w");
        spg_num_international = spg_get_international(spg_symbol_international, lattice4spg, atom_positions, atom_types, number_of_atoms_total, 1e-5);
        spg_num_schoenflies = spg_get_schoenflies(spg_symbol_schoenflies, lattice4spg, atom_positions, atom_types, number_of_atoms_total, 1e-5);
        fprintf(fsymm, "space group infomation:\n");
        fprintf(fsymm, "    International: %s (%d)\n",spg_symbol_international, spg_num_international);
        fprintf(fsymm, "    schoenflies: %s (%d)\n",spg_symbol_schoenflies, spg_num_schoenflies);
        print_flag(fsymm, "global time-reversal symmetry", *p2flag_global_trsymm);
        if(magnetic_type != 0){
            fprintf(fsymm, "magnetic order detected\n");
            fprintf(fsymm, "symmetries in the corresponding magnetic group:\n");
        }
        fprintf(fsymm, "nsymm = %d\n", *p2number_of_symm);
        for (i = 0; i < *p2number_of_symm; i++){
            fprintf(fsymm,"--- %d ---\n",i+1); // number of symmetry
            for (j = 0; j < 3; j++){
                    fprintf(fsymm,"%2d %2d %2d\n", irotations[i][j][0],
                                                   irotations[i][j][1],
                                                   irotations[i][j][2]);
            }
            fprintf(fsymm, "%f %f %f", translations[i][0],
                                        translations[i][1],
                                        translations[i][2]);
            if(*p2flag_global_trsymm == 0){
                if(TR[i] == 1) {
                    fprintf(fsymm, " T");
                } else {
                    fprintf(fsymm, " F");
                }
            }
            fprintf(fsymm,"\n");
        }
        fclose(fsymm);
    }
    free(atom_types);
}

void parsestr_find_remove( char strout[MAXLEN], char strin[MAXLEN], char target){
    // remove $target from strin, output is strout, strin and strout can be same
    char tmpout[MAXLEN];
    strcpy(tmpout,"NULL");
    int i;
    int j=0;

    for(i=0; strin[i]!='\0' && strin[i] != '\n' && i<MAXLEN; i++){
        if( strin[i]==target){
            continue;
        }
        tmpout[j++]=strin[i];
    }
    tmpout[j]='\0';
    strcpy(strout, tmpout);
}

void parseline( char tag[MAXLEN], char arg[MAXLEN], char line[MAXLEN], int ignorecase){
    strcpy(tag,"NULL");
    strcpy(arg,"NULL");
    int metdivid=0;         //met '=' or ':'
    int metarg=0;           //met first nonspace argument after '=' or ':'
    int flag_inquota=0;
    int i;
    int j=0;
    int k=0;

    for(i=0; line[i]!='\0' && line[i]!='\n' && i<MAXLEN; i++){
        if( flag_inquota==0 && (line[i]=='#' || line[i]=='!' || ( line[i]=='/' && line[i+i]=='/'))){
            break;
        }
        if( line[i]=='\'' || line[i]==34){
            flag_inquota = 1-flag_inquota;
            continue;
        }
        if( flag_inquota==0 && line[i]=='=' || line[i]==':'){
            metdivid++;
            continue;
        }
        if( metarg==0 && line[i] == ' ') continue;
        if( flag_inquota==0  && metdivid==0 && ignorecase !=0 && line[i]>=65 && line[i]<=90){
            line[i]+=32;    //lower-case
        }
        if( flag_inquota==0 && line[i] <= 31 ) continue;    //Control characters, include '\t'
        //if( flag_inquota==0 && line[i] ==42 ) continue;   // '*'
        //if( flag_inquota==0 && line[i] ==47 ) continue;   // '/'
        //if( flag_inquota==0 && line[i] ==58 ) continue;   // ':'
        //if( flag_inquota==0 && line[i] >=60 && line[i] <=64 ) continue; // '<' '=' '>' '?' '@'
        //if( flag_inquota==0 && line[i] >=91 && line[i] <=96 ) continue; // '[' '\' ']' '^' '_' '`'
        if( flag_inquota==0 && line[i] >=127) continue; 

        if( metdivid ==0){
            tag[j++] = line[i];
        }else if( metdivid > 0){
            arg[k++] = line[i];
            if(line[i] != ' ') metarg++;
        }else{
            fprintf(stderr,"ERROR: parseline() at readinput.c");
            exit(1);
        }
    }
    if(flag_inquota==1){
        fprintf(stderr, "quotation marks not complete in input file\n");
        exit(1);
    }
    if(j!=0)  tag[j] = '\0';
    if(k!=0){
        arg[k] = '\0';
        while(arg[k-1]==' '){
            arg[ k - 1 ] = '\0';
            k--;
        }
    }
#ifdef __DEBUG
        printf( "tag:%s. arg:%s.\n", tag, arg);
#endif
}


void setup_codetype(int * p2code_type, char dftcode[MAXLEN]){
    if( strcmp(dftcode, "VASP") == 0 ){
        *p2code_type=1;
    } else if( strcmp(dftcode, "QE") == 0){
        *p2code_type=2;
        fprintf(stderr, "\n\nNote: DFTcode=QE is a testing feature\n\n");
    } else {
        fprintf(stderr, "ERROR!!!: not supported code type, %s,", dftcode);
        fprintf(stderr, "Only supporting VASP now\n");
        exit(1);
    }
}

void read_pos_info(double lattice[3][3], 
                   double atom_positions[MAX_NUM_of_atoms][3],
                   char   name_of_atoms_each[128][4],
                   int    number_of_atoms_each[128],
                   int *  p2number_of_atoms_total,
                   int *  p2number_of_atomtypes,
                   char   fn_pos[MAXLEN],
                   FILE * fin,
                   int begin_tag_exist){
    FILE * fin_pos;
    char line[MAXLEN];
    char tag[MAXLEN];
    char arg[MAXLEN];
    char msg[MAXLEN];
    int i,j,k;
    double scaling;

    if( fn_pos!=NULL){
        if( ! file_exists(fn_pos)){
            sprintf(msg,"ERROR!!! the file \"%s\" can not be found\n", fn_pos);
            print_error(msg);
            exit(1);
        }
        fin_pos=fopen(fn_pos, "r");
        if(fin_pos==NULL){
            sprintf(msg,"ERROR!!! unable to open FILE: %s.\n", fn_pos);
            print_error(msg);
            exit(1);
        }
    }
    else if(fin!=NULL){
        fin_pos=fin;
    }
    else{
        sprintf(msg,"BUG: can not find any position info\n");
        print_error(msg);
        exit(1);
    }

    fgets(line, MAXLEN, fin_pos);
    parseline(tag, arg, line, 0);
    while( isletter(tag[0]) ){ 
        // if the first line in POSCAR is name, skip first line
        fgets(line, MAXLEN, fin_pos);
        parseline(tag, arg, line, 0);
    }
    sscanf(line, "%lf", &scaling);      // scale
    for(i=0;i<3;i++){
        fgets(line, MAXLEN, fin_pos);   // lattice
        sscanf(line, "%lf%lf%lf", &lattice[i][0], &lattice[i][1], &lattice[i][2]);
        lattice[i][0] *= scaling; lattice[i][1] *= scaling; lattice[i][2] *= scaling;
    }

    j=0;k=0;
    fgets(line, MAXLEN, fin_pos);       //name of atoms
    for(i=0; line[i]!='\0' && line[i]!='\n'; i++){
        if( ! isletter(line[i]) ){
            if(i==0) j--;
            continue;
        }
        else {
            if( i!=0 && ! isletter(line[i-1])){
                name_of_atoms_each[j][k] = '\0';
                j++;
                k=0;
            }
            name_of_atoms_each[j][k++] = line[i];
        }
    }
    name_of_atoms_each[j][k] = '\0';
    if(j==-1 || (j==0 && k==0)){
        fprintf(stderr, "ERROR!!!: The 6th line of POSCAR must contain elements' name.\n");
        exit(1);
    }

    j=0;
    fgets(line, MAXLEN, fin_pos);       //number of atoms
    for(i=0;line[i] != '\0'; i++){
        if( line[i]<48 || line[i]>57  ){
            if(i==0) j--;
            continue;
        }
        else {
            if( i!=0 && (line[i-1]<48 || line[i-1]>57)){
                    j++;
            }
            number_of_atoms_each[j] = number_of_atoms_each[j] * 10 + (line[i] - 48);
        }
    }

    *p2number_of_atomtypes=j+1;
    *p2number_of_atoms_total=0;
    for(i=0;i<*p2number_of_atomtypes;i++){
        *p2number_of_atoms_total += number_of_atoms_each[i];
    }

    int flag_cart=0;

    fgets(line, MAXLEN, fin_pos);          //Direct or Cartisian or Seletive dynamics
    parseline(tag, arg, line, 0);
    if( tag[0] == 'S'){
        fgets(line, MAXLEN, fin_pos);       
        parseline(tag, arg, line, 0);
    }
    if( tag[0] == 'D' || tag[0] == 'd') flag_cart=0;
    else if( tag[0] == 'C' || tag[0] == 'c' || tag[0] == 'K' || tag[0] == 'k') flag_cart=1;

    for(i=0;i<*p2number_of_atoms_total;i++){
        fgets(line, MAXLEN, fin_pos);   //atoms' position
        sscanf(line, "%lf%lf%lf",&atom_positions[i][0],
                                 &atom_positions[i][1],
                                 &atom_positions[i][2]);
    }
    if(flag_cart==1){
        sprintf(msg, "ERROR!!! Cartesian not support now, use Direct instead.\n Cartesian will be support in future.\n");
        print_error(msg);
        exit(1);
    }
    if( fn_pos!=NULL){
        // skip other lines (such as lines of initial velocities for MD)
        fclose(fin_pos);
    } else if(fin!=NULL && begin_tag_exist == 1){
        // skip lines before "end"
        parseline(tag, arg, line, 1);
        while( ! (tag[0] == 'e' && tag[1] == 'n' && tag[2] == 'd')){
            if(fgets(line, MAXLEN, fin_pos)!=NULL){
                parseline(tag, arg, line, 1);
            }
            else{
                sprintf(msg, "ERROR in intputfile: Can not find \"end\" tag for \"(begin) Structure_in_format_of_POSCAR\"\n");
                print_error(msg);
                exit(1);
            }
        }
    }
}

void read_projection_info(projgroup pjgroup[MAXLEN],
                          int *  p2num_pjgroup,
                          FILE * fin){
    int i,j,k;
    int ii,jj,kk;
    char line[MAXLEN];
    char tag[MAXLEN]="trivial initial string";
    char tag_lowercase[MAXLEN]="trivial initial string";
    char arg[MAXLEN];
    char arg_lowercase[MAXLEN];

    i=0;j=0;k=0;
    fgets(line, MAXLEN, fin);
    parseline(tag, arg, line, 0);
    parsestr_find_remove(arg, arg, ' ');
    parseline(tag_lowercase, arg_lowercase, line, 1);
    parsestr_find_remove(arg_lowercase, arg_lowercase, ' ');
    for(i=0;strcmp(tag_lowercase, "endprojections")!=0 && strcmp(tag_lowercase, "endprojection")!=0; i++){
        strcpy(pjgroup[i].element, tag);
        for(ii=0; arg[ii] != '\0'; ii++){
            if(arg[ii]==',' || arg[ii]==59){
                pjgroup[i].orbname[j++][k] = '\0';
                k=0;
            }
            else{
                pjgroup[i].orbname[j][k++] = arg[ii];
            }
        }
        if( !isletter( pjgroup[i].orbname[j][0]) ){
            fprintf(stderr, "ERROR!!! projections defined wrong.\n");
            exit(1);
        }
        pjgroup[i].orbname[j][k] = '\0';
        pjgroup[i].npj=j+1;
        j=k=0;
        fgets(line, MAXLEN, fin);
        parseline(tag, arg, line, 0);
        parsestr_find_remove(arg, arg, ' ');
        parseline(tag_lowercase, arg_lowercase, line, 1);
        parsestr_find_remove(arg_lowercase, arg_lowercase, ' ');
    }
    *p2num_pjgroup = i;
#ifdef __DEBUG
    for(i=0;i<*p2num_pjgroup;i++){
        printf("orb define: %s : ", pjgroup[i].element);
        for(j=0;j<pjgroup[i].npj;j++)
            printf("%s, ", pjgroup[i].orbname[j]);
        printf("\n");
    }
#endif
}

void derive_projection_info(int *  p2num_wann, 
                            wannorb ** p2orb_info,
                            projgroup pjgroup[],
                            int    num_pjgroup,
                            double lattice[3][3],
                            double atom_positions[MAX_NUM_of_atoms][3],
                            char   name_of_atoms_each[128][4], 
                            int    number_of_atoms_each[128], 
                            int    number_of_atoms_total, 
                            int    number_of_atomtypes,
                            int    code_type,
                            int    flag_soc){

    int i,j,k;
    int ii,jj,kk;
    int isite, iproj;
    int site_first;
    vector vx,vy,vz;
    vector site;
    int l, mr, r;

    r=1;    // set to 1 by deafult
    init_vector(&vz, 0, 0, 1);  //  Local orbital set to default now,
    init_vector(&vx, 1, 0, 0);  //  will add support for self defined 
    vy=cross_product(vz, vx);
    //init_vector(&vy, 0, 1, 0);  //  local orb someday.

    *p2num_wann=0;

    for(i=0;i<num_pjgroup;i++){
        for(ii=0;ii<number_of_atomtypes;ii++){
            if(strcmp(pjgroup[i].element, name_of_atoms_each[ii])==0)
                break;
        }
        for(j=0;j<pjgroup[i].npj; j++){
            if(      pjgroup[i].orbname[j][1] != '\0') *p2num_wann += 1*number_of_atoms_each[ii];
            else if( pjgroup[i].orbname[j][0] == 's')  *p2num_wann += 1*number_of_atoms_each[ii];
            else if( pjgroup[i].orbname[j][0] == 'p')  *p2num_wann += 3*number_of_atoms_each[ii];
            else if( pjgroup[i].orbname[j][0] == 'd')  *p2num_wann += 5*number_of_atoms_each[ii];
            else if( pjgroup[i].orbname[j][0] == 'f')  *p2num_wann += 7*number_of_atoms_each[ii];
        }
    }
    if(flag_soc==1) *p2num_wann *= 2;
    *p2orb_info=(wannorb *) malloc(sizeof(wannorb)*(*p2num_wann));


    iproj=0;
    for(i=0;i<num_pjgroup;i++){
        site_first=0;
        for(ii=0;ii<number_of_atomtypes;ii++){
            if(strcmp(pjgroup[i].element, name_of_atoms_each[ii])==0)
                break;
            site_first += number_of_atoms_each[ii];
        }
        for(isite=site_first; isite < (site_first+number_of_atoms_each[ii]); isite++){
            init_vector(&site, atom_positions[isite][0], 
                               atom_positions[isite][1], 
                               atom_positions[isite][2]);
            for(j=0; j < pjgroup[i].npj; j++){
                if( pjgroup[i].orbname[j][1] == '\0'){
                    if(      strcmp(pjgroup[i].orbname[j],"s")==0) l=0;
                    else if( strcmp(pjgroup[i].orbname[j],"p")==0) l=1;
                    else if( strcmp(pjgroup[i].orbname[j],"d")==0) l=2;
                    else if( strcmp(pjgroup[i].orbname[j],"f")==0) l=3;
                    else{
                        fprintf(stderr, "unknown orbital definition\n");
                        exit(1);
                    }
                    for(mr=1;mr<(2*l+2);mr++){
                        init_wannorb((*p2orb_info + (iproj++)), site, l, mr, 0, r, vz, vx);
                        if(flag_soc==1 && code_type!=1)
                            init_wannorb((*p2orb_info + (iproj++)), site, l, mr, 1, r, vz, vx);
#ifdef __DEBUG
                        printf("ln487: orbinfo: %10.5lf%10.5lf%10.5lf%5d%5d%5d\n", site.x,site.y,site.z, l, mr, r);
#endif
                    }
                }else{
                    if(      strcmp(pjgroup[i].orbname[j],"pz")==0){ l=1; mr=1;}
                    else if( strcmp(pjgroup[i].orbname[j],"px")==0){ l=1; mr=2;}
                    else if( strcmp(pjgroup[i].orbname[j],"py")==0){ l=1; mr=3;}

                    else if( strcmp(pjgroup[i].orbname[j],"dz2")==0){    l=2; mr=1;}
                    else if( strcmp(pjgroup[i].orbname[j],"dxz")==0){    l=2; mr=2;}
                    else if( strcmp(pjgroup[i].orbname[j],"dyz")==0){    l=2; mr=3;}
                    else if( strcmp(pjgroup[i].orbname[j],"dx2-y2")==0){ l=2; mr=4;}
                    else if( strcmp(pjgroup[i].orbname[j],"dxy")==0){    l=2; mr=5;}

                    else if( strcmp(pjgroup[i].orbname[j],"fz3")==0){        l=3; mr=1;}
                    else if( strcmp(pjgroup[i].orbname[j],"fxz2")==0){       l=3; mr=2;}
                    else if( strcmp(pjgroup[i].orbname[j],"fyz2")==0){       l=3; mr=3;}
                    else if( strcmp(pjgroup[i].orbname[j],"fz(x2-y2)")==0){  l=3; mr=4;}
                    else if( strcmp(pjgroup[i].orbname[j],"fxyz")==0){       l=3; mr=5;}
                    else if( strcmp(pjgroup[i].orbname[j],"fx(x2-3y2)")==0){ l=3; mr=6;}
                    else if( strcmp(pjgroup[i].orbname[j],"fy(3x2-y2)")==0){ l=3; mr=7;}
                    else{
                        fprintf(stderr, "wrong defined orbital: \"%s\"\n", pjgroup[i].orbname[j] );
                        exit(1);
                    }

                    init_wannorb((*p2orb_info+(iproj++)), site, l, mr, 0, r, vz, vx);
                    if(flag_soc==1 && code_type!=1)
                        init_wannorb((*p2orb_info+(iproj++)), site, l, mr, 1, r, vz, vx);
                }
            }
        }
    }
    if(flag_soc==1 && code_type==1){
        if( iproj*2 != *p2num_wann){
            fprintf(stderr, "ERROR!!! error when deriving orbital info at readinput.c\n");
            exit(1);
        }
        for(iproj=0;iproj<*p2num_wann/2;iproj++){
            site = (*p2orb_info+iproj)->site;
            l    = (*p2orb_info+iproj)->l;
            mr   = (*p2orb_info+iproj)->mr;
            r    = (*p2orb_info+iproj)->r;
            vz   = (*p2orb_info+iproj)->axis[2];
            vx   = (*p2orb_info+iproj)->axis[0];
            init_wannorb((*p2orb_info + *p2num_wann/2 + iproj), site, l, mr, 1, r, vz, vx);
        }
    }
}
    

static void show_spg_dataset(double lattice[3][3],
			     const double origin_shift[3],
			     double position[][3],
			     const int num_atom,
			     const int types[],
                 FILE * fout){
  SpglibDataset *dataset;
  char ptsymbol[6];
  int pt_trans_mat[3][3];

  int i, j, size;
  const char *wl = "abcdefghijklmnopqrstuvwxyz";

  for ( i = 0; i < num_atom; i++ ) {
    for ( j = 0; j < 3; j++ ) {
      position[i][j] += origin_shift[j];
    }
  }

  dataset = spg_get_dataset(lattice,
			    position,
			    types,
			    num_atom,
			    1e-5);
  
  fprintf(fout,"International: %s (%d)\n", dataset->international_symbol, dataset->spacegroup_number );
  fprintf(fout,"Hall symbol:   %s\n", dataset->hall_symbol );
  spg_get_pointgroup(ptsymbol,
		     pt_trans_mat,
		     dataset->rotations,
		     dataset->n_operations);
  fprintf(fout,"Point group:   %s\n", ptsymbol);
  fprintf(fout,"Transformation matrix:\n");
  for ( i = 0; i < 3; i++ ) {
    fprintf(fout,"%f %f %f\n",
	   dataset->transformation_matrix[i][0],
	   dataset->transformation_matrix[i][1],
	   dataset->transformation_matrix[i][2]);
  }
  fprintf(fout,"Wyckoff letters:\n");
  for ( i = 0; i < dataset->n_atoms; i++ ) {
    fprintf(fout,"%c ", wl[dataset->wyckoffs[i]]);
  }
  fprintf(fout,"\n");
  fprintf(fout,"Equivalent atoms:\n");
  for (i = 0; i < dataset->n_atoms; i++) {
    fprintf(fout,"%d ", dataset->equivalent_atoms[i]);
  }
  fprintf(fout,"\n");
  
  for (i = 0; i < dataset->n_operations; i++) {
    fprintf(fout,"--- %d ---\n", i + 1);
    for (j = 0; j < 3; j++) {
      fprintf(fout,"%2d %2d %2d\n",
	     dataset->rotations[i][j][0],
	     dataset->rotations[i][j][1],
	     dataset->rotations[i][j][2]);
    }
    fprintf(fout,"%f %f %f\n",
	   dataset->translations[i][0],
	   dataset->translations[i][1],
	   dataset->translations[i][2]);
  }

  spg_free_dataset(dataset);

}

int isletter(char candidate){
    if(candidate>=65 && candidate<=90){
            return 1;
    }
    else if(candidate>=97 && candidate<=122){
            return 1;
    }
    else{
            return 0;
    }
}

//vector find_pos(double pos[MAX_NUM_of_atoms][3], double lattice[3][3], int natoms, vector site){
//    //find the exact position writed in POSCAR, with the site of low accuracy 
//    //return the position closest to the input site
//    int i,j,k;
//    double eps;
//    double dis, dis_min;
//    vector lattice_v[3];
//    vector pos_v;
//    vector pos_v_result;
//    for(i=0;i<3;i++){
//        init_vector(&lattice_v[i], lattice[i][0], lattice[i][1], lattice[i][2]);
//    }
//    dis_min=10;
//    for(i=0;i<natoms;i++){
//        init_vector(&pos_v, pos[i][0], pos[i][1], pos[i][2]);
//        dis = distance(pos_v, site, lattice_v);
//        if(dis < dis_min){
//            dis_min = dis;
//            pos_v_result=pos_v;
//        }
//    }
//    return pos_v_result;
//}



void derive_magmom_from_string(double * magmom_array, char * magmom_string, int number_of_atoms_total, int flag_soc){
    char * token;
    int i=0,j;
    int magmom_len;
    FILE * ferr;
    char rep_s[50], mag_s[50];
    int  rep;
    double mag;
    char * magtmp;
    char msg[1024];

    token = strtok(magmom_string, " ");
    while(token != NULL){
        if( (magtmp = strchr(token, '*')) != NULL){
            *magtmp = '\0';
            strcpy(rep_s, token);
            sscanf(rep_s, "%d", &rep);
            strcpy(mag_s, magtmp+1);
            sscanf(mag_s, "%lf", &mag);
            for(j=0;j<rep;j++){
                magmom_array[i++] = mag;
            }
        }
        else{
            sscanf(token, "%lf", magmom_array+i++ );
        }
        token = strtok(NULL, " ");
    }
    magmom_len = i;
    if( magmom_len != (number_of_atoms_total*(1+flag_soc*2)) ){
        sprintf(msg, "ERROR: length of MAGMOM mismatch with the atoms. MAGMOM length in file: %d, while needed length: %d\n", 
                magmom_len, number_of_atoms_total*(1+flag_soc*2) );
        print_error(msg);
        exit(1);
    }
}

void derive_magmom_from_array(vector * magmom, double * magmom_array,int natom, int flag_soc, vector SAXIS){
    int ii,jj,kk,i,j,k;
    double alpha, beta;
    vector maxis;       // magnetic moment under the axis of SAXIS
    if(SAXIS.x<1E-10){
        if(SAXIS.y<1E-10){
            alpha = 0;
        }
        else{
            alpha = PI/2;
        }
    }
    else{
        alpha = atan(SAXIS.y/SAXIS.x);
    }
    if(SAXIS.z<1E-10){
        beta = PI/2;
    }
    else{
        beta = atan(sqrt(SAXIS.x*SAXIS.x + SAXIS.y*SAXIS.y)/SAXIS.z);
    }
    for(ii=0;ii<natom;ii++){
        if(flag_soc==0){
            maxis.x=0;
            maxis.y=0;
            maxis.z=magmom_array[ii];
        }
        else{
            maxis.x=magmom_array[3*ii+0];
            maxis.y=magmom_array[3*ii+1];
            maxis.z=magmom_array[3*ii+2];
        }
        magmom[ii].x = cos(beta)*cos(alpha)*maxis.x - sin(alpha)*maxis.y + sin(beta)*cos(alpha)*maxis.z;
        magmom[ii].y = cos(beta)*sin(alpha)*maxis.x + cos(alpha)*maxis.y + sin(beta)*cos(alpha)*maxis.z;
        magmom[ii].z =           -sin(beta)*maxis.x                      +            cos(beta)*maxis.z;
    }
}

int get_magnetic_type(vector magmom[], int natom , double symm_magnetic_tolerance)
{
    // 0: non-magnetic
    // 1: other magnetic type
    // 2: FM
    // -1: AFM, spiral, SDW etc.
    int i;
    vector origin;
    //vector last_nonzero;
    int count_nonzero;
    //int count_same_direction;
    vector sumall;  // total magnetic moment
    
    init_vector(&origin, 0,0,0);
    init_vector(&sumall, 0,0,0);
    count_nonzero=0;
    for(i=0;i<natom;i++){
        sumall = vector_add(sumall, magmom[i]);
        if( distance(magmom[i],origin,NULL) > symm_magnetic_tolerance){
            count_nonzero++; // with magnetism
        }
    }
    if( distance(sumall, origin, NULL) < symm_magnetic_tolerance){
        // to be more accurate, the tolerance here can be changed to
        // symm_magnetic_tolerance * Num_of_magnetic_atom
        if( count_nonzero == 0 )
            return 0;   // non-magnetic
        else 
            return -1;  // AFM
    }
    else{
        return 2;       // Ferro- or Ferri- magentic
    }
    return 1;
}

void derive_symm_for_magnetic_materials(double rotations[][3][3],
                                        double translations[][3],
                                        int TR[],
                                        int * p2number_of_symm,
                                        double lattice[3][3],
                                        double atom_positions[][3],
                                        int *atom_types,
                                        int natom,
                                        vector *magmom,
                                        double eps)
{
    int i,j,k;
    int iatom, jatom, isymm;
    vector atpos[MAX_NUM_of_atoms]; //atoms' position in vector form
    vector atpos_symmed, mag_roted;
    vector origin;
    char msg[MAXLEN];

    int nsymm_crystal=*p2number_of_symm;
    int nsymm_mag=0;
    SymmetryOperator * magsymm;
    magsymm = (SymmetryOperator *) malloc(sizeof(SymmetryOperator)*nsymm_crystal);

    double rot_cart[3][3];
    double latt_trans[3][3];
    double latt_trans_inv[3][3];

    matrix3x3_transpose(latt_trans, lattice);
    matrix3x3_inverse(latt_trans_inv, latt_trans);

    init_vector(&origin, 0,0,0);

    int keep_symm_vote, tr_vote, tot_votes;

    for(iatom=0;iatom<natom;iatom++)
    {
        atpos[iatom]=array2vector(atom_positions[iatom]);
    }
    for(isymm=0;isymm<*p2number_of_symm;isymm++){
        keep_symm_vote = tr_vote = 0;
        tot_votes = natom;

        // rot_cart = trans(latt) * rot * inv(trans(latt))
        matrix3x3_dot(rot_cart, latt_trans, rotations[isymm]);
        matrix3x3_dot(rot_cart, rot_cart, latt_trans_inv);

        for(iatom=0;iatom<natom;iatom++)
        {
            //zero magmom have no right to vote
            if( distance(magmom[iatom], origin, NULL) < fabs(eps) ) {
                tot_votes--;
                continue;
            }
            atpos_symmed = vector_rotate(atpos[iatom], rotations[isymm]);
            atpos_symmed = vector_add(atpos_symmed, array2vector(translations[isymm]));
            mag_roted    = vector_rotate(magmom[iatom], rot_cart);
            if( matrix3x3_determinant(rotations[isymm]) < 0 ){
                //inversion do not have influence on spin space
                mag_roted = vector_scale(-1, mag_roted);
            }
            jatom        = find_roted_atom(atpos_symmed, atpos, natom, lattice);
            if(atom_types[iatom] != atom_types[jatom]){
                print_error("Error in finding the symmetries related to magnetism. \n");
                sprintf(msg, "The atom No.%3d rotated to atom No.%3d \n", iatom+1, jatom+1);
                print_error(msg);
                sprintf(msg, "atpos[%3d]=%16.9lf%16.9lf%16.9lf\n", iatom+1,atpos[iatom].x, atpos[iatom].y, atpos[iatom].z);
                print_error(msg);
                sprintf(msg, "atpos[%3d]=%16.9lf%16.9lf%16.9lf\n", jatom+1,atpos[jatom].x, atpos[jatom].y, atpos[jatom].z);
                print_error(msg);
                sprintf(msg, "atpos_symmed=%16.9lf%16.9lf%16.9lf\n", atpos_symmed.x, atpos_symmed.y, atpos_symmed.z);
                print_error(msg);
                sprintf(msg, "The related symmetry is No.%3d\n", isymm+1);
                print_error(msg);
                //exit(1);
            }
            if(  distance( magmom[jatom], mag_roted, NULL) < fabs(eps) )
            {
                keep_symm_vote++;
            }
            else if( distance( magmom[jatom], vector_scale(-1, mag_roted), NULL) < fabs(eps) )
            {
                tr_vote++;
            }
        }
        if( tot_votes == 0 && isymm == 0){
            sprintf(msg, "Warning: no magnetic ion found, please make sure the \"MAGMOM\" in input file is correctly set.");
            print_error(msg);
        }
        if( keep_symm_vote == tot_votes || tr_vote == tot_votes ){
            // this symmetry belongs to the magnetic group
            matrix3x3_copy(magsymm[nsymm_mag].rot, rotations[isymm]);
            matrix3x1_copy(magsymm[nsymm_mag].trans, translations[isymm]);
            if( keep_symm_vote == tot_votes ) {
                // symmetry without TRS
                magsymm[nsymm_mag].TR = 0;
            }
            else{
                // symmetry with TRS
                magsymm[nsymm_mag].TR = 1;
            }
            nsymm_mag++;
        }
    }
    for(isymm=0; isymm<nsymm_mag; isymm++){
        matrix3x3_copy(rotations[isymm], magsymm[isymm].rot);
        matrix3x1_copy(translations[isymm], magsymm[isymm].trans);
        TR[isymm] = magsymm[isymm].TR;
    }
    *p2number_of_symm = nsymm_mag;
}

int find_roted_atom(vector atom_roted, vector * atpos, int natom, double lattice[3][3])
{
    int i;
    vector dis, dis_int, dis_rem, dis_rem_cartesian;
    for(i=0;i<natom;i++)
    {
        dis = vector_sub(atom_roted, *(atpos+i));
        dis_int = vector_round(dis);
        dis_rem = vector_sub(dis, dis_int);
        dis_rem_cartesian = vector_rotate(dis_rem, lattice);
        if( vector_norm(dis_rem_cartesian) < 1E-3){
            return i;
        }
    }
    // if no corresponding atom found, we meet error
    print_error("ERROR in finding symmetry : can not find the rotated atom.");
    exit(1);
}

void remove_a_symm(double rotations[][3][3], double translations[][3], int TR[], int r)
{
    // remove the r th symmetry, r may be 0, 1, 2...
    int i,j,k;
    for(i=r;i<MAX_NUM_of_SYMM-1; i++)
    {
        for(j=0;j<3;j++){
            for(k=0;k<3;k++){
                rotations[i][j][k] = rotations[i+1][j][k];
            }
            translations[i][j] = translations[i+1][j];
        }
        TR[i] = TR[i+1];
    }
}

int file_exists(const char * fn){
    FILE * file;
    if ((file = fopen(fn, "r"))){
        fclose(file);
        return 1;
    }
    return 0;
}

void print_error(char * msg){
    FILE * ferr;
    ferr = fopen("wannsymm.err", "a");
    fprintf(ferr, "%s", msg);
    fprintf(ferr, "\n");
    fclose(ferr);
    print_msg("ERROR: see 'wannsymm.err' file for detail.\n");
}

void print_msg(char * msg){
    FILE * fstdout;
    fstdout = fopen("wannsymm.out", "a");
    fprintf(fstdout, "%s", msg);
    //fprintf(fstdout, "\n");   // "\n" is aborbed in the variable msg
    fclose(fstdout);
}

void print_flag(FILE * fout, char * tag, int arg){
    char msg[MAXLEN];
    char * str_true="True";
    char * str_false="False";

    if(arg == 1) {
        sprintf(msg, "%s = %s", tag, str_true);
    } else {
        sprintf(msg, "%s = %s", tag, str_false);
    }
    fprintf(fout, "%s\n", msg);
}

void print_symmetry(const char * fnout, double lattice[3][3], int isymm, double rot[3][3], double trans[3], int TR, int flag_showtrans, int flag_showmirror, int flag_dbgrp){
    // flag_dbgrp means flag_double_group_enabled
    // flag_showtrans 1: always show; 0: show if nonzero; -1: never show translation part of symmetry
    double rot_axis[3];
    double rot_angle;
    int inv_flag;
    int ii;
    FILE * fout;
    fout = fopen(fnout, "a");
    get_axis_angle_of_rotation(rot_axis, &rot_angle, &inv_flag, rot, lattice);

    //if (flag_dbgrp == 0 && rot_angle > PI + 1E-3){
    //    rot_angle -= 2*PI; // rot_angle range (-PI, PI]
    //}
    
    if(rot_angle < 0 && fabs(rot_angle) < 1E-3 ){
        rot_angle = 0;     // display -0.00 deg as 0.00 deg.
    }

    if(isymm < 0){
        //double group
        fprintf(fout, "symm%4d: ", isymm);
        fprintf(fout, "symm%4d with additional 360 deg rotation for spinor\n", abs(isymm));
    } else {
        fprintf(fout, "symm%4d:", isymm+1);
        if( fabs(rot_angle) < 1e-3 
            && TR != 1 
            && (trans == NULL || fabs(trans[0]) + fabs(trans[1]) + fabs(trans[2]) < eps4) )
        {
            if( inv_flag == 0 ){
                fprintf(fout, " Identity\n");
            } else {
                fprintf(fout, " Inversion\n");
            }
            fclose(fout);
            return;
        }
        fprintf(fout, "%6.1lf deg rot",rot_angle/PI*180);
        for(ii=0; ii<3; ii++){
            if(fabs(rot_axis[ii]) < eps8 )
                rot_axis[ii] = 0;  // display -0.0000 as 0.0000
        }
        fprintf(fout, " around (%7.4lf,%7.4lf,%7.4lf)",rot_axis[0], rot_axis[1], rot_axis[2]);
        if(trans != NULL && flag_showtrans != -1){
            for(ii=0; ii<3; ii++){
                if(trans[ii] < 0 && fabs(trans[ii]) < eps8 )
                    trans[ii] = 0;  // display -0.0000 as 0.0000
            }
            if( flag_showtrans == 1 || fabs(trans[0]) + fabs(trans[1]) + fabs(trans[2]) > eps4 ){
                fprintf(fout, " with trans (%7.4lf %7.4lf %7.4lf)",trans[0], trans[1], trans[2]);
            }
        }
        
        if(inv_flag){
            fprintf(fout," with inv");
            if(flag_showmirror==1 && ( fabs(rot_angle/PI*180-180) < 1.0 || fabs(rot_angle/PI*180+180) < 1.0 ))
                fprintf(fout," (Mirror)");
        }
        if(TR == 1)
            fprintf(fout, " with TRS");
        fprintf(fout,"\n");
    }    
    fclose(fout);
}

int str2boolean(char * arg){
    // allowed false: F* f* 0* .F* .f*  e.g. F  0  .False. .fALSE.
    // allowed true:  T* t* 1* .T* .t*  e.g. T  1  .True.  .TRUE.
    if( arg[0] == 'F' || arg[0]=='f' || arg[0]=='0' || (arg[0]=='.' && (arg[1]=='F' || arg[1]=='f')) )
        return 0;
    else if( arg[0] == 'T' || arg[0]=='t' || arg[0]=='1' || (arg[0]=='.' && (arg[1]=='T' || arg[1]=='t')) )
        return 1;
    else{
        return -1;
    }
}

void read_kpath_info(vec_llist ** p2kpaths, char klabels[][SHORTLEN], int * p2nkpath, FILE * fin){
    char line[MAXLEN];
    char tag[MAXLEN]="trivial initial string";
    char arg[MAXLEN]="trivial initial string";
    char * tmpstr;
    int i,j;
    int ipath;
    vector tmpkpt;
    double array_kpt[2][3];
    char * fgets_state;

    fgets_state = fgets(line, MAXLEN, fin);
    parseline(tag, arg, line, 1);
    for(ipath=0; ! (tag[0] == 'e' && tag[1] == 'n' && tag[2] == 'd') && fgets_state != NULL; ipath++){
        if(ipath > MAXLEN-1) continue;  // skip lines beyond MAXLEN to prevent memory buffer overflow.
        tmpstr = strtok(line, " ");
        for(i=0; i<8 && tmpstr!=NULL; i++){
            if(i%4 == 0){
                strcpy(klabels[2*ipath + i/4], tmpstr );
            } else {
                sscanf(tmpstr, "%lf", array_kpt[i/4]+i%4-1);
            }
            tmpstr = strtok(NULL, " ");
        }
        for(i=0; i<2; i++){
            tmpkpt = array2vector(array_kpt[i]);
            vec_llist_add(p2kpaths, tmpkpt);
        }
        fgets_state = fgets(line, MAXLEN, fin);
        parseline(tag, arg, line, 1);
    }
    *p2nkpath = ipath;
}

void read_kpts_info(vec_llist ** p2kpts, int * p2nkpt, FILE * fin){
    int ikpt;
    char line[MAXLEN];
    char tag[MAXLEN]="trivial initial string";
    char arg[MAXLEN]="trivial initial string";
    char * fgets_state;
    vector tmpkpt;
    
    fgets_state = fgets(line, MAXLEN, fin);
    parseline(tag, arg, line, 1);
    for(ikpt=0; ! (tag[0] == 'e' && tag[1] == 'n' && tag[2] == 'd') && fgets_state != NULL; ikpt++){
        sscanf(line, "%lf%lf%lf", &(tmpkpt.x), &(tmpkpt.y), &(tmpkpt.z));
        vec_llist_add(p2kpts, tmpkpt);
        (*p2nkpt)++;
        fgets_state = fgets(line, MAXLEN, fin);
        parseline(tag, arg, line, 1);
    }
}

