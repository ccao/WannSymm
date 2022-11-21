## Use WannSymm with Wien2K

There are several issues we need to pay special attention to when use WannSymm with Wien2K. 

1. The Wannier orbital definition order is different in Wien2K (more precisely, in Wien2Wannier) than in Wannier90. As WannSymm default to Wannier90 order, this must be changed. This part is done by modifying the write_inwf_lapw code in Wien2K distribution. Use your favorite editor to open the write_inwf_lapw file and find the line begin with _ang_func_aliases. In Wien2K 19.1, it is located at line 261, looks like:

    ```python
    _angfunc_aliases = {
        'p' : ['px', 'py', 'pz'],
        'p-x' : ['px'],
        'p-y' : ['py'],
        'p-z' : ['pz'],

    'd' : ['dxy', 'dxz', 'dyz', 'dx^2-y^2', 'dz^2'],
    ……
    ……
    ```

    Change the order of p, d, f aliases into:

    ```python
    ……
        'p' : ['pz', 'px', 'py'],
    ……
        'd' : ['dz^2', 'dxz', 'dyz', 'dx^2-y^2', 'dxy'],
    ……
        'f' : ['fz^3', 'fxz^2', 'fyz^2', 'fz(x^2-y^2)', 'fxyz', 'fx(x^2-3y^2)', 'fy(3x^2-y^2)'],
    ……
    ```

    So that the output order of Wien2Wannier would be the same as Wannier90 default.

2. The Wannier orbitals generated with Wien2Wannier has “local axis” by default. In addition, if we specify “local axis” in write_inwf, it rotates the spin as well (in the SOC case), which is not desired in Wannier90. Therefore, it is highly recommended NOT to use any local axis specification in write_inwf, and use the utility mkw2kinput.x provided to generate appropriate WannSymm input file with local-axis.
A typical symmetrize procedure goes as following:
    1. Normal Wien2K calculation. Keep the case.output0 file. It will be used by  mkw2kinput.x
    2. Use the Wien2Wannier code to generate the files required in Wannier90. Remember to change the write_inwf_lapw file before use it to generate correctly ordered Wannier orbitals. Keep the case.inwf(up) file, which will also be used by mkw2kinput.x
    3. Perform Wannierize using Wannier90. In order to keep maximum symmetry, do not minimize the spreading (num_iter = 0). Keep the seed.wout file, it will also be used by mkw2kinput.x
    Use mkw2kinput.x to generate input template file for WannSymm. The code will spit out input file template with appropriate specification of local axis, for example:

    ```bash
    [ccao@everest:mprj_local]$ ./mkw2kinput.x seed.wout case.output0 case.inwfup
    ######### WannSymm INPUT File Follows ########
    # Modify according to your system
    # template input file of wannsymm.x
    # anything following '#', '!' or '//' in a line will be regard as comments
    # tag names are case insensitive( SeedName and seednAme are equivalent)

    DFTcode  = VASP

    Spinors  = T

    SeedName ='wannier90'

    Use_POSCAR = 'POSCAR'

    # Projections_in_Format_of_wannier90

    begin projections
    f=   0.01004,  0.25000,  0.69361 : d : z=   0.00000,  1.00000,  0.00000 : x=   0.00000,  0.00000,  1.00000 : y=   1.00000,  0.00000,  0.00000
    f=  -0.01004, -0.25000, -0.69361 : d : z=   0.00000, -1.00000,  0.00000 : x=   0.00000,  0.00000, -1.00000 : y=  -1.00000,  0.00000,  0.00000
    f=   0.51004,  0.25000, -0.19361 : d : z=   0.00000, -1.00000,  0.00000 : x=   0.00000,  0.00000, -1.00000 : y=   1.00000,  0.00000,  0.00000
    f=   0.48996, -0.25000,  1.19361 : d : z=   0.00000, -1.00000,  0.00000 : x=   0.00000,  0.00000,  1.00000 : y=  -1.00000,  0.00000,  0.00000
    f=   0.35703, -0.25000,  0.55810 : d : z=   0.00000, -1.00000,  0.00000 : x=   0.00000,  0.00000,  1.00000 : y=  -1.00000,  0.00000,  0.00000
    f=   0.14297,  0.25000,  0.05810 : d : z=   0.00000,  1.00000,  0.00000 : x=   0.00000,  0.00000,  1.00000 : y=   1.00000,  0.00000,  0.00000
    f=  -0.14297, -0.25000, -0.05810 : d : z=   0.00000, -1.00000,  0.00000 : x=   0.00000,  0.00000, -1.00000 : y=  -1.00000,  0.00000,  0.00000
    f=   0.64297,  0.25000,  0.44190 : d : z=   0.00000, -1.00000,  0.00000 : x=   0.00000,  0.00000, -1.00000 : y=   1.00000,  0.00000,  0.00000
    end projections

    restart = F
    # Kpoint for calculating band's character of every symmetry
    #kpt = 0 0 0
    #chaeig_in_kpath = T

    nk_per_kpath = 101

    beginkpath
    endkpath
    ```

    4. Copy the screen output to the input file. Modify the input file according to your need. Most likely the only thing you need to change is between beginkpath and endkpath. But please check the projection specification. The orbital centers are calculated from Wannier90 output file, and are most likely different from the atomic positions. There are two reasons for the difference. Firstly, the atomic positions in Wien2K may differ from input by integer number of lattice vectors. Secondly, the projection itself may introduce small deviation from atomic positions. Therefore, one needs to check the orbital center positions has the correct symmetry. In addition, it is highly recommended to use these orbital center positions in the POSCAR file. 

    5. Symmetrize and enjoy.
