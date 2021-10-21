# WannSymm

This program was developed by Chao Cao and Guo-Xiang Zhi @ Zhejiang University. Please cite the following paper:

- G.-X. Zhi, C.-C. Xu, S.-Q. Wu, F.-L. Ning and Chao Cao, WannSymm: A symmetry analysis code for Wannier orbitals, Computer Physics Communications, 271 (2022) 108196, doi: https://doi.org/10.1016/j.cpc.2021.108196.

Bugs and comments can be submitted via GitHub or email to:

ccao@zju.edu.cn

or

guoxiang.zhi@gmail.com

## Installation

1. install spglib
        https://atztogo.github.io/spglib/install.html

2. A "make.sys" file must exist in the package root directory, there are two templates in the root directory, copy an appropriate one as e.g.:

```
        cp make.sys.static make.sys
```

3. update make.sys with the path of spglib
  - static link:
```
        SPGINCLUDES=${your spglib installation dir}/include     
        SPGLIBS=${your spglib installation dir}/lib/libsymspg.a
```
  - dynamic link:
```
        SPGINCLUDES=${your spglib installation dir}/include     
        SPGLIBS=${your spglib installation dir}/lib -lsymspg
```

4. update make.sys with the path of MKLROOT
```
        MKLROOT=${your MKLROOT's path}
```

5. compile the program
```
        cd src/
        make all
```

6. the executable file is "bin/wannsymm.x"


## Run the Program
```
        mpirun -np ${Num_of_process} wannsymm.x ${InputFile}
```

##     Notes
- Please use mpicc or mpiicc for compiling.
- If no inputfile specified, a default input file "wannsymm.in" must be provided.
- Please check examples to find the format of "wannsymm.in".
- The output file "${seedname}_symmed_hr.dat" is the symmetrized Hamtonian
- The output file "symmetries.dat" contains the symmetries found by spglib
- The output file "wannsymm.out" contains important information of the calculation.
- If you want to know the progress of thread N while calculating:
```bash
        tail -f .progress-of-threadN
```
- If spglib is dynamically linked, before executing wannsymm.x, you need run this first
```
        export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${your spglib installation dir}/lib
```



## An Example of Input File

This example is for anti-ferromagnetic MnF<sub>2</sub> with SOC considered.
In the working directory, in this case, there should be a structure file 'POSCAR' and a Hamiltonian file 'MnF2_hr.dat'
```
DFTcode  = VASP
Spinors  = T

SeedName ='MnF2'
# file 'MnF2_hr.dat' must exist in the working directory

Use_POSCAR = 'POSCAR'
# file 'POSCAR' must exist in the working directory

begin projections
  Mn: s;  d
   F: p
end projections

MAGMOM = 0 0 5 0 0 -5 12*0
```




## Input File Description

Anything following '#', '!' or '//' in a line will be regard as comments.

#### Necessary Tags and Blocks

##### DFTcode tag
 ```
 Tag Name:      DFTcode   
 Type:          String   
 Description:   Type of DFT code used prior to the Wannier90 calculation.
                Available options are:
                'VASP' for VASP and WIEN2k 
                'QE'   for Quantum ESPRESSO 
 ```

##### Spinors tag
 ```
 Tag Name:      Spinors
 Type:          Logical
 Description:   T or True  for spin-orbit coupling (SOC) calculations
                F or False for non-SOC calculations
 ```

##### SeedName tag
 ```
 Tag Name:      SeedName
 Type:          String   
 Description:   Specify the input and output real-space Hamiltonian 
                Input:  SeedName_hr.dat
                Output: SeedName_symmed_hr.dat
 ```

##### Use_POSCAR tag
 ```
 Tag Name:      Use_POSCAR
 Type:          String of filename
 Description:   File in POSCAR-format consisting of the crystal structure 
                information.
 ```


##### Projections Block
This block follows Wannier90 convention. 
An example is given here:
 ```
 begin projections 
 Mn: s;d
 F:  pz;px;py
 end projections
 ```

#### Optional Tags

##### MAGMOM tag
 ```
 Tag Name:      MAGMOM
                (optional)
 Type:          Array of numbers
 Description:   magnetic moment for each atom. 
                - for SOC considered calculation, Each magnetic moment 
                  consists of 3 numbers which are x, y and z components of 
                  the moment.  
                - for non-SOC calculation, Each magnetic moment consists 1 
                  number which is the magnitude of the moment.
                Note: 
                - non-magnetic calculations should avoid this tag.
                - This tag conflicts with Global_TRsymm tag.
 ```

##### symm_magnetic_tolerance tag
 ```
 Tag Name:      symm_magnetic_tolerance 
                (optional) 
 Type:          Number
 Default:       1E-3
 Description:   if the norm of difference of two magnetic moments is smaller
                than symm_magnetic_tolerance, the two moments are considered
                equal.
                Used for determined the magnetic group.
 ```

##### Global_TRsymm tag
 ```
 Tag Name:      Global_TRsymm
                (optional)
 Type:          Logical
 Default:       True  for non-magnetic calculation
                False for magnetic calculation
 Description:   Specify whether the global time-reversal symmetry is
                considered or not.
 ```

##### kpt tag
 ```
 Tag Name:      kpt
                (optional) 
 Type:          Array of numbers
 Description:   the x, y and z components of the k-point used for calculation
                of characters.
 ```

##### degenerate_tolerance tag
 ```
 Tag Name:      degenerate_tolerance
                (optional) 
 Type:          Number
 Default:       1E-6 
 Description:   degenerate tolerance (in eV) of two band used for calculation of
                characters.
 ```

##### restart tag
 ```
 Tag Name:      restart
                (optional)
 Type:          Logical
 Default:       False
 Description:   if True the program will read the symmetrized Hamiltonian
                directly and skip the symmetrization procedure. 
                Useful for calculation of characters.
 ```

##### use_symmetry tag
 ```
 Tag Name:      use_symmetry
                (optional)
 Type:          String of filename
 Description:   File consisting of information for symmetry operations.
                If not provided, WannSymm will find the symmetry operations 
                automatically.
                
 An example of this file reads:
nsymm=6
--- 1 ---
 1  0  0
 0  1  0
 0  0  1
0.000000 0.000000 0.000000 F
--- 2 ---
-1  1  0
-1  0  0
 0  0 -1
0.000000 0.000000 0.000000 F
--- 3 ---
...
...
1st line:       number of symmetries in this file.
2nd line:       symmetry operation number
3rd-5th lines:  rotation part of the symmetry.
6th line:       translation part of the symmetry  and a tag indicating 
                 time-reversal symmetry (TRS).
 ```

##### EverySymm tag
 ```
 Tag Name:      EverySymm
                (optional)
 Type:          Logical
 Default:       False
 Description:   if T or True, the program will output the Hamiltonians 
                corresponding to every symmetry operation to files:
                symm01_hr.dat
                symm02_hr.dat
                ...
 ```

##### ExpandRvec tag
 ```
 Tag Name:      ExpandRvec
                (optional)
 Type:          Logical
 Default:       True
 Description:   if F or False, when writing  the  symmetrized  Hamiltonian, 
                the program will ommit the R-vectors not  included  in  the 
                original Hamiltonian. So the symmetrized Hamiltonian   will
                have the same dimension as the original one. And the program
                will calculate the difference of symmetrized Hamiltonian and
                the original one. If the difference of some components of    
                the two Hamiltonians are larger than ham_tolerance, warnings  
                will appear in the wannsymm.out file.
 ```

##### ham_tolerance tag
 ```
 Tag Name:      ham_tolerance
                (optional)
 Type:          Number
 Default:       0.1
 Description:   Only make sense when ExpandRvec=False.
 ```

#### Experimental Tags and Blocks (are not fully tested)

##### SAXIS tag
 ```
 Tag Name:      SAXIS
                (optional, experimental, not fully implemented)
 Type:          Array of Numbers
 Default:       0 0 1
 Description:   quantisation axis for spin 
 ```

##### structure_in_format_of_POSCAR Block
Block that providing crystal structure information, conflicts with 'Use_POSCAR' tag, and can be used to substitute for 'Use_POSCAR' tag.

An example of this block reads:

 ```
begin structure_in_format_of_POSCAR
some comments
   1.00000000000000
     3.8380000000000000    0.0000000000000000    0.0000000000000000
     0.0000000000000000    3.8380000000000000    0.0000000000000000
     0.0000000000000000    0.0000000000000000    3.8380000000000000
   Sr   V    O
     1     1     3
Direct
  0.5000000000000000  0.5000000000000000  0.5000000000000000
  0.0000000000000000  0.0000000000000000  0.0000000000000000
  0.5000000000000000  0.0000000000000000  0.0000000000000000
  0.0000000000000000  0.5000000000000000  0.0000000000000000
  0.0000000000000000  0.0000000000000000  0.5000000000000000
end structure_in_format_of_POSCAR
 ```

##### Output_Mem_Usage Tag
```
 Tag Name:      Output_Mem_Usage
                (optional)
 Type:          Logical
 Default:       False
 Description:   if True output memory usage of every thread
```
