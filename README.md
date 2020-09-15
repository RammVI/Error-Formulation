# Error-Formulation

Refines grids obtained from MSMS or NanoShaper, based on adjoint error formulations of the Solvation Energy.

## Requirements
  The following softwares are needed to run this code.
  * python2.7
  * [BEMPP](https://bempp.com/) v2.0 
  * [FETK](http://fetk.org/codes/gamer/) package installed in the Software folder
  * MSMS and NanoShaper installed
## Usage
  
  * Example: 
  ```
  python2 main.py -M=meth_test -D=2.0 -N_it=3 -E=E_u
  ```
  
  For more options run 
  ```
  python2 main.py -h
  ```
  
  The code allows to change a lot of parameters of the process, which are listed on the help section.

## Output data

  The main.py rutine requires that at least the user has a .pdb or .pqr file in Molecule/{name}/{name}.pdb or Molecule/{name}/{name}.pqr, then two grids with {dens} and 40.0 elements are created. This last one is used for the smoothing option, and may be replaced in the future by a geometry extrapolation.

  This process will generate a sequence of grids in .msh format, and a sequence of local error in .vtk format. Also will save the total electrostatic potential, in both spaces, the marked elements to be refined.
  
  The local error will be saved in a file named as \Molecule\{mol_name}\{mol_name}\_{dens}-{N_it}\_N_ref-{N_ref}.vtk
  
  ### Warnings
  Avoid using -N_ref grather than 4. RAM memory usage will be arround x4 times higher each time and the process can crash the computer, also, the number of iterations for solving phi when increasing this factor raises significally.
