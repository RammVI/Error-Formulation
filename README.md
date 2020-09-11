# Error-Formulation

Refines grids, based on adjoint error formulations of the Solvation Energy.

## Requirements
  The following softwares are needed to run this code.
  * python2.7
  * bempp v2.0 
  * fetk package installed in the Software folder
  * MSMS and NanoShaper installed
## Usage
  
  * Example: 
  > python2 main.py -M=meth_test -D=2.0 -N_it=3 -E=E_u
  
  For more options run 
  > python2 main.py -h
  
  The code allows to change a lot of parameters of the process, which are listed on the help section.

## Output data

  The main.py rutine requires at least that the user has a .pdb or .pqr file in Molecule/{name}/{name}.pdb or Molecule/{name}/{name}.pqr, then two grids with {dens} and 40.0 elements are created. This last one is used for the smoothing option, and may be replaced in the future by a geometry extrapolation.

