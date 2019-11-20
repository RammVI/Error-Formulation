% Description and Overview of MCLite = < MC-Lite >
%
% Copyright (C) 1994  Michael Holst
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
% See the GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program; if not, write to the Free Software Foundation, Inc.,
% 675 Mass Ave, Cambridge, MA 02139, USA.
%
% The MCLite package has two primary functions:
%
%   1. MCLite implements piecewise linear simplex-based finite elements in
%      the plane, a posteriori error estimation, adaptive mesh refinement
%      (quadrasection and bisection), continuation methods, and global
%      inexact-newton methods for nonlinear elliptic equations on arbitrary
%      polyhedral domains in R^2.  It is essentially a 2D MATLAB prototype
%      of a more involved 2D/3D ANSI-C code called "MC".  Both MC and MCLite
%      were developed by Michael Holst at Caltech and UC San Diego.  They
%      both share the same core geometry datastructure referred to as the
%      RIVER (RInged VERtex) structure, and the same geometry algorithms.
%      Development of a new feature in MC is usually preceeded by a
%      prototype implementation in MCLite, so that the two packages usually
%      remain close in structure, functionality, as well as user-interface.
%
%   2. MCLite can be used as an interactive tool to create, modify, and
%      examine meshes of 2-simplices and 3-simplices in "MCSF" format
%      (MC-Simplex-Format) for later use by the package MC (the MCSF file
%      format is the same mesh file format used by the code MC).  MCSF file
%      format is actually legal MATLAB script syntax.  The compatibility
%      between MC and MCLite allows for the use of the MC utilities "MCsg"
%      and "MCbridge" to replace MATLAB's builtin graphics in favor of
%      Geomview or MCsg (a minimal Geomview clone).
%
%      Doug Arnold has also recently written a supporting library for MCLite
%      specifically for examining the quality of 3-simplex (tetrahedral)
%      meshes; his library is included with MCLite (see below).
%
% If you find MCLite useful and would like to contribute to its development,
% you retain the copyrights to your contributions through the GNU license;
% see the COPYING file described below.  The latest version of MCLite can
% always be found at the following website:
%
%      http://www.scicomp.ucsd.edu/~mholst/
%
% The MCLite package can be used to solve the following general second order
% nonlinear elliptic equation posed in a weak formulation:
%
%      Find [u-g] in H^1_{0,D}(M) s.t.  <F(u),v>=0, for all v in H^1_{0,D}(M),
%
% where <F(u),v> is the nonlinear weak form generated (by appropriate
% integration by parts) from a second-order, divergence-form, nonlinear 
% elliptic equation in a polyhedral domain in R^2.  Given the general nature
% of the above formulation, quite general boundary conditions can be employed.
% The user must specify the Dirichlet function g, the domain M, the nonlinear
% form <F(u),v>, and also a linearization form <DF(u),v>.
%
% The libraries which form MCLite are as follows:
%
%      ADAPT         Refinement and error estimation library
%
%      GEOM          Geometry (RIVER=RInged VERtex datastructure) library
%
%      IO            I/O library
%
%      MASTER        Master element data generation library
%
%      MCSF          The simplex mesh library (MCSF=MC Simplex Format)
%
%      PDE           PDE library
%
%      SOLV          Finite element assembly and solver library
%
%      UTIL          Misc utilities
%
% A sample working directory called "work" is also provided, initially
% containing the following support files for the complete definition of
% a test problem:
%
%      go.m          A sample MCLite driver program
%      gomin.m       A sample MCLite driver program
%      goo.m         A sample MCLite driver program
%      mclitego.m    The main MCLite initialization routine
%
%      fu.m          A nonlinear strong form defining the pde
%      fu_v.m        A nonlinear weak form defining the pde
%      dfu_wv.m      A linearization bilinear form for the pde
%      rho.m         A mass function for building mass matrices
%      u_d.m         A dirichlet boundary function
%      u_t.m         An optional analytical solution for testing purposes
%
%      mcin.m        A sample domain mesh file in MCSF format
%
%      edgsplit.m    A user-defined procedure to bisect edges in the mesh
%
%      startsg       A Bourne shell script for starting socket graphics
%      getp          A Bourne shell script for getting PDE files
%      putp          A Bourne shell script for putting PDE files
%      difp          A Bourne shell script for diffing PDE files
%
% Author:   Michael Holst and contributing authors
% rcsid="$Id: Contents.m,v 1.1.1.1 2007/04/27 08:28:05 hrg Exp $"
