% RInged VERtex geometry datastructure library
%
% This library implements the RInged VERtex geometry datastructure,
% along with a number of traversal algorithms and datastructure
% manipulation procedures.
%
% The available routines are as follows:
%
%      ADDERING      add a single edge to all of its edge rings
%      ADDSRING      add a single simplex to all of its simplex rings
%      BLDSRING      build simplex rings for all simplices in the mesh
%      CHKMESH       extensive consistency/conformity check of mesh
%      DELSRING      delete a single simplex from all of its rings
%      DIRIFIND      find the dirichlet edges
%      GETEDGE       get the unique edge (if exists) between two vertices
%      GETNABOR      get the unique simplex (if exists) sharing a face
%      KILRING       destroy simplex rings for all simplices in the mesh
%      LONGEDGE      determine the longest edge of a single simplex
%      SINRING       is a particular simplex in a particular simplex ring
%      VBND          yes/no answer whether tri/edge/vert is boundary
%      VDIRI         yes/no answer whether tri/edge/vert is dirichlet
%      VINT          yes/no answer whether tri/edge/vert is interior
%      VNEUM         yes/no answer whether tri/edge/vert is neumann
%
% Author:   Michael Holst
% rcsid="$Id: Contents.m,v 1.1.1.1 2007/04/27 08:28:06 hrg Exp $"
