function [rc] = write(dim,VERT,SIMP);
%WRITE  Write out a 2D or 3D simplex mesh (coordinates always in R^3)
%
% Usage: [rc] = write(dim,VERT,SIMP);
%
% Input:
%
%    The input mesh is in internal MCLite format.
%
% Output:
%
%    The output mesh is in MCSF format.
%
%    (N = number of vertices in mesh)
%    (L = number of simplices in mesh)
%    (K = number of edges in mesh)
%
%    VERT = size(N,8) ==> vertices
%
%        VERT(:,1)  == vertex identifier
%        VERT(:,2)  == vertex processor id or color
%        VERT(:,3)  == vertex type
%        VERT(:,4)  == x-coordinates of all vertices
%        VERT(:,5)  == y-coordinates of all vertices
%        VERT(:,6)  == z-coordinates of all vertices (3d case & manifold)
%
%    SIMP = size(L,9) ==> 2-simplices
%
%        SIMP(:,1)  == simplex identifier
%        SIMP(:,2)  == simplex processor id or color
%        SIMP(:,3)  == simplex type
%        SIMP(:,4)  == 1st-face type; face opposite vertex 1
%        SIMP(:,5)  == 2nd-face type; face opposite vertex 2
%        SIMP(:,6)  == 3rd-face type; face opposite vertex 3
%        SIMP(:,7)  == 1st-vertex of the simplex
%        SIMP(:,8)  == 2nd-vertex of the simplex
%        SIMP(:,9)  == 3rd-vertex of the simplex
%
%    SIMP = size(L,11) ==> 3-simplices
%
%        SIMP(:,1)  == simplex identifier
%        SIMP(:,2)  == simplex processor id or color
%        SIMP(:,3)  == simplex type
%        SIMP(:,4)  == 1st-face type; face opposite vertex 1
%        SIMP(:,5)  == 2nd-face type; face opposite vertex 2
%        SIMP(:,6)  == 3rd-face type; face opposite vertex 3
%        SIMP(:,7)  == 4th-face type; face opposite vertex 4 (3d case)
%        SIMP(:,8)  == 1st-vertex of the simplex
%        SIMP(:,9)  == 2nd-vertex of the simplex
%        SIMP(:,10) == 3rd-vertex of the simplex
%        SIMP(:,11) == 4th-vertex of the simplex (3d case)
%
% Author:   Michael Holst
% rcsid="$Id: write.m,v 1.1.1.1 2007/04/27 08:28:07 hrg Exp $"

%%% change to external MCSF-style storage

   [N,dummy]=size(VERT);
   [L,dummy]=size(SIMP);
   ONEN=ones(N,1);
   ONEL=ones(L,1);
   VERT_F=zeros(N,5);
   assertw( ((dim==2) | (dim==3)), 'write 1' );
   if (dim==2)
       SIMP_F=zeros(L,9);
   else
       SIMP_F=zeros(L,11);
   end

   VERT_F(:,1:2) = VERT(:,4:5);
   VERT_F(:,3:5) = VERT(:,1:3);
   VERT_F(:,1)   = VERT_F(:,1) - ONEN;

   if (dim == 2)
       SIMP_F(:,1:3)  = SIMP(:,9:11);
       SIMP_F(:,4:6)  = SIMP(:,5:7);
       SIMP_F(:,7:9)  = SIMP(:,1:3);
       SIMP_F(:,1)    = SIMP_F(:,1) - ONEL;
       SIMP_F(:,7)    = SIMP_F(:,7) - ONEL;
       SIMP_F(:,8)    = SIMP_F(:,8) - ONEL;
       SIMP_F(:,9)    = SIMP_F(:,9) - ONEL;
   else
       SIMP_F(:,1:3)  = SIMP(:,9:11);
       SIMP_F(:,4:7)  = SIMP(:,5:8);
       SIMP_F(:,8:11) = SIMP(:,1:4);
       SIMP_F(:,1)    = SIMP_F(:,1)  - ONEL;
       SIMP_F(:,8)    = SIMP_F(:,8)  - ONEL;
       SIMP_F(:,9)    = SIMP_F(:,9)  - ONEL;
       SIMP_F(:,10)   = SIMP_F(:,10) - ONEL;
       SIMP_F(:,11)   = SIMP_F(:,11) - ONEL;
   end

   dimii = 3;
   writemcs(dim,dimii,VERT_F,SIMP_F);

