function [VERT_F,SIMP_F] = read(key);
%READ  Define a 2D or 3D simplicial finite element (coordinates always in R^3)
%
% Usage: [VERT,SIMP] = read(n);
%
% Input:
%
%    An "MCSF" (MC-Simplex-Format) mesh definition file.
%
% Output:
%
%    (N = number of vertices in mesh)
%    (L = number of simplices in mesh)
%    (K = number of edges in mesh)
%
%    VERT = size(N,8) ==> vertices
%
%              VERT(:,1)  == x-coordinates of all vertices
%              VERT(:,2)  == y-coordinates of all vertices
%              VERT(:,3)  == z-coordinates of all vertices (3d case & manifold)
%              VERT(:,4)  == vertex identifier
%              VERT(:,5)  == vertex processor id or color
%              VERT(:,6)  == vertex type
%              VERT(:,7)  == blank; later used as first simplex ptr
%              VERT(:,8)  == blank; later used as first edge ptr
%
%    SIMP = size(L,17) ==> simplices
%
%              SIMP(:,1)  == 1st-vertex of the simplex
%              SIMP(:,2)  == 2nd-vertex of the simplex
%              SIMP(:,3)  == 3rd-vertex of the simplex
%              SIMP(:,4)  == 4th-vertex of the simplex (3d case)
%              SIMP(:,5)  == 1st-face type; face opposite vertex 1
%              SIMP(:,6)  == 2nd-face type; face opposite vertex 2
%              SIMP(:,7)  == 3rd-face type; face opposite vertex 3
%              SIMP(:,8)  == 4th-face type; face opposite vertex 4 (3d case)
%              SIMP(:,9)  == simplex identifier
%              SIMP(:,10) == simplex processor id or color
%              SIMP(:,11) == simplex type
%              SIMP(:,12) == blank; later used as simplex ring ptr 1
%              SIMP(:,13) == blank; later used as simplex ring ptr 2
%              SIMP(:,14) == blank; later used as simplex ring ptr 3
%              SIMP(:,15) == blank; later used as simplex ring ptr 4 (3d case)
%              SIMP(:,16) == first queue marker
%              SIMP(:,17) == second queue marker
%              SIMP(:,18) == marked marker
%              SIMP(:,19) == creation type (0=bisection,1=quad,2=quad interior)
%              SIMP(:,20) == simplex generation number
%
%    EDGE = size(K,8) ==> edges
%
%              EDGE(:,1) == 1st vertex making up edge
%              EDGE(:,2) == 2nd vertex making up edge
%              EDGE(:,3) == midpoint vertex number
%              EDGE(:,4) == edge identifier
%              EDGE(:,5) == edge processor id or color
%              EDGE(:,6) == edge type
%              EDGE(:,7) == next edge number in ring of edges about vertex 1
%              EDGE(:,8) == next edge number in ring of edges about vertex 2
%
% Author:   Michael Holst
% rcsid="$Id: read.m,v 1.1.1.1 2007/04/27 08:28:07 hrg Exp $"

%%% read in the mcf file

    mcin

%%% this mapping allows us to grab the other two vertices given the third

    vmapOV3 = [ 2, 3;
                1, 3;
                1, 2 ];

%%% change to internal storage

    N=vertices;
    L=simplices;
    VERT_F=zeros(N,8);
    SIMP_F=zeros(L,20);
    ZERON=zeros(N,1);
    ZEROL=zeros(L,1);
    ONEN=ones(N,1);
    ONEL=ones(L,1);

    VERT_F(:,1:3) = vert(:,3:5);
    VERT_F(:,4:5) = vert(:,1:2);
    VERT_F(:,6)   = ZERON;
    VERT_F(:,4)   = VERT_F(:,4) + ONEN;

    SIMP_F(:,18)  = ONEL;
    SIMP_F(:,19)  = 2*ONEL;
    SIMP_F(:,20)  = ONEL;

    if (dim == 2)
        SIMP_F(:,1:3)  = simp(:,7:9);
        SIMP_F(:,5:7)  = simp(:,4:6);
        SIMP_F(:,9:11) = simp(:,1:3);
        SIMP_F(:,1)    = SIMP_F(:,1)  + ONEL;
        SIMP_F(:,2)    = SIMP_F(:,2)  + ONEL;
        SIMP_F(:,3)    = SIMP_F(:,3)  + ONEL;
        SIMP_F(:,9)    = SIMP_F(:,9)  + ONEL;
 
        %%% calculate the vertex types from the face types

        for element=1:L
            for vertex=1:3
                % vertex id and type (type is possibly incomplete)
                vid = SIMP_F(element,vertex);
                vtp = VERT_F(vid,6);

                % the two faces that use the vertex and their types
                f1   = vmapOV3(vertex,1);
                f2   = vmapOV3(vertex,2);
                f1tp = SIMP_F(element,4+f1);
                f2tp = SIMP_F(element,4+f2);

                % if anything touching vertex is diri, mark vertex as such
                if ( vdiri(vtp) | vdiri(f1tp) | vdiri(f2tp) )
                    tmp = max( vdiri(f1tp)*f1tp, vdiri(f2tp)*f2tp );
                    VERT_F(vid,6) = max( vdiri(vtp)*vtp, tmp );

                % else if anything touching vertex is neum, mark vertex as such
                else
                    tmp = max( vneum(f1tp)*f1tp, vneum(f2tp)*f2tp );
                    VERT_F(vid,6) = max( vneum(vtp)*vtp, tmp );
                end
            end
        end

    else
        SIMP_F(:,1:4)  = simp(:,8:11);
        SIMP_F(:,5:8)  = simp(:,4:7);
        SIMP_F(:,9:11) = simp(:,1:3);
        SIMP_F(:,1)    = SIMP_F(:,1)  + ONEL;
        SIMP_F(:,2)    = SIMP_F(:,2)  + ONEL;
        SIMP_F(:,3)    = SIMP_F(:,3)  + ONEL;
        SIMP_F(:,4)    = SIMP_F(:,4)  + ONEL;
        SIMP_F(:,9)    = SIMP_F(:,9)  + ONEL;

        %%% WARNING: we don't calculate vertex types in 3D case
        assertw( (0), 'read: refusing to calulate 3D vertex types' );

    end

