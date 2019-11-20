function [rc] = writemcs(dim,dimii,VERT,SIMP);
%WRITEMCS  Write out a 2D or 3D simplex mesh (coordinates always in R^3)
%
% Usage: [rc] = writemcs(dim,dimii,VERT,SIMP);
%
% Input:
%
%    The input mesh is in MC format (NOT MCLite format).
%
% Output:
%
%    The output mesh is in MC format.
%
%    (N = number of vertices in mesh)
%    (L = number of simplices in mesh)
%    (K = number of edges in mesh)
%
%    VERT = size(N,8) ==> vertices
%
%        VERT(:,1)  == vertex identifier
%        VERT(:,2)  == vertex chart number
%        VERT(:,3)  == x-coordinates of all vertices
%        VERT(:,4)  == y-coordinates of all vertices
%        VERT(:,5)  == z-coordinates of all vertices (3d case & manifold)
%
%    SIMP = size(L,9) ==> 2-simplices
%
%        SIMP(:,1)  == simplex identifier
%        SIMP(:,2)  == simplex group number
%        SIMP(:,3)  == simplex material type
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
%        SIMP(:,2)  == simplex group number
%        SIMP(:,3)  == simplex material type
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
% rcsid="$Id: writemcs.m,v 1.1.1.1 2007/04/27 08:28:08 hrg Exp $"

%%% change to external MCSF-style storage

   fid = fopen('mcout.m','W');

   [N,dummy]=size(VERT);
   [L,dummy]=size(SIMP);

fprintf(fid,'mcsf_begin=1;\n');
fprintf(fid,'\n');
fprintf(fid,'      dim=%d;         %% intrinsic manifold dimension\n',dim);
fprintf(fid,'    dimii=%d;         %% imbedding manifold dimension\n',dimii);
fprintf(fid,' vertices=%d;         %% number of vertices\n',N);
fprintf(fid,'simplices=%d;         %% number of simplices\n',L);
fprintf(fid,'\n');
fprintf(fid,'vert=[\n');
fprintf(fid, '%%-------- ---- ----------------- ----------------- -----------------\n');
fprintf(fid, ...
    '%% Vert-ID Chrt X-Coordinate      Y-Coordinate      Z-Coordinate\n');
fprintf(fid, '%%-------- ---- ----------------- ----------------- -----------------\n');
for i=1:N
    fprintf(fid,'%-9d %-4d %17.10e %17.10e %17.10e\n', ...
    VERT(i,1), VERT(i,2), VERT(i,3), VERT(i,4), VERT(i,5) );
end
fprintf(fid,'];\n');
fprintf(fid,'\n');
fprintf(fid,'simp=[\n');
fprintf(fid, '%%-------- ---- ---- ------------------- ---------------------------------------\n');
fprintf(fid, ...
  '%% Simp-ID Grp  Mat  Face-Types          Vertex-Numbers\n');
fprintf(fid, '%%-------- ---- ---- ------------------- ---------------------------------------\n');
for i=1:L
    if (dim==2)
        fprintf(fid,'%-9d %-4d %-4d %-4d %-4d %-4d      %d %d %d\n', ...
            SIMP(i,1),SIMP(i,2),SIMP(i,3),SIMP(i,4), ...
            SIMP(i,5),SIMP(i,6),SIMP(i,7),SIMP(i,8),SIMP(i,9) );
    else
        fprintf(fid,'%-9d %-4d %-4d %-4d %-4d %-4d %-4d %d %d %d %d\n', ...
            SIMP(i,1),SIMP(i,2),SIMP(i,3),SIMP(i,4), ...
            SIMP(i,5),SIMP(i,6),SIMP(i,7), ...
            SIMP(i,8),SIMP(i,9),SIMP(i,10),SIMP(i,11) );
    end
end
fprintf(fid,'];\n');
fprintf(fid,'\n');
fprintf(fid,'mcsf_end=1;\n');

fclose(fid);

