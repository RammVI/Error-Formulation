function [rc] = writemce(dim,dimii,VERT,EDGE);
%WRITEMCE  Write out a 2D or 3D edge mesh (coordinates always in R^3)
%
% Usage: [rc] = writemce(dim,dimii,VERT,EDGE);
%
% Input:
%
%    The input mesh is in MCEF format (NOT MCLite format).
%
% Output:
%
%    The output mesh is in MCEF format.
%
%    (N = number of vertices in mesh)
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
%    EDGE = size(K,9) ==> edges
%
%        EDGE(:,1)  == edge identifier
%        EDGE(:,2)  == edge group number
%        EDGE(:,3)  == edge material type
%        EDGE(:,4)  == 1st-vertex of the edge
%        EDGE(:,5)  == 2nd-vertex of the edge
%
% Author:   Michael Holst
% rcsid="$Id: writemce.m,v 1.1.1.1 2007/04/27 08:28:07 hrg Exp $"

%%% change to external MCEF-style storage

   fid = fopen('mcout.m','W');

   [N,dummy]=size(VERT);
   [K,dummy]=size(EDGE);

fprintf(fid,'mcef_begin=1;\n');
fprintf(fid,'\n');
fprintf(fid,'      dim=%d;         %% intrinsic manifold dimension\n',dim);
fprintf(fid,'    dimii=%d;         %% imbedding manifold dimension\n',dimii);
fprintf(fid,' vertices=%d;         %% number of vertices\n',N);
fprintf(fid,'    edges=%d;         %% number of edges\n',K);
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
fprintf(fid,'edge=[\n');
fprintf(fid, '%%-------- ---- ---- ---------------------------------------\n');
fprintf(fid, ...
  '%% Edge-ID Grp  Mat  Vertex-Numbers\n');
fprintf(fid, '%%-------- ---- ---- ---------------------------------------\n');
for i=1:K
    fprintf(fid,'%-9d %-4d %-4d %d %d\n', ...
        EDGE(i,1),EDGE(i,2),EDGE(i,3),EDGE(i,4), EDGE(i,5) );
end
fprintf(fid,'];\n');
fprintf(fid,'\n');
fprintf(fid,'mcef_end=1;\n');

fclose(fid);

