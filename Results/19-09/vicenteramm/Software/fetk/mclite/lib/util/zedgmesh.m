%ZEDGMESH  build [vert,edge] input from vertices and an adjacency matrix
%
% Usage: zedgmesh;
%
%    A quick hack to build [vert,edge] input for writemce.m,
%    from an adjacency matrix and a set of vertex coordinate.
%
%    We assume that the adjacency matrix is "A", and coordinates
%    in R2 or R3 is "xy", and that they are already loaded in matlab.
%
% Author:   Michael Holst
% rcsid="$Id: zedgmesh.m,v 1.1.1.1 2007/04/27 08:28:19 hrg Exp $"

% first grab the strict lower triangle of adjacency matrix for edge repres
B=tril(A,-1);

% now get the vertex pairs for each edge
[II,JJ,AA]=find(B);

% sizes
[n,m]=size(B);
[nedg,one]=size(AA);

% build the MCEF-compatible vertex array
vert=zeros(n,5);
vert(:,1)=[ 1:n ]' - ones(n,1);
vert(:,2)=zeros(n,1);
vert(:,3)=xy(:,1);
vert(:,4)=xy(:,2);
% vert(:,5)=xy(:,3);
vert(:,5)=zeros(n,1);

% build the MCEF-compatible edge array
edge=zeros(nedg,5);
edge(:,1)=[ 1:nedg ]' - ones(nedg,1);
edge(:,2)=zeros(nedg,1);
edge(:,3)=zeros(nedg,1);
edge(:,4)=II-ones(nedg,1);
edge(:,5)=JJ-ones(nedg,1);

% write the MCEF file
writemcef(2,3,vert,edge);
% writemcef(3,3,vert,edge);

