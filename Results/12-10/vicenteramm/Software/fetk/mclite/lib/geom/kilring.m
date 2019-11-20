function [VERT_F,SIMP_F] = kilring(VERT,SIMP);
%KILRING  Destroys all simplex rings associated with all vertices in the mesh
%
% Usage: [VERT_F,SIMP_F] = kilring(VERT,SIMP);
%
% Input:
%
%    VERT   = <see "read.m" for description of the datastructure>
%    SIMP   = <see "read.m" for description of the datastructure>
%
% Output:
%
%    VERT_F = list of vertices in the mesh with the new ring info
%    SIMP_F = list of simplices in the mesh with the new ring info
%
% Author:   Michael Holst
% rcsid="$Id: kilring.m,v 1.1.1.1 2007/04/27 08:28:07 hrg Exp $"

%%% recover array dimensions

   [N,eight]     = size(VERT);
   [L,seventeen] = size(SIMP);

%%% remove all simplex ring info

   VERT(:,7:8)   = zeros(N,2);
   SIMP(:,12:15) = zeros(L,4);

%%% return the results

   VERT_F = VERT;
   SIMP_F = SIMP;

