function [VERT_F,SIMP_F] = bldsring(VERT,SIMP);
%BLDSRING  Builds all simplex rings around all vertices in the mesh
%
% Usage: [VERT_F,SIMP_F] = bldsring(VERT,SIMP);
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
% rcsid="$Id: bldsring.m,v 1.1.1.1 2007/04/27 08:28:07 hrg Exp $"

%%% recover array dimensions

   [N,eight]     = size(VERT);
   [L,seventeen] = size(SIMP);

%%% add each simplex one at a time to the associated rings

   [VERT,SIMP] = kilring(VERT,SIMP);
   for which=1:L
       [VERT,SIMP] = addsring(VERT,SIMP,which);
   end

%%% return the results

   VERT_F = VERT;
   SIMP_F = SIMP;

