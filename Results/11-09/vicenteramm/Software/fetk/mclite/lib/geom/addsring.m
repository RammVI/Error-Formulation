function [VERT_F,SIMP_F] = addsring(VERT,SIMP,which);
%ADDSRING   Add a simplex to the simplex rings associated with each vertex
%
% Usage: [VERT_F,SIMP_F] = addsring(VERT,SIMP,which);
%
% Input:
%
%    which  = the simplex number to add to rings
%    VERT   = <see "read.m" for description of the datastructure>
%    SIMP   = <see "read.m" for description of the datastructure>
%    
% Output:
%
%    VERT_F = list of vertices in the mesh with the new ring info
%    SIMP_F = list of simplices in the mesh with the new ring info
%
% Author:   Michael Holst
% rcsid="$Id: addsring.m,v 1.1.1.1 2007/04/27 08:28:06 hrg Exp $"

%%% recover array dimensions

   T=3;

%%% add the "which" simplex to the rings

   for i=1:T
      oldFirstS = VERT( SIMP(which,i) , 7 );
      VERT( SIMP(which,i) , 7 ) = which;
      SIMP(which,i+11) = oldFirstS;
   end

%%% return the results

   VERT_F = VERT;
   SIMP_F = SIMP;

