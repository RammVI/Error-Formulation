function [dVec] = dirifind (VERT) ;
%FINDDIRI  Finds the dirichlet edges
%
% Usage: [dVec]=dirifind(VERT);
%
% Input:
%
%    VERT   = <see "read.m" for description of the datastructure>
%
% Output:
%
%    dVect  = stores the dirichlet edge index
%
% Author:   Michael Holst
% rcsid="$Id: dirifind.m,v 1.1.1.1 2007/04/27 08:28:07 hrg Exp $"

   [N,eight] = size(VERT);     
   dIndex = 0;
   dVec = [];
   for i=1:N   
      if (VERT(i,6)==1) 
         dIndex = dIndex + 1;
         dVec(dIndex) = i;
      end;
   end;
