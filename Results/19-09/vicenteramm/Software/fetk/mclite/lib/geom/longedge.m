function [canon] = longedge(VERT,SIMP,which);
%LONGEDGE  Determine longest in simplex; unique integer label breaks ties
%
% Usage: [canon] = longedge(VERT,SIMP,which);
%
% Input:
%
%    which  = the simplex number to determine the longest edge of
%    VERT   = <see "read.m" for description of the datastructure>
%    SIMP   = <see "read.m" for description of the datastructure>
%
% Output:
%
%    canon = the edge number of the longest edge
%
% A useful mapping:
%
%    mapE[3][3] ==> Convenient (invertible, pairs-to-singles) edge mapping:
%
%                   Edge (vtx pair)  Edge (number)
%                   --------------   -----------
%                   (1,2) or (2,1)   1
%                   (1,3) or (3,1)   2
%                   (2,3) or (3,2)   3
%
% Author:   Michael Holst
% rcsid="$Id: longedge.m,v 1.1.1.1 2007/04/27 08:28:07 hrg Exp $"

%%% this mapping allows us to map vertex pairs onto single edge numbers

   mapE  = [ 0,  1,  2;
             1,  0,  3;
             2,  3,  0 ];

%%% recover array dimensions

   T=3;

%%% get the element info

   xname  = zeros(T);    %%% vertex numbers of this simp
   xcoord = zeros(T,3);  %%% the x,y,z coordinates of the vertices
   for i=1:T
      xname(i)    = SIMP( which, i );
      xcoord(i,:) = VERT( xname(i), 1:3 );
   end

%%% determin the longest edge

   maxI   = -1;
   maxJ   = -1;
   maxLen = -1;
   maxId  = -1;
   for i=1:T
      xi = xcoord(i,:);
      for j=i+1:T
         xj = xcoord(j,:);

         %%% first get unique global id for tie-breaking
         myId = 1000000 * xname(i) + xname(j);

         len = norm(xi-xj);
         if ( len >= maxLen ) 
             if ( (len > maxLen) | (myId < maxId) )
                maxLen = len;
                maxId  = myId;
                maxI = i;
                maxJ = j;
             end
         end
      end
   end
   assertw( (maxI >= 0), 'longedge 1' );
   assertw( (maxJ >= 0), 'longedge 2' );

%%% return result

   canon = mapE(maxI,maxJ);

