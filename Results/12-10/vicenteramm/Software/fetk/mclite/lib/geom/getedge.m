function [edg] = getedge(VERT,EDGE,g1,g2);
%GETEDGE  Find the edge (if it exists) connecting g1-g2
%
% Usage: [edg] = getedge(VERT,EDGE,g1,g2);
%
% Input:
%
%    g1     = first vertex of the edge
%    g2     = second vertex of the edge
%    VERT   = <see "read.m" for description of the datastructure>
%    EDGE   = <see "read.m" for description of the datastructure>
%    
% Output:
%
%    edg    = the number of the edge (or zero if not found)
%
% Author:   Michael Holst
% rcsid="$Id: getedge.m,v 1.1.1.1 2007/04/27 08:28:07 hrg Exp $"

%%% recover array dimensions

   T=3;

%%% fine the edge using the rings

   %%% get the edge rings for the first vertex
   ring1 = [];
   nextE = VERT(g1,8);
   if (nextE ~= 0)
      ring1 = [ ring1 ; nextE ];
   end
   while ( nextE ~= 0 )
      locv = EDGE(nextE,1:2);
      lfind = 0;
      if (locv(1)==g1) lfind = 1; end;
      if (locv(2)==g1) lfind = 2; end;
      assertw( (lfind>0), 'getedge 1' );
      nextE = EDGE(nextE,lfind+6);
      if (nextE ~= 0)
         ring1 = [ ring1 ; nextE ];
      end
   end

   %%% get the edge rings for the first vertex
   ring2 = [];
   nextE = VERT(g2,8);
   if (nextE ~= 0)
      ring2 = [ ring2 ; nextE ];
   end
   while ( nextE ~= 0 )
      locv = EDGE(nextE,1:2);
      lfind = 0;
      if (locv(1)==g2) lfind = 1; end;
      if (locv(2)==g2) lfind = 2; end;
      assertw( (lfind>0), 'getedge 2' );
      nextE = EDGE(nextE,lfind+6);
      if (nextE ~= 0)
         ring2 = [ ring2 ; nextE ];
      end
   end

   %%% see if there is a common edge
   [n1,one] = size(ring1);
   [n2,one] = size(ring2);
   eg = 0;
   for i1=1:n1
      for i2=1:n2
         if (ring1(i1) == ring2(i2))
            assertw( (eg==0), 'getedge 3' );
            eg = ring1(i1);
         end
      end
   end

%%% return the results

   edg = eg;

