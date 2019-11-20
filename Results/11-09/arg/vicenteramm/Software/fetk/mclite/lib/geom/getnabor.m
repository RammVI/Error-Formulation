function [nab,nface] = getnabor(VERT,SIMP,me,face);
%GETNABOR  Find the nabor simplex (if it exists) sharing a face
%
% Usage: [nab,nface] = getnabor(VERT,SIMP,me,face);
%
% Input:
%
%    face   = the particular face of the simplex we want the nabor of
%    me     = the particular simplex we are working with
%    VERT   = <see "read.m" for description of the datastructure>
%    SIMP   = <see "read.m" for description of the datastructure>
%    
% Output:
%
%    nab    = number of the nabor simplex (zero if not found)
%    nface  = (local) number of corresponding face in nabor (zero if no nabor)
%
% Author:   Michael Holst
% rcsid="$Id: getnabor.m,v 1.1.1.1 2007/04/27 08:28:07 hrg Exp $"

%%% recover array dimensions

   T=3;

%%% get the local vertex numbers forming this face

   if (face==1) v1=2; v2=3; end;
   if (face==2) v1=3; v2=1; end;
   if (face==3) v1=2; v2=1; end;

%%% fine the nabor simplex using the rings

   g1 = SIMP(me,v1);
   g2 = SIMP(me,v2);
   
   %%% get the simplex rings for the first vertex
   ring1 = [];
   nextS = VERT(g1,7);
   if ((nextS ~= me) & (nextS ~= 0))
      ring1 = [ ring1 ; nextS ];
   end
   while ( nextS ~= 0 )
      locv = SIMP(nextS,1:3);
      lfind = 0;
      if (locv(1)==g1) lfind = 1; end;
      if (locv(2)==g1) lfind = 2; end;
      if (locv(3)==g1) lfind = 3; end;
      assertw( (lfind>0), 'getnabor 1' );
      nextS = SIMP(nextS,lfind+11);
      if ((nextS ~= me) & (nextS ~= 0))
         ring1 = [ ring1 ; nextS ];
      end
   end

   %%% get the simplex rings for the second vertex
   ring2 = [];
   nextS = VERT(g2,7);
   if ((nextS ~= me) & (nextS ~= 0))
      ring2 = [ ring2 ; nextS ];
   end
   while ( nextS ~= 0 )
      locv = SIMP(nextS,1:3);
      lfind = 0;
      if (locv(1)==g2) lfind = 1; end;
      if (locv(2)==g2) lfind = 2; end;
      if (locv(3)==g2) lfind = 3; end;
      assertw( (lfind>0), 'getnabor 2' );
      nextS = SIMP(nextS,lfind+11);
      if ((nextS ~= me) & (nextS ~= 0))
         ring2 = [ ring2 ; nextS ];
      end
   end

   %%% see if there is a common simplex
   [n1,one] = size(ring1);
   [n2,one] = size(ring2);
   nb = 0;
   for i1=1:n1
      for i2=1:n2
         if (ring1(i1) == ring2(i2))
            assertw( (nb==0), 'getnabor 3' );
            nb = ring1(i1);
         end
      end
   end

   %%% if there is a common simplex, determine local face number
   nf = 0;
   if (nb ~= 0)
       h1 = SIMP(nb,1);
       h2 = SIMP(nb,2);
       h3 = SIMP(nb,3);
       if     ( ((g1 == h1) & (g2 == h2)) | ((g1 == h2) & (g2 == h1)) )
           nf = 3;
       elseif ( ((g1 == h1) & (g2 == h3)) | ((g1 == h3) & (g2 == h1)) )
           nf = 2;
       elseif ( ((g1 == h2) & (g2 == h3)) | ((g1 == h3) & (g2 == h2)) )
           nf = 1;
       end;
   end;

%%% return the results

   nab   = nb;
   nface = nf;

