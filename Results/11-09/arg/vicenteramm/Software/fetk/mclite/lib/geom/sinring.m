function [value] = sinring(VERT,SIMP,me,gv);
%SINRING  Check to see if a simplex is already on a simplex ring
%
% Usage: [value] = sinring(VERT,SIMP,me,gv);
%
% Input:
%
%    gv     = vertex holding the ring
%    me     = the particular simplex we are working with
%    VERT   = <see "read.m" for description of the datastructure>
%    SIMP   = <see "read.m" for description of the datastructure>
%    
% Output:
%
%    value  = 1 if in ring, 0 if not in ring
%
% Author:   Michael Holst
% rcsid="$Id: sinring.m,v 1.1.1.1 2007/04/27 08:28:07 hrg Exp $"

%%% recover array dimensions

   T=3;

%%% fine the nabor simplex using the rings

   %%% go through the simplex ring
   found = 0;
   nextS = VERT(gv,7);
   if (nextS == me) 
      found = 1;
   end
   while ((nextS ~= 0) & (~found))
      locv = SIMP(nextS,1:3);
      lfind = 0;
      if (locv(1)==gv) lfind = 1; end;
      if (locv(2)==gv) lfind = 2; end;
      if (locv(3)==gv) lfind = 3; end;
      assertw( (lfind>0), 'sinring 1' );
      nextS = SIMP(nextS,lfind+11);
      if (nextS == me) 
         found=1;
      end
   end

%%% return the result

   value = found;

