function [VERT_F,SIMP_F] = delsring(VERT,SIMP,which);
%DELSRING   Delete simplex from simplex rings associated with all its vertices
%
% Usage: [VERT_F,SIMP_F] = delsring(VERT,SIMP,which);
%
% Input:
%
%    which  = the simplex number to delete from rings
%    VERT   = <see "read.m" for description of the datastructure>
%    SIMP   = <see "read.m" for description of the datastructure>
%    
% Output:
%
%    VERT_F = list of vertices in the mesh with the new ring info
%    SIMP_F = list of simplices in the mesh with the new ring info
%
% Author:   Michael Holst
% rcsid="$Id: delsring.m,v 1.1.1.1 2007/04/27 08:28:07 hrg Exp $"

%%% recover array dimensions

   T=3;

%%% delete the "which" simplex from the rings

   for i=1:T
      vn    = SIMP( which, i    );  %%% the current vertex number
      vr    = SIMP( which, i+11 );  %%% the next simplex in the ring
      nextS = VERT( vn, 7 );        %%% first simp in ring of current vertex
      loop = 0;
      while ( nextS ~= 0 )
         loop=loop+1;
         if ( nextS == which )
            if (loop == 1)
                VERT( vn, 7 ) = vr;
            else
                SIMP( prevS, prevI+11 ) = vr;
            end
            SIMP(which,i+11) = 0;
            nextS = 0;
         else
            ii = -1;
            for j=1:T
               if ( SIMP(nextS,j) == vn ) ii = j; end
            end
            prevI = ii;
            prevS = nextS;
            nextS = SIMP( prevS, prevI+11 );
            assertw( (ii>0), 'delsring 1' );
         end
      end
   end

%%% return the results

   VERT_F = VERT;
   SIMP_F = SIMP;

