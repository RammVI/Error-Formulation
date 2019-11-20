function [value] = chkmesh(VERT,SIMP);
%CHKMESH  Check the mesh for correctness, including conformity
%
% Usage:  [value] = chkmesh(VERT,SIMP);
%
% Input:
%
%    VERT   = <see "read.m" for description of the datastructure>
%    SIMP   = <see "read.m" for description of the datastructure>
% 
% Output:
%    
%    value  = 1 if no error, 0 otherwise
%
% Author:   Michael Holst
% rcsid="$Id: chkmesh.m,v 1.1.1.1 2007/04/27 08:28:07 hrg Exp $"

   fprintf('chkmesh: starting conformity and consistency check of mesh.\n');

%%% recover array dimensions

   [N,eight]     = size(VERT);
   [L,seventeen] = size(SIMP);
   T=3;

%%% create the structures

   [VERT,SIMP] = bldsring(VERT,SIMP);

%%% check the structures

   %%% check the vertex structures
   fprintf('chkmesh: checking vertices..');
   t0 = clock;

   vertBndD = 0;
   vertBndN = 0;
   vertInt = 0;
   for i=1:N

      %%% check id for consistency
      assertw( (VERT(i,4)==i), 'chkmesh 1' );

      %%% every vertex should have at least one simplex in ring
      assertw( (VERT(i,7)~=0), 'chkmesh 2' );

      %%% count up the interior points
      if ( vint(VERT(i,6)) )
         vertInt = vertInt+1;
      end

      %%% count up the dirichlet points
      if ( vdiri(VERT(i,6)) )
         vertBndD = vertBndD+1;
      end

      %%% count up the neumann points
      if ( vneum(VERT(i,6)) )
         vertBndN = vertBndN+1;
      end
   end
   the_time = etime(clock,t0);
   fprintf('..done.  [time=%g].\n', the_time);

   %%% check the simplex structures
   fprintf('chkmesh: checking simplices..');
   t0 = clock;

   bndFaceD = 0;
   bndFaceN = 0;
   for i=1:L

      %%% check id for consistency
      assertw( (SIMP(i,9)==i), 'chkmesh 3' );

      %%% no simplex should be marked for refinement
      assertw( (SIMP(i,16)==0), 'chkmesh 4' );
      assertw( (SIMP(i,17)==0), 'chkmesh 5' );

      %%% facetypes
      nvert = SIMP(i,1:4);
      ftype = SIMP(i,5:8);

      %%% check vertex data and simplex rings for consistency
      for j=1:T
         assertw( (nvert(j)~=0), 'chkmesh 10' );
         assertw( (sinring(VERT,SIMP,i,nvert(j))>0), 'chkmesh 11' );
      end

      %%% check simplex rings and nabors
      for j=1:T
         if (j==1) k=2; l=3; end;
         if (j==2) k=3; l=1; end;
         if (j==3) k=2; l=1; end;

         %%% check dirichlet points on dirichlet faces
         if ( vdiri(ftype(j)) )
            bndFaceD = bndFaceD+1;
            assertw( vdiri(VERT(SIMP(i,k),6)), 'chkmesh 6' );
            assertw( vdiri(VERT(SIMP(i,l),6)), 'chkmesh 7' );
         end

         %%% check neumann points on neumann faces
         if ( vneum(ftype(j)) )
            bndFaceN = bndFaceN+1;
            assertw( vbnd(VERT(SIMP(i,k),6)), 'chkmesh 8' );
            assertw( vbnd(VERT(SIMP(i,l),6)), 'chkmesh 9' );
         end

         %%% check face nabors; verify we have a conforming mesh
         [nabor,cface] = getnabor(VERT,SIMP,i,j);
         cond1 = (vbnd(ftype(j)) & (nabor==0));
         cond2 = (vint(ftype(j)) & (nabor>0));
         assertw( (cond1 | cond2), 'chkmesh 12' );
      end
   end
   the_time = etime(clock,t0);
   fprintf('..done.  [time=%g]\n', the_time);

%%% destroy the structures

   [VERT,SIMP] = kilring(VERT,SIMP);

%%% return

   fprintf('chkmesh: finished conformity and consistency check of mesh.\n');
   fprintf('chkmesh:    point summary:        [T=%d,I=%d,D=%d,N=%d]\n', ...
      N,vertInt,vertBndD,vertBndN);
   fprintf('chkmesh:    bndry face summary:   [D=%d,N=%d]\n', ...
      bndFaceD,bndFaceN);

   value = 1;

