function drawf(VERT,SIMP,USOL);
%DRAWF  Draw the finite element function
%
% Usage: drawf(VERT,SIMP,USOL);
%
% Author:   Michael Holst
% rcsid="$Id: drawf.m,v 1.1.1.1 2007/04/27 08:28:07 hrg Exp $"

%%% recover various problem dimensions

   [N,eight]     = size(VERT);
   [L,seventeen] = size(SIMP);
   T = 3;

%%% setup for plot

   hold off;

%%% cycle through the simplices and draw the polygons

   C = [0.5 0.5 0.5];
   X = zeros(T,L);
   Y = zeros(T,L);
   Z = zeros(T,L);
   for simplex = 1 : L
      X(1:T,simplex) = VERT(SIMP(simplex,1:T),1);
      Y(1:T,simplex) = VERT(SIMP(simplex,1:T),2);
      Z(1:T,simplex) = USOL(SIMP(simplex,1:T));
   end;
   fill3(X,Y,Z,C);
   axis equal;
   axis off;
   hold on;

%%% plot the essential boundary vertices with circles

   for simplex = 1 : L
      %%% plot the essential boundary vertices with circles
      for i = 1:T
         vtype = VERT(SIMP(simplex,i),6);
         if ( vdiri(vtype) )
            plot3(VERT(SIMP(simplex,i),1), ...
                  VERT(SIMP(simplex,i),2), ...
                  USOL(SIMP(simplex,i)), 'or');
         end
      end
   end

