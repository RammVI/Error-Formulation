function drawd(VERT,SIMP,USOL);
%DRAW  Draw the finite element solution as a mesh displacement
%
% Usage: drawd(VERT,SIMP,USOL);
%
% Author:   Michael Holst
% rcsid="$Id: drawd.m,v 1.1.1.1 2007/04/27 08:28:07 hrg Exp $"

%%% recover various problem dimensions

   [N,eight]     = size(VERT);
   [L,seventeen] = size(SIMP);
   T = 3;

%%% displace the mesh by the solution

   NODE = zeros(N,3);
   NODE(:,1:2) = VERT(:,1:2) + reshape(USOL,N,2);

%%% setup for plot

   hold off;

%%% plot all vertices with dots

   plot(NODE(1:N,1),NODE(1:N,2), '.w')
   axis equal;
   axis off;
   hold on;

%%% cycle through the simplices and draw the edges

   for ell = 1:L
      Xil(1:T,1:3) = NODE(SIMP(ell,1:T),1:3);
      plot([Xil(1:T,1)' Xil(1,1)], [Xil(1:T,2)' Xil(1,2)], 'w');

      %%% face info for the T face making up this simplex
      ftype = SIMP(ell,5:7);

      %%% plot the essential boundary vertices with circles
      for i = 1:T
         j=i+1;
         if (j>T) j=1; end;
         k=j+1;
         if (k>T) k=1; end;

         vtype = VERT(SIMP(ell,i),6);
         if ( vdiri(vtype) )
            plot(NODE(SIMP(ell,i),1), NODE(SIMP(ell,i),2), 'oy');
         end
         if ( vbnd(ftype(i)) )
            Xil(1,1:3) = NODE(SIMP(ell,j),1:3);
            Xil(2,1:3) = NODE(SIMP(ell,k),1:3);
            if ( vdiri(ftype(i)) )
               plot(Xil(1:2,1), Xil(1:2,2), 'r');
            else
               plot(Xil(1:2,1), Xil(1:2,2), 'g');
            end
         end
      end
   end

