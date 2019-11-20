function [thetam] = dfu_wv(evalKey, parm, L, xm,ym,type, u,ux,uy, ...
   ps,psx,psy, pr,prx,pry);
%DFU_WV  evaluate bilinear linearization form of a PDE
%
% Usage: [thetam] = dfu_wv(evalKey, parm, L, xm,ym,type, u,ux,uy, ...
%                   ps,psx,psy, pr,prx,pry);
%
%     Evaluate the bilinear form by providing the value of the quantity
%     under the integral sign, at a particular quadrature point.
%
%     Because we do this for all of the elements at one time to make things
%     reasonably fast, we pass some arrays around.
%
% Input:
%
%     evalKey     ==> area or surface integral
%     L           ==> number of elements to do all at once
%     xm(1:L,1)   ==> x-coordindinates of quadrature points, one per element
%     ym(1:L,1)   ==> y-coordindinates of quadrature points, one per element
%     type(1:L,1) ==> point type (interior/diri/neum), one per element
%
%     u(1:L,1)    ==> u, evaluated on master element
%     ux(1:L,1)   ==> partial w.r.t. x of u, one per element
%     uy(1:L,1)   ==> partial w.r.t. x of u, one per element
%
%     ps          ==> phi(s), evaluated on master element
%     psx(1:L,1)  ==> partial w.r.t. x of phi(s), one per element
%     psy(1:L,1)  ==> partial w.r.t. x of phi(s), one per element
%
%     pr          ==> phi(r), evaluated on master element
%     prx(1:L,1)  ==> partial w.r.t. x of phi(r), one per element
%     pry(1:L,1)  ==> partial w.r.t. x of phi(r), one per element
%
% Output:
%
%     thetam      ==> the integrand for the area or surface integral,
%                     not including the jacobian of transformation
%
% Author:   Michael Holst
% rcsid="$Id: dfu_wv_1.m,v 1.1.1.1 2007/04/27 08:28:18 hrg Exp $"

   [L,n] = size(u);
   assertw( n==2, 'dfu_wv 1' );

   %%% elasticity parameters:  specify (E,nu) to yield (lambda,mu)
   E        = 21.0;     %%% Young's modulus
   nu       = 0.28;     %%% Poisson ratio
   lambda   = (E * nu) / ( (1+nu) * (1-2*nu) );
   mu       = E / ( 2 * (1+nu) );

   thetam = zeros(L,2*2);

   %%% area integral
   if (evalKey == 0)

      thetam(1:L,1) = ...
          (2*mu+lambda)*psx(1:L,1).*prx(1:L,1) + mu*psy(1:L,1).*pry(1:L,1);
      thetam(1:L,4) = ...
          (2*mu+lambda)*psy(1:L,1).*pry(1:L,1) + mu*psx(1:L,1).*prx(1:L,1);
      thetam(1:L,3) = ...
          mu*psy(1:L,1).*prx(1:L,1) + lambda*psx(1:L,1).*pry(1:L,1);
      thetam(1:L,2) = ...
          mu*psx(1:L,1).*pry(1:L,1) + lambda*psy(1:L,1).*prx(1:L,1);

   %%% edge integral (for natural boundary condition)
   elseif (evalKey == 1)

   %%% error
   else
      assertw( 0, 'dfu_wv 2' );
   end

