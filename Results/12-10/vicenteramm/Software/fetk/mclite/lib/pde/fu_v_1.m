function [thetam] = fu_v(evalKey, parm, L, xm,ym,type, u,ux,uy, pr,prx,pry);
%FU_V  Evaluate the nonlinear weak form of a PDE
%
% Usage: [thetam] = fu_v(evalKey, parm, L, xm,ym,type, u,ux,uy, pr,prx,pry);
%
%     Evaluate the nonlinear form by providing the value of the quantity
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
% rcsid="$Id: fu_v_1.m,v 1.1.1.1 2007/04/27 08:28:18 hrg Exp $"

   [L,n] = size(u);
   assertw( n==2, 'fu_v 1' );

   %%% elasticity parameters:  specify (E,nu) to yield (lambda,mu)
   E        = 21.0;     %%% Young's modulus
   nu       = 0.28;     %%% Poisson ratio
   lambda   = (E * nu) / ( (1+nu) * (1-2*nu) );
   mu       = E / ( 2 * (1+nu) );

   thetam = zeros(L,2);

   %%% area integral
   if (evalKey == 0)

      ars = zeros(L,4);
      ars(1:L,1) = ...
          (2*mu+lambda)*ux(1:L,1).*prx(1:L,1) + mu*uy(1:L,1).*pry(1:L,1);
      ars(1:L,4) = ...
          (2*mu+lambda)*uy(1:L,2).*pry(1:L,1) + mu*ux(1:L,2).*prx(1:L,1);
      ars(1:L,3) = ...
          mu*uy(1:L,1).*prx(1:L,1) + lambda*ux(1:L,1).*pry(1:L,1);
      ars(1:L,2) = ...
          mu*ux(1:L,2).*pry(1:L,1) + lambda*uy(1:L,2).*prx(1:L,1);

      ftilde = zeros(L,2);
      % ftilde(1:L,1) = zeros(L,1);
      % ftilde(1:L,2) = - 1.0 * (type==1*ones(L,1));

      thetam(1:L,1) = ars(1:L,1) + ars(1:L,2) - ftilde(1:L,1)*pr;
      thetam(1:L,2) = ars(1:L,3) + ars(1:L,4) - ftilde(1:L,2)*pr;

   %%% edge integral (for natural boundary condition)
   elseif (evalKey == 1)

      ftilde = zeros(L,2);
      ftilde(1:L,1) = zeros(L,1);
      ftilde(1:L,2) = - 1.0 * (type==4*ones(L,1));

      thetam(1:L,1) = - ftilde(1:L,1)*pr;
      thetam(1:L,2) = - ftilde(1:L,2)*pr;

   %%% error
   else
      assertw( 0, 'fu_v 2' );
   end

