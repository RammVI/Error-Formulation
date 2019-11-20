function [thetam] = fu(evalKey, parm, xm,ym,type, nvec, u,ux,uy);
%FU  Evaluate the nonlinear strong form of a PDE
%
% Usage: [thetam] = fu(evalKey, parm, xm,ym,type, u,ux,uy, nvec);
%
%     Evaluate the nonlinear form by providing the value of the quantity
%     at a particular point.
%
% Input:
%
%     evalKey     ==> particular piece of the strong form that is desired
%     xm          ==> x-coordindinates of the point
%     ym          ==> y-coordindinates of the point
%     type        ==> point type (interior/diri/neum)
%
%     u           ==> u, evaluated on master element
%     ux          ==> partial w.r.t. x of u
%     uy          ==> partial w.r.t. x of u
%
%     nvec(1:2)   ==> unit outward normal vector if an element face
%
% Output:
%
%     thetam      ==> the value of the strong form at the point
%                     not including the jacobian of transformation
%
% Author:   Michael Holst
% rcsid="$Id: fu_1.m,v 1.1.1.1 2007/04/27 08:28:18 hrg Exp $"

   [n,one] = size(u);
   assertw( n==2, 'fu 1' );

   %%% elasticity parameters:  specify (E,nu) to yield (lambda,mu)
   E        = 21.0;     %%% Young's modulus
   nu       = 0.28;     %%% Poisson ratio
   lambda   = (E * nu) / ( (1+nu) * (1-2*nu) );
   mu       = E / ( 2 * (1+nu) );

   thetam = zeros(n,one);

   %%% element residual case:  b(u)^i - a(u)^{iq}_{~;q}
   if (evalKey == 0)

   %%% neumann face residual case:  c(u)^i + a(u)^{iq} n_q
   elseif (evalKey == 1)

   %%% interior face residual case:  a(u)^{iq} n_q
   elseif (evalKey == 2)

      thetam(1) = ux(1) * nvec(1) + uy(1) * nvec(2);
      thetam(2) = ux(2) * nvec(1) + uy(2) * nvec(2);

   %%% error
   else
      assertw( 0, 'fu 2' );
   end

