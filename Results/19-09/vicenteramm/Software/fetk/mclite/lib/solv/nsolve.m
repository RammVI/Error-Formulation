function [U,Ug,A,F]=nsolve(VERT, SIMP);
%NEWTON  Discretize and solve the nonlinear pde with a newton iteration
%
% Usage: [U,Ug,A,F]=nsolve(VERT, SIMP);
%
% Input:
%
%    VERT   = <see "read.m" for description of the datastructure>
%    SIMP   = <see "read.m" for description of the datastructure>
% Output:
%
%    U    = solution to nonlinear problem
%    Ug   = essential/dirichlet boundary function
%    A    = last jacobian matrix
%    F    = last nonlinear residual vector
%
% Author:   Michael Holst
% rcsid="$Id: nsolve.m,v 1.1.1.1 2007/04/27 08:28:19 hrg Exp $"   

   %%% sizes of things
   [N,eight]     = size(VERT);
   [L,seventeen] = size(SIMP);

   %%% initialize the dirichlet function
   Ugg = u_d( VERT(1:N,1), VERT(1:N,2), VERT(1:N,6) );
   [NN,n] = size(Ugg);
   diriMask = vdiri(VERT(1:N,6));
   theMask = [];
   for i=1:n
      theMask = [ theMask , diriMask ];
   end
   Ugg = Ugg .* theMask;
   Ug = reshape( Ugg, N*n, 1 );

   %%% initialize the solution
   U = zeros(N*n,1);

   %%% parameters
   parm = zeros(10,1);

   %%% newton iteration
   itmax = 30;
   errtol = 1.0e-7;
   iter = 0;
   done = 0;
   while ((done == 0) & (iter < itmax))
      iter = iter + 1;

      %%% assemble jacobian system
      assKey = 2;
      parm(1) = 0;
      [A,F]=assem(assKey, parm, U, Ug, VERT, SIMP);

      %%% initial residual
      if (iter == 1)
         normF = norm(F,2);
         fprintf('NEWTON: iter = %g, ||F(u)|| = %g\n', 0, normF);
      end

      %%% solve jacobian system -- direct method
      E = A \ F;

      %%% correct the solution
      U = U + E;

      %%% nonlinear residual
      assKey = 3;
      parm(1) = 0;
      [Adummy,F]=assem(assKey, parm, U, Ug, VERT, SIMP);

      %%% solution quality
      normF = norm(F,2);
      fprintf('NEWTON: iter = %g, ||F(u)|| = %g\n', iter, normF);
      if (normF < errtol)
         done = 1;
      end;
   end

