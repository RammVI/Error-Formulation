function [U,Ug,A,F]=lsolve(VERT, SIMP);
%NEWTON  Discretize and solve the linear pde
%
% Usage: [U,Ug,A,F]=lsolve(VERT, SIMP);
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
% rcsid="$Id: lsolve.m,v 1.1.1.1 2007/04/27 08:28:19 hrg Exp $"   

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

   %%% assemble linear system
   assKey = 2;
   parm(1) = 0;
   [A,F]=assem(assKey, parm, U, Ug, VERT, SIMP);

   %%% solve system -- direct method
   U = A \ F;

   %%% residual
   R = F - A*U;

   %%% solution quality
   normR = norm(R,2);
   fprintf('LSOLVE: ||R|| = %g\n', normR);

