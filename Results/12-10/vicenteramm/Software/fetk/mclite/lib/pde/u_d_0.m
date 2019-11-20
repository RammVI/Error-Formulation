function [gxy] = u_d(x,y,type);
%U_D  Dirichlet boundary condition function
%
% Usage: [gxy] = u_d(x,y,type);
%
% Input:
%
%     xm(1:L,1)   ==> x-coordindinates of evaluation points
%     ym(1:L,1)   ==> y-coordindinates of evaluation points
%     type(1:L,1) ==> point type (interior/diri/neum)
%
% Output:
%
%     gxy         ==> value of the dirichlet function at evaluation points
%
% Author:   Michael Holst
% rcsid="$Id: u_d_0.m,v 1.1.1.1 2007/04/27 08:28:18 hrg Exp $"

   [M,one] = size(x);

   gxy = zeros(M,1);
   gxy = u_t(x,y,type);

