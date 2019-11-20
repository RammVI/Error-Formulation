function [rhoxy] = rho(x,y);
%RHO  Mass function for computing the mass matrix.
%
% Usage: [rhoxy] = rho(x,y);
%
% Author:   Michael Holst
% rcsid="$Id: rho_0.m,v 1.1.1.1 2007/04/27 08:28:18 hrg Exp $"

   [L,one] = size(x);
   rhoxy = ones(L,1);

