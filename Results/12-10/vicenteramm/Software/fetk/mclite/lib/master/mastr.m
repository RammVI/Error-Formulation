function [INTPT,  W,  PHI,  PHIX,  PHIY,   ...
          INTPTE0,WE0,PHIE0,PHIXE0,PHIYE0, ...
          INTPTE1,WE1,PHIE1,PHIXE1,PHIYE1, ...
          INTPTE2,WE2,PHIE2,PHIXE2,PHIYE2] = mastr(basis,Qorder);
%MASTR  Define the master element and master element basis function info
%
% Usage: [INTPT,  W,  PHI,  PHIX,  PHIY,   ...
%         INTPTE0,WE0,PHIE0,PHIXE0,PHIYE0, ...
%         INTPTE1,WE1,PHIE1,PHIXE1,PHIYE1, ...
%         INTPTE2,WE2,PHIE2,PHIXE2,PHIYE2] = mastr(basis,Qorder);
%
% Input:
%
%    basis  = basis function choice (string; name of script)
%    Qorder = quadrature order choice (exact for polys of degree "Qorder")
%
% Output:
%
%    (T  = number of nodes in a single element)
%    (Q  = number of quadrature points in the interior of the master element)
%    (QE = number of quadrature points along the master element edges)
%
%    INTPT   = zeros(Q,2)  = master element volume integration points
%    W       = zeros(Q,1)  = master element volume integration weights
%    PHI     = zeros(T,Q)  = values of local basis funcs at integ pts
%    PHIX    = zeros(T,Q)  = values of x-derivs of local basis funcs integ pts
%    PHIY    = zeros(T,Q)  = values of y-derivs of local basis funcs integ pts
%
%    INTPTE0 = zeros(QE,2) = type (1,2) EDGE integration points
%    WE0     = zeros(QE,1) = type (1,2) EDGE integration weights
%    PHIE0   = zeros(T,QE) = type (1,2) basis func values at EDGE integ pts
%    PHIXE0  = zeros(T,QE) = type (1,2) basis func x-derivs at EDGE integ pts
%    PHIYE0  = zeros(T,QE) = type (1,2) basis func y-derivs at EDGE integ pts
%
%    INTPTE1 = zeros(QE,2) = type (2,0) EDGE integration points
%    WE1     = zeros(QE,1) = type (2,0) EDGE integration weights
%    PHIE1   = zeros(T,QE) = type (2,0) basis function values at EDGE integ pts
%    PHIXE1  = zeros(T,QE) = type (2,0) basis func x-derivs at EDGE integ pts
%    PHIYE1  = zeros(T,QE) = type (2,0) basis func y-derivs at EDGE integ pts
%
%    INTPTE2 = zeros(QE,2) = type (0,1) EDGE integration points
%    WE2     = zeros(QE,1) = type (0,1) EDGE integration weights
%    PHIE2   = zeros(T,QE) = type (0,1) basis function values at EDGE integ pts
%    PHIXE2  = zeros(T,QE) = type (0,1) basis func x-derivs at EDGE integ pts
%    PHIYE2  = zeros(T,QE) = type (0,1) basis func y-derivs at EDGE integ pts
%
% Author:   Michael Holst
% rcsid="$Id: mastr.m,v 1.1.1.1 2007/04/27 08:28:08 hrg Exp $"

   [INTPT,  W,  PHI,  PHIX,  PHIY,   ...
    INTPTE0,WE0,PHIE0,PHIXE0,PHIYE0, ...
    INTPTE1,WE1,PHIE1,PHIXE1,PHIYE1, ...
    INTPTE2,WE2,PHIE2,PHIXE2,PHIYE2] = feval(basis,Qorder);

