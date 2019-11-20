function [result] = vdiri(type);
%VDIRI  Yes or no answer to whether triang/edge/vertex is "dirichlet" type
%
% Usage: [result] = vdiri(type);
%
% Input:
%
%    type   = the type in question
% 
% Output:
%    
%    result = yes or no
%
% Author:   Michael Holst
% rcsid="$Id: vdiri.m,v 1.1.1.1 2007/04/27 08:28:07 hrg Exp $"

   [m,n] = size(type);

   odd = type ./ ( 2 * ones(m,n) );
   oddint = floor(odd);
   oddtst = abs(odd - oddint);

   result = (oddtst~=0);

