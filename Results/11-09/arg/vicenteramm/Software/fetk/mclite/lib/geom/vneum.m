function [result] = vneum(type);
%VNEUM  Yes or no answer to whether triang/edge/vertex is "neumann" type
%
% Usage: [result] = vneum(type);
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
% rcsid="$Id: vneum.m,v 1.1.1.1 2007/04/27 08:28:07 hrg Exp $"

   [m,n] = size(type);

   even = type ./ ( 2 * ones(m,n) );
   evenint = floor(even);
   eventst = ( (type > 0) & (~abs(even - evenint)) );

   result = (eventst~=0);

