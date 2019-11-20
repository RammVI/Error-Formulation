function [result] = vbnd(type);
%VBND  Yes or no answer to whether triang/edge/vertex is "boundary" type
%
% Usage: [result] = vbnd(type);
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
% rcsid="$Id: vbnd.m,v 1.1.1.1 2007/04/27 08:28:07 hrg Exp $"

   [m,n] = size(type);

   result = (type~=0);

