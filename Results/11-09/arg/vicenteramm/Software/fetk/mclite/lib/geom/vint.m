function [result] = vint(type);
%VINT  Yes or no answer to whether triang/edge/vertex is "interior" type
%
% Usage: [result] = vint(type);
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
% rcsid="$Id: vint.m,v 1.1.1.1 2007/04/27 08:28:07 hrg Exp $"

   [m,n] = size(type);

   result = (type==0);

