function [value] = assertw(condition,name);
%ASSERT  Assertion macro for debugging MATLAB codes
%
% Usage: [value] = assertw(condition);
%
% Input:
%
%    condition = zero or nonzero
%
% Output:
%
%    value = condition
%    AND a print statement about an error of condition!=1
%
% Author:   Michael Holst
% rcsid="$Id: assertw.m,v 1.1.1.1 2007/04/27 08:28:19 hrg Exp $"

if (condition == 0)
   fprintf('***** ASSERTION FAILURE AT: <%s> *****\n', name);
end

value = condition;

