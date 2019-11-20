function [value] = mymod(x,y)
%MYMOD  Simple arithmetic "modulus" operation
%
% Usage: mymod(x,y);
%
% Author:   Michael Holst
% rcsid="$Id: mymod.m,v 1.1.1.1 2007/04/27 08:28:19 hrg Exp $"

if (y == 0)
    value = x;
else
    value = x - y.*floor(x./y);
end

