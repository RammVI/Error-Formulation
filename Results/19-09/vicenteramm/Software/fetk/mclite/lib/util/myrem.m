function [value] = myrem(x,y)
%MYMOD  Simple arithmetic "remainder" operation
%
% Usage: myrem(x,y);
%
% Author:   Michael Holst
% rcsid="$Id: myrem.m,v 1.1.1.1 2007/04/27 08:28:19 hrg Exp $"

if (y == 0)
    value = NaN;
else
    value = x - y.*fix(x./y);
end

