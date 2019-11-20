function d(deg);
%D  Rotate the current plot downward a specified number of degrees
%
% Usage: d(deg);
%
% Author:   Michael Holst
% rcsid="$Id: dd.m,v 1.1.1.1 2007/04/27 08:28:19 hrg Exp $"

[AZ,EL] = view;
EL = EL + deg;
view(AZ,EL);

