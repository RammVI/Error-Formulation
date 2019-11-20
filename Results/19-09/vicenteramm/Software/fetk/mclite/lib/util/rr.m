function r(deg);
%R  Rotate the current plot rightward a specified number of degrees
%
% Usage: r(deg);
%
% Author:   Michael Holst
% rcsid="$Id: rr.m,v 1.1.1.1 2007/04/27 08:28:19 hrg Exp $"

[AZ,EL] = view;
AZ = AZ - deg;
view(AZ,EL);

