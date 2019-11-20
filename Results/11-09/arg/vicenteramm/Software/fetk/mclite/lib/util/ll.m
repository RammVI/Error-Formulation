function l(deg);
%L  Rotate the current plot leftward a specified number of degrees
%
% Usage: l(deg);
%
% Author:   Michael Holst
% rcsid="$Id: ll.m,v 1.1.1.1 2007/04/27 08:28:19 hrg Exp $"

[AZ,EL] = view;
AZ = AZ + deg;
view(AZ,EL);

