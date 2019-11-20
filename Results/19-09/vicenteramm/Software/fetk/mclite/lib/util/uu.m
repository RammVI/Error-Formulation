function u(deg);
%U  Rotate the current plot upward a specified number of degrees
%
% Usage: u(deg);
%
% Author:   Michael Holst
% rcsid="$Id: uu.m,v 1.1.1.1 2007/04/27 08:28:19 hrg Exp $"

[AZ,EL] = view;
EL = EL - deg;
view(AZ,EL);

