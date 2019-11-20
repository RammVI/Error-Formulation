function [V,D] = myeig(A)
%MYEIG  A simple QR itation
%
% Usage: [V,D] = myeig(A)
%
% Author:   Michael Holst
% rcsid="$Id: myeig.m,v 1.1.1.1 2007/04/27 08:28:19 hrg Exp $"

    Qsum = eye(size(A));
    for i=1:1000
        [Q,R] = qr(A);
        A = R * Q;
        Qsum = Qsum * Q;
    end;

    V = Qsum;
    D = A;

end;

