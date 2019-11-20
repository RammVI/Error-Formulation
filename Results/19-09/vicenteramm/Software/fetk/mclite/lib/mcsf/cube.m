% tetrahedral mesh file (MCSF format)
%
% Domain is unit cube.  Faces in the plane x=0 are labelled 1, those
% in the plane x=1 are labelled 2, and other boundary faces are labelled 3.
%
% Author:   Doug Arnold
% rcsid="$Id: cube.m,v 1.1.1.1 2007/04/27 08:28:08 hrg Exp $"

mcsf_begin=1;

dim=           3;    % intrinsic manifold dimension
dimii=         3;    % imbedding manifold dimension
vertices=      9;    % number of vertices
simplices=    12;    % number of simplices

vert=[
%-------- ---- ----------------- ----------------- -----------------
% Vert-ID Chrt X-Coordinate      Y-Coordinate      Z-Coordinate
%-------- ---- ----------------- ----------------- -----------------
        0    0  0.0000000000e+00  0.0000000000e+00  0.0000000000e+00
        1    0  1.0000000000e+00  0.0000000000e+00  0.0000000000e+00
        2    0  1.0000000000e+00  1.0000000000e+00  1.0000000000e+00
        3    0  1.0000000000e+00  0.0000000000e+00  1.0000000000e+00
        4    0  0.0000000000e+00  1.0000000000e+00  1.0000000000e+00
        5    0  0.0000000000e+00  0.0000000000e+00  1.0000000000e+00
        6    0  0.0000000000e+00  1.0000000000e+00  0.0000000000e+00
        7    0  1.0000000000e+00  1.0000000000e+00  0.0000000000e+00
        8    0  5.0000000000e-01  5.0000000000e-01  5.0000000000e-01
];

simp=[
%-------- ---- ---- ------------------- ---------------------------------------
% Simp-ID Grp  Mat  Face-Types          Vertex-Numbers
%-------- ---- ---- ------------------- ---------------------------------------
        0    0    0    0    0    0    3     6     2     4     8
        1    0    0    0    0    0    3     1     5     3     8
        2    0    0    0    0    3    0     1     5     8     0
        3    0    0    1    0    0    0     8     4     0     5
        4    0    0    0    0    0    1     6     4     0     8
        5    0    0    3    0    0    0     8     2     5     3
        6    0    0    0    0    3    0     6     2     8     7
        7    0    0    3    0    0    0     8     1     0     7
        8    0    0    0    3    0    0     6     8     0     7
        9    0    0    2    0    0    0     8     1     7     3
       10    0    0    0    0    0    3     4     2     5     8
       11    0    0    0    2    0    0     2     8     7     3
];

mcsf_end=1;
