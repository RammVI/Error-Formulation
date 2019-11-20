% Domain number: 2
%
% Dimension:     2
% Domain:        Two-dimensional bar.
% Boundary:      One dirichlet face at each end of the bar.
%
% Notes:         The elements are split into three material types (0,1,2).
%                This is basically to tag the two triangles in the middle
%                section (type 1) and the two at the right end (type 2),
%                for applying (body) forces.  All of the other elements
%                are marked as the third type (type 0).
%
%                The boundary faces are split into four types.
%                The left end of the bar is Dirichlet (type 1), so it fixes 
%                the bar to an imaginary wall on the left.  The rest of the 
%                boundary faces are normal Neumann (type 2), except for
%                the face in the middle of the bar on the top (type 4),
%                and the face at the right end on the bottom (type 6).
%                These two special faces are for applying 
%                (traction) forces.
%
% Author:   Michael Holst
% rcsid="$Id: mcin_1.m,v 1.1.1.1 2007/04/27 08:28:18 hrg Exp $"

mcsf_begin=1;

      dim=2;         % intrinsic manifold dimension
    dimii=3;         % imbedding manifold dimension
 vertices=22;        % number of vertices
simplices=20;        % number of simplices

vert=[
%-------- ---- ----------------- ----------------- -----------------
% Vert-ID Chrt X-Coordinate      Y-Coordinate      Z-Coordinate
%-------- ---- ----------------- ----------------- -----------------
0         0     0.0000000000e+00  0.0000000000e+00  0.0000000000e+00
1         0     1.0000000000e+00  0.0000000000e+00  0.0000000000e+00
2         0     2.0000000000e+00  0.0000000000e+00  0.0000000000e+00
3         0     3.0000000000e+00  0.0000000000e+00  0.0000000000e+00
4         0     4.0000000000e+00  0.0000000000e+00  0.0000000000e+00
5         0     5.0000000000e+00  0.0000000000e+00  0.0000000000e+00
6         0     6.0000000000e+00  0.0000000000e+00  0.0000000000e+00
7         0     7.0000000000e+00  0.0000000000e+00  0.0000000000e+00
8         0     8.0000000000e+00  0.0000000000e+00  0.0000000000e+00
9         0     9.0000000000e+00  0.0000000000e+00  0.0000000000e+00
10        0     1.0000000000e+01  0.0000000000e+00  0.0000000000e+00
11        0     0.0000000000e+00  1.0000000000e+00  0.0000000000e+00
12        0     1.0000000000e+00  1.0000000000e+00  0.0000000000e+00
13        0     2.0000000000e+00  1.0000000000e+00  0.0000000000e+00
14        0     3.0000000000e+00  1.0000000000e+00  0.0000000000e+00
15        0     4.0000000000e+00  1.0000000000e+00  0.0000000000e+00
16        0     5.0000000000e+00  1.0000000000e+00  0.0000000000e+00
17        0     6.0000000000e+00  1.0000000000e+00  0.0000000000e+00
18        0     7.0000000000e+00  1.0000000000e+00  0.0000000000e+00
19        0     8.0000000000e+00  1.0000000000e+00  0.0000000000e+00
20        0     9.0000000000e+00  1.0000000000e+00  0.0000000000e+00
21        0     1.0000000000e+01  1.0000000000e+00  0.0000000000e+00
];

simp=[
%-------- ---- ---- ------------------- ---------------------------------------
% Simp-ID Grp  Mat  Face-Types          Vertex-Numbers
%-------- ---- ---- ------------------- ---------------------------------------
0         0    0    0    1    2         0 1 11
1         0    0    2    0    0         1 12 11
2         0    0    0    0    2         1 2 12
3         0    0    2    0    0         2 13 12
4         0    0    0    0    2         2 3 13
5         0    0    2    0    0         3 14 13
6         0    0    0    0    2         3 4 14
7         0    0    2    0    0         4 15 14
8         0    0    0    0    2         4 5 15
9         0    0    2    0    0         5 16 15
10        0    1    0    0    2         5 6 16
11        0    1    4    0    0         6 17 16
12        0    0    0    0    2         6 7 17
13        0    0    2    0    0         7 18 17
14        0    0    0    0    2         7 8 18
15        0    0    2    0    0         8 19 18
16        0    0    0    0    2         8 9 19
17        0    0    2    0    0         9 20 19
18        0    2    0    0    6         9 10 20
19        0    2    2    0    1         10 21 20
];

mcsf_end=1;

