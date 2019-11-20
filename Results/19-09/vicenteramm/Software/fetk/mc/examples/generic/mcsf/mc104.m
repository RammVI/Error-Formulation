%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File:          mc104.m
%
% Dimension:     3
% Domain:        3D version of an L-shaped domain.
% Boundary:      All bndry nodes are dirichlet-type.
%
% rcsid="$Id: mc104.m,v 1.4 2010/08/12 05:17:16 fetk Exp $"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MC = < Manifold Code >
%%% Copyright (C) 1994-- Michael Holst 
%%%
%%% This library is free software; you can redistribute it and/or 
%%% modify it under the terms of the GNU Lesser General Public 
%%% License as published by the Free Software Foundation; either  
%%% version 2.1 of the License, or (at your option) any later version. 
%%%
%%% This library is distributed in the hope that it will be useful, 
%%% but WITHOUT ANY WARRANTY; without even the implied warranty of             
%%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
%%% Lesser General Public License for more details. 
%%%
%%% You should have received a copy of the GNU Lesser General Public 
%%% License along with this library; if not, write to the Free Software  
%%% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA   
%%%
%%% rcsid="$Id: mc104.m,v 1.4 2010/08/12 05:17:16 fetk Exp $"
%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mcsf_begin=1;

      dim=3;         % intrinsic manifold dimension
    dimii=3;         % imbedding manifold dimension
 vertices=152;       % number of vertices
simplices=404;       % number of simplices

vert=[
%-------- ---- ----------------- ----------------- -----------------
% Vert-ID Chrt X-Coordinate      Y-Coordinate      Z-Coordinate
%-------- ---- ----------------- ----------------- -----------------
0         0     0.0000000000e+00  0.0000000000e+00  0.0000000000e+00
1         0     1.0000000000e+00  0.0000000000e+00  0.0000000000e+00
2         0     0.0000000000e+00  1.0000000000e+00  0.0000000000e+00
3         0     0.0000000000e+00  0.0000000000e+00  1.0000000000e+00
4         0    -1.0000000000e+00  0.0000000000e+00  0.0000000000e+00
5         0     0.0000000000e+00 -1.0000000000e+00  0.0000000000e+00
6         0     0.0000000000e+00  0.0000000000e+00 -1.0000000000e+00
7         0     1.0000000000e+00  0.0000000000e+00 -1.0000000000e+00
8         0    -1.0000000000e+00  0.0000000000e+00  1.0000000000e+00
9         0     1.0000000000e+00 -1.0000000000e+00  0.0000000000e+00
10        0     0.0000000000e+00 -1.0000000000e+00  1.0000000000e+00
11        0     0.0000000000e+00  1.0000000000e+00 -1.0000000000e+00
12        0    -1.0000000000e+00  1.0000000000e+00  0.0000000000e+00
13        0     1.0000000000e+00 -1.0000000000e+00  1.0000000000e+00
14        0     1.0000000000e+00  1.0000000000e+00 -1.0000000000e+00
15        0    -1.0000000000e+00  1.0000000000e+00  1.0000000000e+00
16        0     1.0000000000e+00  1.0000000000e+00  1.0000000000e+00
17        0     0.0000000000e+00  1.0000000000e+00  1.0000000000e+00
18        0     1.0000000000e+00  1.0000000000e+00  0.0000000000e+00
19        0     3.3333330000e-01  0.0000000000e+00  0.0000000000e+00
20        0     6.6666670000e-01  0.0000000000e+00  0.0000000000e+00
21        0     1.0000000000e+00  0.0000000000e+00 -3.3333330000e-01
22        0     1.0000000000e+00  0.0000000000e+00 -6.6666670000e-01
23        0     6.6666670000e-01  0.0000000000e+00 -1.0000000000e+00
24        0     3.3333330000e-01  0.0000000000e+00 -1.0000000000e+00
25        0     0.0000000000e+00  0.0000000000e+00 -6.6666670000e-01
26        0     0.0000000000e+00  0.0000000000e+00 -3.3333330000e-01
27        0     0.0000000000e+00  3.3333330000e-01 -1.0000000000e+00
28        0     0.0000000000e+00  6.6666670000e-01 -1.0000000000e+00
29        0     0.0000000000e+00  1.0000000000e+00 -6.6666670000e-01
30        0     0.0000000000e+00  1.0000000000e+00 -3.3333330000e-01
31        0     0.0000000000e+00  6.6666670000e-01  0.0000000000e+00
32        0     0.0000000000e+00  3.3333330000e-01  0.0000000000e+00
33        0    -3.3333330000e-01  1.0000000000e+00  0.0000000000e+00
34        0    -6.6666670000e-01  1.0000000000e+00  0.0000000000e+00
35        0    -1.0000000000e+00  6.6666670000e-01  0.0000000000e+00
36        0    -1.0000000000e+00  3.3333330000e-01  0.0000000000e+00
37        0    -6.6666670000e-01  0.0000000000e+00  0.0000000000e+00
38        0    -3.3333330000e-01  0.0000000000e+00  0.0000000000e+00
39        0    -1.0000000000e+00  0.0000000000e+00  3.3333330000e-01
40        0    -1.0000000000e+00  0.0000000000e+00  6.6666670000e-01
41        0    -6.6666670000e-01  0.0000000000e+00  1.0000000000e+00
42        0    -3.3333330000e-01  0.0000000000e+00  1.0000000000e+00
43        0     0.0000000000e+00  0.0000000000e+00  6.6666670000e-01
44        0     0.0000000000e+00  0.0000000000e+00  3.3333330000e-01
45        0     0.0000000000e+00 -3.3333330000e-01  1.0000000000e+00
46        0     0.0000000000e+00 -6.6666670000e-01  1.0000000000e+00
47        0     0.0000000000e+00 -1.0000000000e+00  6.6666670000e-01
48        0     0.0000000000e+00 -1.0000000000e+00  3.3333330000e-01
49        0     0.0000000000e+00 -6.6666670000e-01  0.0000000000e+00
50        0     0.0000000000e+00 -3.3333330000e-01  0.0000000000e+00
51        0     3.3333330000e-01 -1.0000000000e+00  0.0000000000e+00
52        0     6.6666670000e-01 -1.0000000000e+00  0.0000000000e+00
53        0     1.0000000000e+00 -6.6666670000e-01  0.0000000000e+00
54        0     1.0000000000e+00 -3.3333330000e-01  0.0000000000e+00
55        0     1.0000000000e+00  3.3333330000e-01 -1.0000000000e+00
56        0     1.0000000000e+00  6.6666670000e-01 -1.0000000000e+00
57        0     6.6666670000e-01  1.0000000000e+00 -1.0000000000e+00
58        0     3.3333330000e-01  1.0000000000e+00 -1.0000000000e+00
59        0    -1.0000000000e+00  1.0000000000e+00  3.3333330000e-01
60        0    -1.0000000000e+00  1.0000000000e+00  6.6666670000e-01
61        0    -1.0000000000e+00  6.6666670000e-01  1.0000000000e+00
62        0    -1.0000000000e+00  3.3333330000e-01  1.0000000000e+00
63        0     1.0000000000e+00 -1.0000000000e+00  3.3333330000e-01
64        0     1.0000000000e+00 -1.0000000000e+00  6.6666670000e-01
65        0     6.6666670000e-01 -1.0000000000e+00  1.0000000000e+00
66        0     3.3333330000e-01 -1.0000000000e+00  1.0000000000e+00
67        0     1.0000000000e+00 -6.6666670000e-01  1.0000000000e+00
68        0     1.0000000000e+00 -3.3333330000e-01  1.0000000000e+00
69        0     1.0000000000e+00 -5.5511150000e-17  1.0000000000e+00
70        0     1.0000000000e+00  3.3333330000e-01  1.0000000000e+00
71        0     1.0000000000e+00  6.6666670000e-01  1.0000000000e+00
72        0     6.6666670000e-01  1.0000000000e+00  1.0000000000e+00
73        0     3.3333330000e-01  1.0000000000e+00  1.0000000000e+00
74        0     1.0000000000e+00  1.0000000000e+00 -6.6666670000e-01
75        0     1.0000000000e+00  1.0000000000e+00 -3.3333330000e-01
76        0     1.0000000000e+00  1.0000000000e+00  6.6666670000e-01
77        0     1.0000000000e+00  1.0000000000e+00  3.3333330000e-01
78        0    -3.3333330000e-01  1.0000000000e+00  1.0000000000e+00
79        0    -6.6666670000e-01  1.0000000000e+00  1.0000000000e+00
80        0     0.0000000000e+00  6.6666670000e-01  1.0000000000e+00
81        0     0.0000000000e+00  3.3333330000e-01  1.0000000000e+00
82        0     0.0000000000e+00  1.0000000000e+00  6.6666670000e-01
83        0     0.0000000000e+00  1.0000000000e+00  3.3333330000e-01
84        0     3.3333330000e-01  1.0000000000e+00  0.0000000000e+00
85        0     6.6666670000e-01  1.0000000000e+00  0.0000000000e+00
86        0     1.0000000000e+00  6.6666670000e-01  0.0000000000e+00
87        0     1.0000000000e+00  3.3333330000e-01  0.0000000000e+00
88        0     3.7151020000e-01  0.0000000000e+00 -6.2848980000e-01
89        0     6.2848980000e-01  0.0000000000e+00 -3.7151020000e-01
90        0     0.0000000000e+00  3.7151020000e-01 -6.2848980000e-01
91        0     0.0000000000e+00  6.2848980000e-01 -3.7151020000e-01
92        0    -3.7151020000e-01  6.2848980000e-01  0.0000000000e+00
93        0    -6.2848980000e-01  3.7151020000e-01  0.0000000000e+00
94        0    -6.2848980000e-01  0.0000000000e+00  3.7151020000e-01
95        0    -3.7151020000e-01  0.0000000000e+00  6.2848980000e-01
96        0     0.0000000000e+00 -6.2848980000e-01  3.7151020000e-01
97        0     0.0000000000e+00 -3.7151020000e-01  6.2848980000e-01
98        0     6.2848980000e-01 -3.7151020000e-01  0.0000000000e+00
99        0     3.7151020000e-01 -6.2848980000e-01  0.0000000000e+00
100       0     3.7151020000e-01  6.2848980000e-01 -1.0000000000e+00
101       0     6.2848980000e-01  3.7151020000e-01 -1.0000000000e+00
102       0    -1.0000000000e+00  6.2848980000e-01  3.7151020000e-01
103       0    -1.0000000000e+00  3.7151020000e-01  6.2848980000e-01
104       0     3.7151020000e-01 -1.0000000000e+00  6.2848980000e-01
105       0     6.2848980000e-01 -1.0000000000e+00  3.7151020000e-01
106       0    -6.2848980000e-01  6.2848980000e-01  1.0000000000e+00
107       0    -3.7151020000e-01  3.7151020000e-01  1.0000000000e+00
108       0    -6.2848980000e-01  1.0000000000e+00  6.2848980000e-01
109       0    -3.7151020000e-01  1.0000000000e+00  3.7151020000e-01
110       0     1.0000000000e+00  6.2848980000e-01 -6.2848980000e-01
111       0     1.0000000000e+00  3.7151020000e-01 -3.7151020000e-01
112       0     3.4528950000e-01  7.1413680000e-01  1.0000000000e+00
113       0     3.2729110000e-01  3.2588100000e-01  1.0000000000e+00
114       0     6.7270890000e-01  4.3885710000e-01  1.0000000000e+00
115       0     3.2729110000e-01 -5.6488060000e-02  1.0000000000e+00
116       0     6.7270890000e-01  5.6488060000e-02  1.0000000000e+00
117       0     3.2729110000e-01 -4.3885710000e-01  1.0000000000e+00
118       0     6.7270890000e-01 -3.2588100000e-01  1.0000000000e+00
119       0     6.5471050000e-01 -7.1413680000e-01  1.0000000000e+00
120       0     6.2848980000e-01  1.0000000000e+00 -3.7151020000e-01
121       0     3.7151020000e-01  1.0000000000e+00 -6.2848980000e-01
122       0    -2.2747690000e-17  3.7151020000e-01  3.7151020000e-01
123       0    -3.8482630000e-17  6.2848980000e-01  6.2848980000e-01
124       0     6.2848980000e-01  1.0000000000e+00  6.2848980000e-01
125       0     3.7151020000e-01  1.0000000000e+00  3.7151020000e-01
126       0     1.0000000000e+00  7.1413680000e-01  3.4528950000e-01
127       0     1.0000000000e+00  3.2588100000e-01  3.2729110000e-01
128       0     1.0000000000e+00  4.3885710000e-01  6.7270890000e-01
129       0     1.0000000000e+00 -5.6488060000e-02  3.2729110000e-01
130       0     1.0000000000e+00  5.6488060000e-02  6.7270890000e-01
131       0     1.0000000000e+00 -4.3885710000e-01  3.2729110000e-01
132       0     1.0000000000e+00 -3.2588100000e-01  6.7270890000e-01
133       0     1.0000000000e+00 -7.1413680000e-01  6.5471050000e-01
134       0     3.7151020000e-01  3.7151020000e-01 -2.2747690000e-17
135       0     6.2848980000e-01  6.2848980000e-01 -3.8482630000e-17
136       0     6.7099540000e-01 -6.8668110000e-01  6.7099540000e-01
137       0     3.8141860000e-01 -4.6376380000e-01  3.8141860000e-01
138       0     6.1858140000e-01 -3.2381840000e-01  6.1858140000e-01
139       0     3.8141860000e-01 -6.9972660000e-02  3.8141860000e-01
140       0     6.1858140000e-01  6.9972660000e-02  6.1858140000e-01
141       0     3.8141860000e-01  3.2381840000e-01  3.8141860000e-01
142       0     6.1858140000e-01  4.6376380000e-01  6.1858140000e-01
143       0     3.2900460000e-01  6.8668110000e-01  3.2900460000e-01
144       0    -6.9524370000e-01  4.3826410000e-01  6.8122670000e-01
145       0    -4.3826410000e-01  6.9524370000e-01  6.8122670000e-01
146       0    -3.0475630000e-01  5.6173590000e-01  3.1877330000e-01
147       0    -5.6173590000e-01  3.0475630000e-01  3.1877330000e-01
148       0     6.8122670000e-01  4.3826410000e-01 -6.9524370000e-01
149       0     6.8122670000e-01  6.9524370000e-01 -4.3826410000e-01
150       0     3.0475630000e-01  3.1877330000e-01 -5.6173590000e-01
151       0     5.6173590000e-01  3.1877330000e-01 -3.0475630000e-01
];

simp=[
%-------- ---- ---- ------------------- ---------------------------------------
% Simp-ID Grp  Mat  Face-Types          Vertex-Numbers
%-------- ---- ---- ------------------- ---------------------------------------
0         0    0    0    0    0    0    19 32 141 134
1         0    0    0    0    0    0    31 32 134 143
2         0    0    0    0    0    0    31 134 135 143
3         0    0    0    1    0    0    2 31 84 83
4         0    0    0    0    0    0    31 84 143 135
5         0    0    0    0    0    1    84 85 125 135
6         0    0    0    0    0    1    86 87 127 135
7         0    0    0    0    0    0    87 134 141 135
8         0    0    0    0    0    0    20 87 134 141
9         0    0    0    0    0    0    19 20 134 141
10        0    0    0    0    0    1    126 127 128 142
11        0    0    0    0    0    1    127 129 130 140
12        0    0    0    0    1    0    127 128 140 130
13        0    0    0    0    0    1    129 131 132 138
14        0    0    0    0    1    0    129 130 138 132
15        0    0    0    0    1    0    131 132 136 133
16        0    0    0    0    1    1    18 77 85 126
17        0    0    0    1    0    0    18 85 86 126
18        0    0    0    1    0    0    1 20 127 87
19        0    0    0    1    0    0    1 20 129 127
20        0    0    0    1    0    1    1 20 54 129
21        0    0    0    1    0    0    54 98 131 129
22        0    0    0    0    1    1    53 54 98 131
23        0    0    0    0    0    1    53 63 131 105
24        0    0    0    0    1    0    63 131 136 133
25        0    0    0    0    0    1    63 64 133 136
26        0    0    0    0    1    0    13 64 136 133
27        0    0    0    0    0    1    13 67 133 136
28        0    0    0    0    0    1    67 132 133 136
29        0    0    0    0    1    1    67 68 132 118
30        0    0    0    0    1    1    68 69 132 118
31        0    0    0    1    0    0    69 118 130 132
32        0    0    0    0    1    1    69 70 130 116
33        0    0    0    0    0    1    70 128 130 140
34        0    0    0    0    0    1    70 71 128 142
35        0    0    0    0    0    1    71 76 128 142
36        0    0    0    1    1    1    16 71 72 76
37        0    0    0    0    0    1    76 126 128 142
38        0    0    0    0    1    1    76 77 126 124
39        0    0    0    0    1    1    72 73 124 112
40        0    0    0    1    0    0    73 112 125 124
41        0    0    0    0    0    1    82 83 125 123
42        0    0    0    0    0    1    83 84 125 143
43        0    0    0    0    1    0    112 113 142 114
44        0    0    0    0    1    0    113 115 140 116
45        0    0    0    0    0    1    113 114 116 140
46        0    0    0    0    1    0    115 117 138 118
47        0    0    0    0    0    1    115 116 118 140
48        0    0    0    0    0    1    117 118 119 136
49        0    0    0    0    1    1    17 73 112 82
50        0    0    0    0    1    0    17 80 82 112
51        0    0    0    0    0    1    80 112 113 123
52        0    0    0    0    1    0    80 81 123 113
53        0    0    0    1    0    0    3 43 81 113
54        0    0    0    1    0    0    3 43 113 115
55        0    0    0    1    0    0    45 97 115 117
56        0    0    0    0    1    1    45 46 97 117
57        0    0    0    0    1    0    46 66 104 117
58        0    0    0    0    0    1    66 117 119 136
59        0    0    0    0    0    1    65 66 119 136
60        0    0    0    0    0    1    13 65 119 136
61        0    0    0    0    1    0    13 67 136 119
62        0    0    0    0    1    0    67 118 136 119
63        0    0    0    0    1    0    70 114 140 116
64        0    0    0    0    1    0    70 71 142 114
65        0    0    0    0    1    0    71 72 142 114
66        0    0    0    0    1    0    72 112 142 114
67        0    0    0    1    1    1    5 48 51 49
68        0    0    0    0    0    1    48 49 96 51
69        0    0    0    0    1    1    49 50 96 99
70        0    0    0    0    1    0    50 96 137 97
71        0    0    0    1    1    0    0 19 44 50
72        0    0    0    0    1    0    44 50 137 97
73        0    0    0    0    1    0    43 44 139 97
74        0    0    0    0    0    1    46 96 97 137
75        0    0    0    1    1    1    10 46 66 47
76        0    0    0    0    0    1    46 47 96 104
77        0    0    0    0    1    1    47 48 96 104
78        0    0    0    0    1    0    19 20 139 98
79        0    0    0    0    1    0    19 98 137 99
80        0    0    0    0    0    1    19 50 99 137
81        0    0    0    0    1    1    51 52 99 105
82        0    0    0    0    0    1    52 98 99 105
83        0    0    0    1    1    1    9 52 63 53
84        0    0    0    0    0    1    52 53 98 105
85        0    0    0    0    1    0    63 64 136 105
86        0    0    0    0    1    0    64 104 136 105
87        0    0    0    0    0    1    13 64 65 136
88        0    0    0    0    1    0    64 65 136 104
89        0    0    0    0    1    0    65 66 136 104
90        0    0    0    0    0    0    32 44 141 122
91        0    0    0    0    0    0    43 44 122 141
92        0    0    0    0    0    0    43 122 123 141
93        0    0    0    0    0    0    83 122 143 123
94        0    0    0    0    0    0    31 83 122 143
95        0    0    0    0    0    0    31 32 143 122
96        0    0    0    0    0    0    112 113 123 142
97        0    0    0    0    0    0    0 19 32 44
98        0    0    0    0    0    0    125 135 142 143
99        0    0    0    0    0    0    31 83 143 84
100       0    0    0    0    0    0    20 127 141 140
101       0    0    0    0    0    0    20 129 140 139
102       0    0    0    0    0    0    105 131 137 136
103       0    0    0    0    0    0    63 105 136 131
104       0    0    0    0    1    0    20 54 129 98
105       0    0    0    0    0    0    20 127 140 129
106       0    0    0    0    0    0    127 135 141 142
107       0    0    0    0    0    0    129 138 140 139
108       0    0    0    0    0    0    20 98 129 139
109       0    0    0    0    0    0    98 129 139 138
110       0    0    0    0    0    0    131 132 138 136
111       0    0    0    0    0    0    66 104 117 136
112       0    0    0    0    0    0    43 44 141 139
113       0    0    0    0    0    0    87 127 135 141
114       0    0    0    0    0    0    118 130 138 140
115       0    0    0    0    0    0    98 129 138 131
116       0    0    0    1    0    0    52 53 105 63
117       0    0    0    0    0    0    117 118 136 138
118       0    0    0    0    0    0    98 131 138 137
119       0    0    0    0    0    0    118 130 132 138
120       0    0    0    0    0    0    70 128 140 142
121       0    0    0    0    0    0    112 123 125 142
122       0    0    0    0    0    0    43 113 141 123
123       0    0    0    0    0    0    72 112 124 142
124       0    0    0    0    0    0    112 124 142 125
125       0    0    0    0    0    0    80 82 112 123
126       0    0    0    0    0    1    73 82 125 112
127       0    0    0    0    0    0    116 118 140 130
128       0    0    0    0    0    0    124 126 142 135
129       0    0    0    0    0    0    85 124 135 126
130       0    0    0    0    0    0    76 124 126 142
131       0    0    0    0    0    0    19 44 137 139
132       0    0    0    0    0    0    19 44 139 141
133       0    0    0    0    0    0    43 115 139 140
134       0    0    0    0    0    0    71 72 124 142
135       0    0    0    1    1    0    3 43 115 45
136       0    0    0    0    0    0    97 115 117 138
137       0    0    0    0    0    0    97 115 138 139
138       0    0    0    0    0    0    46 104 137 117
139       0    0    1    0    0    0    46 47 104 66
140       0    0    1    0    0    0    71 72 76 124
141       0    0    0    0    0    0    71 76 142 124
142       0    0    0    0    0    0    67 118 132 136
143       0    0    0    0    1    0    69 116 130 118
144       0    0    0    0    0    0    70 114 142 140
145       0    0    0    0    0    0    32 134 143 141
146       0    0    0    0    0    1    49 51 99 96
147       0    0    0    0    1    0    48 51 96 105
148       0    0    0    0    0    0    20 87 141 127
149       0    0    0    0    0    0    113 123 142 141
150       0    0    0    0    0    0    19 44 50 137
151       0    0    0    0    0    1    43 45 97 115
152       0    0    0    0    0    0    43 97 139 115
153       0    0    0    0    0    0    46 97 117 137
154       0    0    0    1    0    0    48 96 104 105
155       0    0    0    0    0    0    131 136 138 137
156       0    0    0    0    0    0    53 98 105 131
157       0    0    0    0    0    0    98 99 105 137
158       0    0    0    0    0    0    51 96 105 99
159       0    0    0    0    0    0    98 105 131 137
160       0    0    0    0    0    0    46 96 137 104
161       0    0    0    0    0    0    104 105 137 136
162       0    0    0    0    0    0    19 20 141 139
163       0    0    0    0    0    0    96 99 137 105
164       0    0    0    0    0    1    85 124 125 135
165       0    0    0    0    0    0    43 81 113 123
166       0    0    0    0    0    0    85 86 126 135
167       0    0    0    0    1    0    86 126 135 127
168       0    0    0    0    0    0    19 32 44 141
169       0    0    0    0    1    0    77 85 126 124
170       0    0    0    0    0    0    82 112 123 125
171       0    0    0    0    0    0    134 135 143 141
172       0    0    0    0    0    0    127 128 142 140
173       0    0    0    0    0    0    50 96 99 137
174       0    0    0    0    0    0    97 137 139 138
175       0    0    0    0    0    0    124 125 135 142
176       0    0    0    0    0    0    43 139 141 140
177       0    0    0    0    0    0    104 117 136 137
178       0    0    0    0    0    0    70 116 140 130
179       0    0    0    0    0    0    129 130 140 138
180       0    0    0    0    0    0    126 127 142 135
181       0    0    0    0    0    0    96 104 105 137
182       0    0    0    0    0    0    43 113 140 141
183       0    0    0    0    0    0    118 132 136 138
184       0    0    0    0    0    0    98 137 138 139
185       0    0    0    0    0    0    43 113 115 140
186       0    0    0    0    0    0    115 118 138 140
187       0    0    0    0    0    0    97 117 137 138
188       0    0    0    0    0    0    117 136 137 138
189       0    0    0    0    0    0    44 97 137 139
190       0    0    0    0    0    0    19 98 139 137
191       0    0    0    0    0    0    113 114 140 142
192       0    0    0    0    0    0    122 123 141 143
193       0    0    0    0    0    0    115 138 139 140
194       0    0    0    0    0    0    20 139 140 141
195       0    0    0    0    0    0    32 122 141 143
196       0    0    0    0    0    0    127 140 142 141
197       0    0    0    0    0    0    113 140 141 142
198       0    0    0    0    0    0    123 141 143 142
199       0    0    0    0    0    0    135 141 142 143
200       0    0    0    0    0    0    123 125 142 143
201       0    0    0    0    0    0    83 123 143 125
202       0    0    0    0    0    0    84 125 143 135
203       0    0    0    0    0    0    32 44 122 146
204       0    0    0    0    0    1    43 44 95 122
205       0    0    0    0    0    0    43 95 123 122
206       0    0    0    0    0    1    80 81 107 123
207       0    0    0    0    0    0    80 82 123 145
208       0    0    0    0    0    0    82 83 123 146
209       0    0    0    0    0    0    83 122 123 146
210       0    0    0    0    0    0    31 83 146 122
211       0    0    0    0    0    0    31 32 122 146
212       0    0    0    0    0    1    33 34 92 146
213       0    0    0    0    1    0    34 92 146 93
214       0    0    0    1    1    1    12 34 59 35
215       0    0    0    0    0    1    34 35 93 102
216       0    0    0    0    1    1    35 36 93 102
217       0    0    0    1    1    1    4 36 39 37
218       0    0    0    0    0    1    36 37 93 147
219       0    0    0    0    0    1    37 38 93 147
220       0    0    0    0    0    1    38 92 93 146
221       0    0    0    1    0    1    0 32 38 44
222       0    0    0    0    1    0    32 38 146 92
223       0    0    0    0    1    0    31 32 146 92
224       0    0    0    1    1    0    2 31 83 33
225       0    0    0    0    0    1    31 33 92 146
226       0    0    0    0    1    1    39 40 94 103
227       0    0    0    0    1    0    40 94 144 95
228       0    0    0    1    1    1    8 40 62 41
229       0    0    0    0    0    1    40 41 95 144
230       0    0    0    0    1    1    41 42 95 107
231       0    0    0    0    1    1    3 42 81 43
232       0    0    0    0    0    1    42 43 95 81
233       0    0    0    0    0    1    44 94 95 147
234       0    0    0    0    1    0    38 44 147 94
235       0    0    0    0    1    0    37 38 147 94
236       0    0    0    0    0    1    37 39 94 147
237       0    0    0    0    1    1    59 60 102 108
238       0    0    0    0    1    0    60 102 144 103
239       0    0    0    1    1    1    15 60 79 61
240       0    0    0    0    0    1    60 61 103 144
241       0    0    0    0    0    1    61 62 103 144
242       0    0    0    0    1    0    40 62 144 103
243       0    0    0    0    0    1    39 102 103 147
244       0    0    0    0    1    0    36 39 147 102
245       0    0    0    0    0    1    61 79 106 108
246       0    0    0    0    1    0    61 62 144 106
247       0    0    0    0    0    1    62 106 107 144
248       0    0    0    0    0    1    41 62 107 144
249       0    0    0    0    1    0    80 106 145 107
250       0    0    0    0    1    1    17 78 82 80
251       0    0    0    0    0    1    78 80 106 145
252       0    0    0    0    1    1    78 79 108 106
253       0    0    0    0    1    0    33 34 146 109
254       0    0    0    0    0    1    33 83 109 146
255       0    0    0    0    1    0    82 83 146 109
256       0    0    0    0    0    1    82 108 109 146
257       0    0    0    0    1    0    78 82 145 108
258       0    0    0    0    0    0    32 38 44 146
259       0    0    0    0    0    0    82 108 146 145
260       0    0    0    0    0    0    95 122 146 123
261       0    0    0    0    0    0    81 95 107 123
262       0    0    0    0    0    0    80 107 145 123
263       0    0    0    0    0    0    78 80 145 82
264       0    0    0    0    0    0    78 106 108 145
265       0    0    0    0    0    0    31 33 146 83
266       0    0    0    0    0    0    93 102 147 146
267       0    0    1    0    0    0    34 35 102 59
268       0    0    0    0    0    0    36 93 102 147
269       0    0    0    0    0    0    34 93 146 102
270       0    0    0    0    0    0    36 37 147 39
271       0    0    0    0    0    0    94 103 144 147
272       0    0    0    0    0    0    40 94 103 144
273       0    0    0    0    0    0    39 94 147 103
274       0    0    0    0    0    0    108 144 146 145
275       0    0    0    0    0    0    94 95 147 144
276       0    0    0    0    0    0    40 41 144 62
277       0    0    0    0    0    0    95 107 123 146
278       0    0    0    0    0    0    41 95 144 107
279       0    0    0    0    1    0    42 81 95 107
280       0    0    0    0    0    0    107 123 146 145
281       0    0    0    0    0    0    95 107 146 144
282       0    0    0    0    0    0    102 103 147 144
283       0    0    0    1    0    0    60 61 108 79
284       0    0    0    0    1    0    34 59 102 109
285       0    0    0    0    0    0    60 61 144 108
286       0    0    0    0    0    0    102 108 146 109
287       0    0    0    0    0    0    106 107 144 145
288       0    0    0    0    0    0    102 108 144 146
289       0    0    0    0    0    0    61 106 144 108
290       0    0    0    0    0    0    60 102 108 144
291       0    0    0    1    0    0    59 102 109 108
292       0    0    0    0    0    0    34 102 146 109
293       0    0    0    0    0    0    82 123 145 146
294       0    0    0    0    0    0    43 81 123 95
295       0    0    0    0    0    0    44 95 122 146
296       0    0    0    0    0    0    107 144 145 146
297       0    0    0    0    0    0    106 108 145 144
298       0    0    0    0    0    0    102 144 147 146
299       0    0    0    0    0    0    38 44 146 147
300       0    0    0    0    0    0    38 93 147 146
301       0    0    0    0    0    0    95 144 146 147
302       0    0    0    0    0    0    44 95 146 147
303       0    0    0    1    0    1    0 19 26 32
304       0    0    0    0    0    0    19 26 32 134
305       0    0    0    0    0    1    31 32 91 134
306       0    0    0    0    0    0    31 91 135 134
307       0    0    0    0    0    0    31 84 135 91
308       0    0    0    0    0    0    84 85 135 149
309       0    0    0    0    0    0    85 86 135 149
310       0    0    0    0    1    0    86 87 135 111
311       0    0    0    0    0    0    87 134 135 151
312       0    0    0    1    1    0    1 20 87 21
313       0    0    0    0    0    0    20 87 151 134
314       0    0    0    0    0    0    19 20 151 134
315       0    0    0    0    1    1    18 75 86 85
316       0    0    0    0    0    1    75 85 120 149
317       0    0    0    0    1    0    84 85 149 120
318       0    0    0    0    1    0    84 120 149 121
319       0    0    0    0    1    1    2 30 84 31
320       0    0    0    0    1    0    30 84 91 121
321       0    0    0    0    1    1    29 30 91 121
322       0    0    0    0    1    1    57 58 100 121
323       0    0    0    0    0    1    57 120 121 149
324       0    0    0    0    0    1    57 74 120 149
325       0    0    0    0    0    1    74 75 120 149
326       0    0    0    0    0    1    27 28 90 150
327       0    0    0    0    1    0    28 90 150 91
328       0    0    0    1    1    1    11 28 58 29
329       0    0    0    0    0    1    32 90 91 150
330       0    0    0    0    1    0    26 32 150 90
331       0    0    0    0    1    0    25 26 150 90
332       0    0    0    0    0    1    25 27 90 150
333       0    0    0    0    1    0    27 28 150 100
334       0    0    0    0    1    0    27 100 150 101
335       0    0    0    1    1    1    6 24 27 25
336       0    0    0    0    1    0    24 27 150 101
337       0    0    0    0    1    1    23 24 88 101
338       0    0    0    0    0    1    23 55 101 148
339       0    0    0    0    0    1    55 56 101 148
340       0    0    0    0    0    1    56 100 101 148
341       0    0    0    1    1    1    14 56 74 57
342       0    0    0    0    0    1    56 57 100 149
343       0    0    0    0    1    0    55 56 148 110
344       0    0    0    0    1    0    55 110 148 111
345       0    0    0    1    1    1    7 22 55 23
346       0    0    0    0    1    0    22 55 148 111
347       0    0    0    0    1    1    21 22 89 111
348       0    0    0    0    0    1    86 110 111 149
349       0    0    0    0    1    0    75 86 149 110
350       0    0    0    0    1    0    74 75 149 110
351       0    0    0    0    1    0    56 74 149 110
352       0    0    0    0    0    1    24 25 88 150
353       0    0    0    0    0    1    25 26 88 150
354       0    0    0    0    1    0    26 88 150 89
355       0    0    0    0    1    0    19 26 150 89
356       0    0    0    0    0    1    19 20 89 151
357       0    0    0    0    0    1    22 88 89 148
358       0    0    0    0    0    1    22 23 88 148
359       0    0    0    0    0    0    19 26 134 150
360       0    0    0    0    0    0    26 32 134 150
361       0    0    0    0    0    0    100 148 149 150
362       0    0    0    0    0    1    30 31 91 84
363       0    0    0    0    0    0    84 91 149 135
364       0    0    0    0    0    0    91 121 149 150
365       0    0    0    0    0    0    75 85 149 86
366       0    0    0    0    0    0    86 111 135 149
367       0    0    0    0    0    0    20 89 151 111
368       0    0    0    0    0    0    84 91 121 149
369       0    0    0    0    1    0    28 58 121 100
370       0    0    1    0    0    0    28 29 121 58
371       0    0    0    0    0    1    28 29 91 121
372       0    0    0    0    0    0    134 149 151 150
373       0    0    0    0    0    0    56 57 149 74
374       0    0    0    0    0    0    111 148 151 149
375       0    0    0    0    0    0    28 91 150 121
376       0    0    0    0    0    0    91 134 149 135
377       0    0    0    0    0    0    24 25 150 27
378       0    0    0    0    0    0    28 100 121 150
379       0    0    0    0    0    0    57 100 149 121
380       0    0    0    0    0    0    100 101 148 150
381       0    0    0    0    0    0    23 88 148 101
382       0    0    0    0    0    0    22 23 148 55
383       0    0    0    0    0    0    111 135 149 151
384       0    0    0    0    0    0    56 100 148 149
385       0    0    0    0    0    0    56 110 149 148
386       0    0    0    0    0    0    110 111 149 148
387       0    0    1    0    0    0    20 21 111 87
388       0    0    0    0    0    0    20 87 111 151
389       0    0    0    0    0    0    134 135 151 149
390       0    0    0    0    0    0    88 89 148 150
391       0    0    0    0    0    0    19 134 151 150
392       0    0    0    0    0    0    22 89 111 148
393       0    0    0    0    0    0    32 91 134 150
394       0    0    0    0    0    0    87 111 151 135
395       0    0    0    0    0    0    19 89 150 151
396       0    0    0    0    0    1    20 21 89 111
397       0    0    0    0    0    0    100 121 150 149
398       0    0    0    0    0    0    24 88 101 150
399       0    0    0    0    0    0    88 101 150 148
400       0    0    0    0    0    0    89 111 148 151
401       0    0    0    0    0    0    89 148 150 151
402       0    0    0    0    0    0    148 149 150 151
403       0    0    0    0    0    0    91 134 150 149
];

mcsf_end=1;

