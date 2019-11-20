function [INTPT,  W,  PHI,  PHIX,  PHIY,   ...
          INTPTE0,WE0,PHIE0,PHIXE0,PHIYE0, ...
          INTPTE1,WE1,PHIE1,PHIXE1,PHIYE1, ...
          INTPTE2,WE2,PHIE2,PHIXE2,PHIYE2] = P1bub(Q);
%MASTR  Define the master element and master element basis function info
%
% Usage: [INTPT,  W,  PHI,  PHIX,  PHIY,   ...
%         INTPTE0,WE0,PHIE0,PHIXE0,PHIYE0, ...
%         INTPTE1,WE1,PHIE1,PHIXE1,PHIYE1, ...
%         INTPTE2,WE2,PHIE2,PHIXE2,PHIYE2] = P1bub(Q);
%
% Input:
%
%    (Q  = number of quadrature points in the interior of the master element)
%
% Output:
%
%    (T  = number of nodes in a single element)
%    (QE = number of quadrature points along the master element edges)
%
%    INTPT   = zeros(Q,2)  = master element volume integration points
%    W       = zeros(Q,1)  = master element volume integration weights
%    PHI     = zeros(T,Q)  = values of local basis funcs at integ pts
%    PHIX    = zeros(T,Q)  = values of x-derivs of local basis funcs integ pts
%    PHIY    = zeros(T,Q)  = values of y-derivs of local basis funcs integ pts
%
%    INTPTE0 = zeros(QE,2) = type (1,2) EDGE integration points
%    WE0     = zeros(QE,1) = type (1,2) EDGE integration weights
%    PHIE0   = zeros(T,QE) = type (1,2) basis func values at EDGE integ pts
%    PHIXE0  = zeros(T,QE) = type (1,2) basis func x-derivs at EDGE integ pts
%    PHIYE0  = zeros(T,QE) = type (1,2) basis func y-derivs at EDGE integ pts
%
%    INTPTE1 = zeros(QE,2) = type (2,0) EDGE integration points
%    WE1     = zeros(QE,1) = type (2,0) EDGE integration weights
%    PHIE1   = zeros(T,QE) = type (2,0) basis function values at EDGE integ pts
%    PHIXE1  = zeros(T,QE) = type (2,0) basis func x-derivs at EDGE integ pts
%    PHIYE1  = zeros(T,QE) = type (2,0) basis func y-derivs at EDGE integ pts
%
%    INTPTE2 = zeros(QE,2) = type (0,1) EDGE integration points
%    WE2     = zeros(QE,1) = type (0,1) EDGE integration weights
%    PHIE2   = zeros(T,QE) = type (0,1) basis function values at EDGE integ pts
%    PHIXE2  = zeros(T,QE) = type (0,1) basis func x-derivs at EDGE integ pts
%    PHIYE2  = zeros(T,QE) = type (0,1) basis func y-derivs at EDGE integ pts
%
% Notes:
%
%    Define the master element and master element basis function info.
%
%                   PHI3(0,1)=1 |\
%                               | \
%                               |  \
%                               | . \ PHI4(1/3,1/3)=1 
%                               |    \
%                   PHI1(0,0)=1 |_____\ PHI2(1,0)=1
%
% Author:   D. Estep, M. Holst, D. Mikulencak
% rcsid="$Id: P1bub.m,v 1.1.1.1 2007/04/27 08:28:08 hrg Exp $"

%%% initialize

   T  = 3+1; %%% number of degrees of freedom (linear + quadratic bubble)
   QE = 2;   %%% number of edge integration points

   INTPT   = zeros(Q,2);
   W       = zeros(Q,1);
   PHI     = zeros(T,Q);
   PHIX    = zeros(T,Q);
   PHIY    = zeros(T,Q);

   INTPTE0 = zeros(QE,2);
   WE0     = zeros(QE,1);
   PHIE0   = zeros(T,QE);
   PHIXE0  = zeros(T,QE);
   PHIYE0  = zeros(T,QE);

   INTPTE1 = zeros(QE,2);
   WE1     = zeros(QE,1);
   PHIE1   = zeros(T,QE);
   PHIXE1  = zeros(T,QE);
   PHIYE1  = zeros(T,QE);

   INTPTE2 = zeros(QE,2);
   WE2     = zeros(QE,1);
   PHIE2   = zeros(T,QE);
   PHIXE2  = zeros(T,QE);
   PHIYE2  = zeros(T,QE);


%%% define volume integration points and weights

   if Q==1      % O(h^2)  % problems with stiffness matrix being singular
      INTPT(1:Q,1:2) = [1/3 1/3];
      W(1:Q,1)     = 1/2;
   elseif Q==3  % O(h^3)  % best results
      INTPT(1:Q,1:2) = [0 1/2; 1/2 0; 1/2 1/2];
      W(1:Q,1) = 1/2*[1/3; 1/3; 1/3];
   elseif Q==4  % O(h^4)  % problems with stiffness matrix being singular
      INTPT(1:Q,1:2) = [1/3 1/3; 1/5 1/5; 3/5 1/5; 1/5 3/5];
      W(1:Q,1) = 1/2*[-27/48; 25/48; 25/48; 25/48];
   elseif Q==7
      INTPT(1:Q,1:2) = [1/3 1/3;                           ...
                        (6+sqrt(15))/21 (6+sqrt(15))/21 ;  ...
                        (9-2*sqrt(15))/21 (6+sqrt(15))/21; ...
                        (6+sqrt(15))/21 (9-2*sqrt(15))/21; ...
                        (6-sqrt(15))/21 (6-sqrt(15))/21;   ...
                        (9+2*sqrt(15))/21 (6-sqrt(15))/21; ...
                        (6-sqrt(15))/21 (9+2*sqrt(15))/21];
      W(1:Q,1) = [9/80;                ...
                  (155+sqrt(15))/2400; ...
                  (155+sqrt(15))/2400; ...
                  (155+sqrt(15))/2400; ...
                  (155-sqrt(15))/2400; ...
                  (155-sqrt(15))/2400; ...
                  (155-sqrt(15))/2400]; 
   end

%%% Evaluate master element basis functions

   for k=1:Q
      x = INTPT(k,1);
      y = INTPT(k,2);
      PHI(1,k)  = 1 - x - y;
      PHI(2,k)  = x;
      PHI(3,k)  = y;
      PHI(4,k)  = 27*PHI(1,k)*PHI(2,k)*PHI(3,k);
      PHIX(1,k) = -1;
      PHIX(2,k) = 1;
      PHIX(3,k) = 0;
      PHIX(4,k) = 27*PHI(3,k)*(PHIX(1,k)*PHI(2,k) + PHI(1,k)*PHIX(2,k));
      PHIY(1,k) = -1;
      PHIY(2,k) = 0;
      PHIY(3,k) = 1;
      PHIY(4,k) = 27*PHI(2,k)*(PHIY(1,k)*PHI(3,k) + PHI(1,k)*PHIY(3,k));
   end

%%% define edge integration pts, weights, and master element basis functions
%%% only difference from P1 is that must add zeros to PHIE for row
%%% 4 since bubble vanishes on edges

   c1 = 1/2;
   c2 = 1/(2*sqrt(3));

   INTPTE0(1:QE,1:2) = [ (1-c1+c2) (c1-c2) ; (1-c1-c2) (c1+c2) ];
   WE0(1:QE,1:1)     = [ c1 ; c1 ];

   INTPTE1(1:QE,1:2) = [ 0 (c1+c2) ; 0 (c1-c2) ];
   WE1(1:QE,1:1)     = [ c1 ; c1 ];

   INTPTE2(1:QE,1:2) = [ (c1-c2) 0 ; (c1+c2) 0 ];
   WE2(1:QE,1:1)     = [ c1 ; c1 ];

   x0=INTPTE0(:,1);
   x1=INTPTE1(:,1);
   x2=INTPTE2(:,1);
   y0=INTPTE0(:,2);
   y1=INTPTE1(:,2);
   y2=INTPTE2(:,2);

   for k=1:QE
      PHIE0(1,k) = 1 - x0(k) - y0(k);
      PHIE0(2,k) = x0(k);
      PHIE0(3,k) = y0(k);
      PHIE0(4,k) = 27*PHIE0(1,k)*PHIE0(2,k)*PHIE0(3,k);
   
      PHIE1(1,k) = 1 - x1(k) - y1(k);
      PHIE1(2,k) = x1(k);
      PHIE1(3,k) = y1(k);
      PHIE1(4,k) = 27*PHIE1(1,k)*PHIE1(2,k)*PHIE1(3,k);
   
      PHIE2(1,k) = 1 - x2(k) - y2(k);
      PHIE2(2,k) = x2(k);
      PHIE2(3,k) = y2(k);
      PHIE2(4,k) = 27*PHIE2(1,k)*PHIE2(2,k)*PHIE2(3,k);
   end

