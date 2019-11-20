function [INTPT,  W,  PHI,  PHIX,  PHIY,   ...
          INTPTE0,WE0,PHIE0,PHIXE0,PHIYE0, ...
          INTPTE1,WE1,PHIE1,PHIXE1,PHIYE1, ...
          INTPTE2,WE2,PHIE2,PHIXE2,PHIYE2] = P2(Q);
%MASTR  Define the master element and master element basis function info
%
% Usage: [INTPT,  W,  PHI,  PHIX,  PHIY,   ...
%         INTPTE0,WE0,PHIE0,PHIXE0,PHIYE0, ...
%         INTPTE1,WE1,PHIE1,PHIXE1,PHIYE1, ...
%         INTPTE2,WE2,PHIE2,PHIXE2,PHIYE2] = P2(Q);
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
% Purpose:
%
%    Define the master element and master element basis function info.
%
%        PHI3(0,1)=1 |\
%                    | \
%                    |  \
%      PHI6(0,1/2)=1 .   . PHI5(1/2,1/2)=1
%                    |    \
%        PHI1(0,0)=1 |__.__\ PHI2(1,0)=1
%                    PHI4(1/2,0)=1
%
% Author:   D. Estep, M. Holst, D. Mikulencak
% rcsid="$Id: P2.m,v 1.1.1.1 2007/04/27 08:28:08 hrg Exp $"

%%% initialize

   T  = 6;   %%% number of degrees of freedom (quadratic)
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

   if Q==1     % O(h^2)     
      INTPT(1:Q,1:2) = [1/3 1/3];
      W(1:Q,1)     = 1/2;
   elseif Q==3 % O(h^3)  
      INTPT(1:Q,1:2) = [0 1/2; 1/2 0; 1/2 1/2];
      W(1:Q,1) = 1/2*[1/3; 1/3; 1/3];
   elseif Q==4 % O(h^4)
      INTPT(1:Q,1:2) = [1/3 1/3; 1/5 1/5; 3/5 1/5; 1/5 3/5];
      W(1:Q,1) = 1/2*[-27/48; 25/48; 25/48; 25/48];
   elseif Q==7
      INTPT(1:Q,1:2) = [1/3 1/3;                            ...
                        (6+sqrt(15))/21 (6+sqrt(15))/21 ;   ...
                        (9-2*sqrt(15))/21 (6+sqrt(15))/21;  ...
                        (6+sqrt(15))/21 (9-2*sqrt(15))/21;  ...
                        (6-sqrt(15))/21 (6-sqrt(15))/21;    ...
                        (9+2*sqrt(15))/21 (6-sqrt(15))/21;  ...
                        (6-sqrt(15))/21 (9+2*sqrt(15))/21];
      W(1:Q,1) = [9/80;               ...
                 (155+sqrt(15))/2400; ...
                 (155+sqrt(15))/2400; ...
                 (155+sqrt(15))/2400; ...
                 (155-sqrt(15))/2400; ...
                 (155-sqrt(15))/2400; ...
                 (155-sqrt(15))/2400]; 
   end

%%% evaluate master element basis functions

   for k=1:Q
      x = INTPT(k,1);
      y = INTPT(k,2);
      L1 = 1 - x - y;
      L2 = x;
      L3 = y;
      PHI(1,k)  = L1*(2*L1 - 1);
      PHI(2,k)  = L2*(2*L2 - 1);
      PHI(3,k)  = L3*(2*L3 - 1);
      PHI(4,k)  = 4*L1*L2;
      PHI(5,k)  = 4*L2*L3;
      PHI(6,k)  = 4*L1*L3;
      
      PHIX(1,k) = 4*x + 4*y - 3;
      PHIX(2,k) = 4*x - 1;
      PHIX(3,k) = 0;
      PHIX(4,k) = -8*x - 4*y + 4;
      PHIX(5,k) = 4*y;
      PHIX(6,k) = -4*y;
   
      PHIY(1,k) = 4*x + 4*y - 3;
      PHIY(2,k) = 0;
      PHIY(3,k) = 4*y - 1;
      PHIY(4,k) = -4*x;
      PHIY(5,k) = 4*x;
      PHIY(6,k) = -8*y - 4*x + 4;
   end

%%% define edge integration pts, weights, and master element basis functions

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
   L1E0=1-x0-y0;
   L1E1=1-x1-y1;
   L1E2=1-x2-y2;
   L2E0=x0;
   L2E1=x1;
   L2E2=x2;
   L3E0=y0;
   L3E1=y1;
   L3E2=y2;

   for k=1:QE
      PHIE0(1,k)  = L1E0(k)*(2*L1E0(k) - 1);
      PHIE0(2,k)  = L2E0(k)*(2*L2E0(k) - 1);
      PHIE0(3,k)  = L3E0(k)*(2*L3E0(k) - 1);
      PHIE0(4,k)  = 4*L1E0(k)*L2E0(k);
      PHIE0(5,k)  = 4*L2E0(k)*L3E0(k);
      PHIE0(6,k)  = 4*L1E0(k)*L3E0(k);
   
      PHIE1(1,k)  = L1E1(k)*(2*L1E1(k) - 1);
      PHIE1(2,k)  = L2E1(k)*(2*L2E1(k) - 1);
      PHIE1(3,k)  = L3E1(k)*(2*L3E1(k) - 1);
      PHIE1(4,k)  = 4*L1E1(k)*L2E1(k);
      PHIE1(5,k)  = 4*L2E1(k)*L3E1(k);
      PHIE1(6,k)  = 4*L1E1(k)*L3E1(k);
   
      PHIE2(1,k)  = L1E2(k)*(2*L1E2(k) - 1);
      PHIE2(2,k)  = L2E2(k)*(2*L2E2(k) - 1);
      PHIE2(3,k)  = L3E2(k)*(2*L3E2(k) - 1);
      PHIE2(4,k)  = 4*L1E2(k)*L2E2(k);
      PHIE2(5,k)  = 4*L2E2(k)*L3E2(k);
      PHIE2(6,k)  = 4*L1E2(k)*L3E2(k);
   end

