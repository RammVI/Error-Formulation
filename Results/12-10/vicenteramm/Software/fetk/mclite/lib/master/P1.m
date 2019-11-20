function [INTPT,  W,  PHI,  PHIX,  PHIY,   ...
          INTPTE0,WE0,PHIE0,PHIXE0,PHIYE0, ...
          INTPTE1,WE1,PHIE1,PHIXE1,PHIYE1, ...
          INTPTE2,WE2,PHIE2,PHIXE2,PHIYE2] = P1(Q);
%MASTR  Define the master element and master element basis function info
%
% Usage: [INTPT,  W,  PHI,  PHIX,  PHIY,   ...
%         INTPTE0,WE0,PHIE0,PHIXE0,PHIYE0, ...
%         INTPTE1,WE1,PHIE1,PHIXE1,PHIYE1, ...
%         INTPTE2,WE2,PHIE2,PHIXE2,PHIYE2] = P1(Q);
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
%                  phi2(0,1)=1 |\
%                              | \
%                              |  \
%                              |   \
%                              |    \
%                  phi0(0,0)=1 |_____\ phi1(1,0)=1
%
%    \phi0(x,y) = 1 - x - y
%    \phi1(x,y) = x
%    \phi2(x,y) = y
%  
%    \phi0_x(x,y) = -1
%    \phi0_y(x,y) = -1
%  
%    \phi1_x(x,y) =  1
%    \phi1_y(x,y) =  0
%  
%    \phi2_x(x,y) =  0
%    \phi2_y(x,y) =  1
%
% Author:   Michael Holst
% rcsid="$Id: P1.m,v 1.1.1.1 2007/04/27 08:28:08 hrg Exp $"

%%% initialize

   T  = 3;   %%% triangles
   QE = 2;   %%% two edge integration points

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

   if Q==1      % O(h^2)    
      INTPT(1:Q,1:2) = [ 1./3. 1./3. ];
      W(1:Q,1)       = 1./2.;
   elseif Q==3  % O(h^3)  
      INTPT(1:Q,1:2) = [ 0. 1./2.; 1./2. 0.; 1./2. 1./2. ];
      W(1:Q,1)       = (1./2.)*[ 1./3.; 1./3.; 1./3. ];
   elseif Q==4  % O(h^4)
      INTPT(1:Q,1:2) = [ 1./3. 1./3.; 1./5. 1./5.; 3./5. 1./5.; 1./5. 3./5. ];
      W(1:Q,1)       = (1./2.)*[ -27./48.; 25./48.; 25./48.; 25./48. ];
   elseif Q==7
      INTPT(1:Q,1:2) = [ 1./3. 1./3.;                                   ...
                         (6.+sqrt(15.))/21. (6.+sqrt(15.))/21. ;        ...
                         (9.-2.*sqrt(15.)) / 21. (6.+sqrt(15.))/21.;    ...
                         (6.+sqrt(15.))    / 21. (9.-2.*sqrt(15.))/21.; ...
                         (6.-sqrt(15.))    / 21. (6.-sqrt(15.))/21.;    ...
                         (9.+2.*sqrt(15.)) / 21. (6.-sqrt(15.))/21.;    ...
                         (6.-sqrt(15.))    / 21. (9.+2.*sqrt(15.))/21.  ];
      W(1:Q,1)       = [ 9./80.;                 ...
                         (155.+sqrt(15.))/2400.; ...
                         (155.+sqrt(15.))/2400.; ...
                         (155.+sqrt(15.))/2400.; ...
                         (155.-sqrt(15.))/2400.; ...
                         (155.-sqrt(15.))/2400.; ...
                         (155.-sqrt(15.))/2400.  ]; 
   end

%%% evaluate master element basis functions at interior quadrature points

   for k=1:Q

      x = INTPT(k,1);
      y = INTPT(k,2);

      PHI(1,k)  = 1. - x - y;
      PHI(2,k)  = x;
      PHI(3,k)  = y;

      PHIX(1,k) = -1.;
      PHIX(2,k) =  1.;
      PHIX(3,k) =  0.;

      PHIY(1,k) = -1.;
      PHIY(2,k) =  0.;
      PHIY(3,k) =  1.;
   end

%%% define edge integration pts, weights, and master element basis functions

   c1 = 1./2.;
   c2 = 1./(2.*sqrt(3.));
   c3 = 1./sqrt(2.);

   INTPTE0(1:QE,1:2) = [ (1.-c1+c2) , (c1-c2) ; (1.-c1-c2) , (c1+c2) ];
   WE0(1:QE,1:1)     = [ c3 ; c3 ];
   PHIE0(1:T,1:QE)   = [ 0. , 0. ; (c1+c2) , (c1-c2) ; (c1-c2) , (c1+c2) ];
   PHIXE0(1:T,1:QE)  = [ -1. , -1. ; 1. , 1. ; 0. , 0. ];
   PHIYE0(1:T,1:QE)  = [ -1. , -1. ; 0. , 0. ; 1. , 1. ];

   INTPTE1(1:QE,1:2) = [ 0. , (c1+c2) ; 0. , (c1-c2) ];
   WE1(1:QE,1:1)     = [ c1 ; c1 ];
   PHIE1(1:T,1:QE)   = [ (c1-c2) , (c1+c2) ; 0. , 0. ; (c1+c2) , (c1-c2) ];
   PHIXE1(1:T,1:QE)  = [ -1. , -1. ; 1. , 1. ; 0. , 0. ];
   PHIYE1(1:T,1:QE)  = [ -1. , -1. ; 0. , 0. ; 1. , 1. ];

   INTPTE2(1:QE,1:2) = [ (c1-c2) , 0. ; (c1+c2) , 0. ];
   WE2(1:QE,1:1)     = [ c1 ; c1 ];
   PHIE2(1:T,1:QE)   = [ (c1+c2) , (c1-c2) ; (c1-c2) , (c1+c2) ; 0. , 0. ];
   PHIXE2(1:T,1:QE)  = [ -1. , -1. ; 1. , 1. ; 0. , 0. ];
   PHIYE2(1:T,1:QE)  = [ -1. , -1. ; 0. , 0. ; 1. , 1. ];

   x0=INTPTE0(:,1);
   x1=INTPTE1(:,1);
   x2=INTPTE2(:,1);
   y0=INTPTE0(:,2);
   y1=INTPTE1(:,2);
   y2=INTPTE2(:,2);

   % error check!
   chk_PHIE0 = zeros(T,QE);
   chk_PHIE1 = zeros(T,QE);
   chk_PHIE2 = zeros(T,QE);
   for k=1:QE
      chk_PHIE0(1,k)  = 1. - x0(k) - y0(k);
      chk_PHIE0(2,k)  = x0(k);
      chk_PHIE0(3,k)  = y0(k);
   
      chk_PHIE1(1,k) = 1. - x1(k) - y1(k);
      chk_PHIE1(2,k) = x1(k);
      chk_PHIE1(3,k) = y1(k);
   
      chk_PHIE2(1,k) = 1. - x2(k) - y2(k);
      chk_PHIE2(2,k) = x2(k);
      chk_PHIE2(3,k) = y2(k);
   end
   assertw( norm((PHIE0 - chk_PHIE0),'fro') < 1.0e-7, 'pt 1' );
   assertw( norm((PHIE1 - chk_PHIE1),'fro') < 1.0e-7, 'pt 1' );
   assertw( norm((PHIE2 - chk_PHIE2),'fro') < 1.0e-7, 'pt 1' );

