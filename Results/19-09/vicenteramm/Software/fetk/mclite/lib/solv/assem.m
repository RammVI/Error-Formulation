function [A,F]=assem(assKey, parm, U, Ug, VERT, SIMP);
%ASSEM  Finite element nonlinear residual and tangent system assembly
%
% Usage: [A,F]=assem(assKey, parm, U, Ug, VERT, SIMP);
%
% Input:
%
%    assKey = assembly path
%    parm   = pde parameters
%    U      = zeros(N*n,1); current finite element solution
%    Ug     = zeros(N*n,1); essential/dirichlet boundary function
%    VERT   = zeros(N,8);  <see "read.m" for description of the datastructure>
%    SIMP   = zeros(L,17); <see "read.m" for description of the datastructure>
%
%       NOTE: n=number of PDE eqns and number of unknowns per spatial point.
%             This is determined by the length of U compared to the length
%             of the vertex array VERT; the length of U is always a positive
%             integer multiple of the length of VERT, because MCLite uses
%             the same (linear) element for each component of an elliptic
%             system.
%
% Output:
%
%    assKey==2 ==> A = zeros(N*n,N*n) = nonlinear tangent stiffness matrix
%                  F = zeros(N*n,1)   = nonlinear residual vector
%
%    assKey==3 ==> A = zeros(N*n,N*n) = (not constructed -- just a zero matrix)
%                  F = zeros(N*n,1)   = nonlinear residual vector
%
%    assKey==4 ==> A = zeros(N*n,N*n) = mass matrix
%                  F = zeros(N*n,1)   = (not constructed -- just a zero vector)
%
% Other Variables:
%
%    T  = number of vertices in a single element
%    Q  = number of quadrature points in the interior of the master element
%    QE = number of quadrature points along the master element edges
%
%    INTPT   = zeros(Q,2)  = master element volume integration points
%    W       = zeros(Q,1)  = master element volume integration weights
%    PHI     = zeros(T,Q)  = values of local basis funcs at integ pts
%    PHIX    = zeros(T,Q)  = values of x-derivs of local basis funcs integ pts
%    PHIY    = zeros(T,Q)  = values of y-derivs of local basis funcs integ pts
%
%    INTPTE0 = zeros(QE,2) = type (1,2) EDGE integration points
%    WE0     = zeros(QE,1) = type (1,2) EDGE integration weights
%    PHIE0   = zeros(T,QE) = type (1,2) basis function values at EDGE integ pts
%
%    INTPTE1 = zeros(QE,2) = type (2,0) EDGE integration points
%    WE1     = zeros(QE,1) = type (2,0) EDGE integration weights
%    PHIE1   = zeros(T,QE) = type (2,0) basis function values at EDGE integ pts
%
%    INTPTE2 = zeros(QE,2) = type (0,1) EDGE integration points
%    WE2     = zeros(QE,1) = type (0,1) EDGE integration weights
%    PHIE2   = zeros(T,QE) = type (0,1) basis function values at EDGE integ pts
%
% Author:   Michael Holst
% rcsid="$Id: assem.m,v 1.1.1.1 2007/04/27 08:28:19 hrg Exp $"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   %%% first recover various problem dimensions.
   [N,eight]     = size(VERT);
   [L,seventeen] = size(SIMP);
   T=3;

   %%% get master element info.
   [INTPT,  W,  PHI,  PHIX,  PHIY,   ...
    INTPTE0,WE0,PHIE0,PHIXE0,PHIYE0, ...
    INTPTE1,WE1,PHIE1,PHIXE1,PHIYE1, ...
    INTPTE2,WE2,PHIE2,PHIXE2,PHIYE2] = mastr('P1',1);
   [Q,two]  = size(INTPT);
   [QE,two] = size(INTPTE0);

   %%% determine how many unknowns per spatial point.
   [Nn,one] = size(U);
   n  = Nn / N;
   n2 = n * n;

   %%% initialize global stiffness matrix and global load vector.
   A = spalloc(Nn,Nn,0);  
   F = zeros(Nn,1);

   %%% get coordinates of vertices of each element.
   Xi0(1:L,1) = VERT(SIMP(1:L,1),1);
   Xi1(1:L,1) = VERT(SIMP(1:L,2),1);
   Xi2(1:L,1) = VERT(SIMP(1:L,3),1);
   Yi0(1:L,1) = VERT(SIMP(1:L,1),2);
   Yi1(1:L,1) = VERT(SIMP(1:L,2),2);
   Yi2(1:L,1) = VERT(SIMP(1:L,3),2);

   %%% build triangle affine transformations (from master el to arbitrary).
   f00(1:L,1) = Xi1(1:L,1) - Xi0(1:L,1);
   f01(1:L,1) = Xi2(1:L,1) - Xi0(1:L,1);
   f10(1:L,1) = Yi1(1:L,1) - Yi0(1:L,1);
   f11(1:L,1) = Yi2(1:L,1) - Yi0(1:L,1);
   b0(1:L,1)  = Xi0(1:L,1);
   b1(1:L,1)  = Yi0(1:L,1);

   %%% build the triangle jacobians.
   D (1:L,1) = abs( f00 .* f11 - f01 .* f10 );
   Di(1:L,1) = 1.0 ./ D(1:L,1);

   %%% build the edge jacobians.
   DE = zeros(L,3);
   DE(1:L,1) = (1/sqrt(2)) ...
             * sqrt( (Xi2(1:L,1) - Xi1(1:L,1)) .* (Xi2(1:L,1) - Xi1(1:L,1)) ...
             + (Yi2(1:L,1) - Yi1(1:L,1)) .* (Yi2(1:L,1) - Yi1(1:L,1)) );
   DE(1:L,2) = sqrt( (Xi0(1:L,1) - Xi2(1:L,1)) .* (Xi0(1:L,1) - Xi2(1:L,1)) ...
             + (Yi0(1:L,1) - Yi2(1:L,1)) .* (Yi0(1:L,1) - Yi2(1:L,1)) );
   DE(1:L,3) = sqrt( (Xi0(1:L,1) - Xi1(1:L,1)) .* (Xi0(1:L,1) - Xi1(1:L,1)) ...
             + (Yi0(1:L,1) - Yi1(1:L,1)) .* (Yi0(1:L,1) - Yi1(1:L,1)) );

   %%% build inverse triangle transformations (from arbitrary el to master).
   g00(1:L,1) =  Di(1:L,1) .* f11(1:L,1);
   g01(1:L,1) = -Di(1:L,1) .* f01(1:L,1);
   g10(1:L,1) = -Di(1:L,1) .* f10(1:L,1);
   g11(1:L,1) =  Di(1:L,1) .* f00(1:L,1);

   %%% simplex and face types
   stype = SIMP(1:L,11);
   ftype = SIMP(1:L,5:7);

   %%% at each r loop iter, compute el stiffness matrix and el load vec entries
   for r=1:T

      %%% initialize element load entry fr
      fr(1:L,1:n) = zeros(L,n);

      %%% at each s loop iter, compute el stiffness matrix entries
      for s=1:T

         %%% initialize element matrix entry ars 
         ars(1:L,1:n2) = zeros(L,n2);

         %%% go thru volume quad pts m=1:Q for contrib to element matrix/vector
         for m=1:Q

            %%% quadrature weights and transformation jacobians
            WM = zeros(L,1);
            WM(1:L,1) = W(m) * D(1:L,1);

            %%% get quad pts by mapping master Gauss pts to our element
            xm(1:L,1)=f00(1:L,1)*INTPT(m,1)+f01(1:L,1)*INTPT(m,2)+b0(1:L,1);
            ym(1:L,1)=f10(1:L,1)*INTPT(m,1)+f11(1:L,1)*INTPT(m,2)+b1(1:L,1);

            %%% evaluate function at the current quadrature point m
            for q=1:T
               p(q) = PHI(q,m);
               px(1:L,q) = g00(1:L,1) * PHIX(q,m) + g10(1:L,1) * PHIY(q,m);
               py(1:L,q) = g01(1:L,1) * PHIX(q,m) + g11(1:L,1) * PHIY(q,m);
            end

            %%% current solution and diri func at quad points in all simps
            UU   = zeros(L,n);
            UUX  = zeros(L,n);
            UUY  = zeros(L,n);
            UUg  = zeros(L,n);
            UUgX = zeros(L,n);
            UUgY = zeros(L,n);
            for i=1:n
                Nii = N * (i-1);
                UU(1:L,i)   = p(1) * U( SIMP(1:L,1)+Nii, 1 ) ...
                            + p(2) * U( SIMP(1:L,2)+Nii, 1 ) ...
                            + p(3) * U( SIMP(1:L,3)+Nii, 1 );
                UUX(1:L,i)  = px(1:L,1) .* U( SIMP(1:L,1)+Nii, 1 ) ...
                            + px(1:L,2) .* U( SIMP(1:L,2)+Nii, 1 ) ...
                            + px(1:L,3) .* U( SIMP(1:L,3)+Nii, 1 );
                UUY(1:L,i)  = py(1:L,1) .* U( SIMP(1:L,1)+Nii, 1 ) ...
                            + py(1:L,2) .* U( SIMP(1:L,2)+Nii, 1 ) ...
                            + py(1:L,3) .* U( SIMP(1:L,3)+Nii, 1 );
                UUg(1:L,i)  = p(1) * Ug( SIMP(1:L,1)+Nii, 1 ) ...
                            + p(2) * Ug( SIMP(1:L,2)+Nii, 1 ) ...
                            + p(3) * Ug( SIMP(1:L,3)+Nii, 1 );
                UUgX(1:L,i) = px(1:L,1) .* Ug( SIMP(1:L,1)+Nii, 1 ) ...
                            + px(1:L,2) .* Ug( SIMP(1:L,2)+Nii, 1 ) ...
                            + px(1:L,3) .* Ug( SIMP(1:L,3)+Nii, 1 );
                UUgY(1:L,i) = py(1:L,1) .* Ug( SIMP(1:L,1)+Nii, 1 ) ...
                            + py(1:L,2) .* Ug( SIMP(1:L,2)+Nii, 1 ) ...
                            + py(1:L,3) .* Ug( SIMP(1:L,3)+Nii, 1 );
            end

            %%% do the load vector while we are here...
            if ((assKey == 2) | (assKey == 3))
               if (s==1)

                  %%% evaluate nonlinear function at quadrature points
                  evalKey = 0;
                  thetam = fu_v(evalKey,parm,L, xm,ym,stype, ...
                     (UU+UUg), (UUX+UUgX), (UUY+UUgY), ...
                     p(r), px(1:L,r),py(1:L,r));

                  %%% do quad using weights W
                  for ii=1:n
                     fr(1:L,ii) = fr(1:L,ii) + WM(1:L,1) .* thetam(1:L,ii);
                  end

               end;
            else
               %%% do nothing
            end;

            %%% evaluate bilinear form at quadrature points
            if (assKey == 2)

               %%% build element stiffness matrix
               evalKey = 0;
               thetam = dfu_wv(evalKey,parm,L, xm,ym,stype, ...
                               (UU+UUg), (UUX+UUgX), (UUY+UUgY), ...
                               p(s),px(1:L,s),py(1:L,s), ...
                               p(r),px(1:L,r),py(1:L,r));

               %%% do quad using weights W
               for ii=1:n2
                  ars(1:L,ii) = ars(1:L,ii) + WM(1:L,1) .* thetam(1:L,ii);
               end

            elseif (assKey == 4)

               %%% build element mass matrix
               thetam = zeros(L,n2);
               for ii=1:n2
                  thetam(1:L,ii) = rho(xm,ym) * p(s) * p(r);
               end

               %%% do quad using weights W
               for ii=1:n2
                  ars(1:L,ii) = ars(1:L,ii) + WM(1:L,1) .* thetam(1:L,ii);
               end

            else
               %%% do nothing
            end;

         end %%% the "m=1:Q" loop

         %%% go thru edge quad pts m=1:QE for contrib to element matrix/vector
         for m=1:QE
            
            %%% handle contributions to each of the possible edges
            for edge=1:T

               %%% use the appropriate quadrature info for particular edge
               if (edge==1)
                  INTE = INTPTE0;
                  PHIE = PHIE0;
                  WE   = WE0;
               elseif (edge==2)
                  INTE = INTPTE1;
                  PHIE = PHIE1;
                  WE   = WE1;
               elseif (edge==3)
                  INTE = INTPTE2;
                  PHIE = PHIE2;
                  WE   = WE2;
               end

               %%% quadrature weights and transformation jacobians
               %%% (mask indicating which edges are neumann edges)
               ftype_edge = ftype(1:L,edge);
               maskit  = vneum(ftype_edge);
               WEM = zeros(L,1);
               WEM(1:L,1) = WE(m) * DE(1:L,edge) .* maskit;

               %%% Get quad pts by mapping master Gauss pts to our element
               xm(1:L,1)=f00(1:L,1)*INTE(m,1)+f01(1:L,1)*INTE(m,2)+b0(1:L,1);
               ym(1:L,1)=f10(1:L,1)*INTE(m,1)+f11(1:L,1)*INTE(m,2)+b1(1:L,1);

               %%% evaluate function at the current quadrature point m
               for q=1:T
                  p(q) = PHIE(q,m);
                  px(1:L,q) = zeros(L,1);
                  py(1:L,q) = zeros(L,1);
               end

               %%% current solution and diri func at quad points in all simps
               UU   = zeros(L,n);
               UUX  = zeros(L,n);
               UUY  = zeros(L,n);
               UUg  = zeros(L,n);
               UUgX = zeros(L,n);
               UUgY = zeros(L,n);
               for ii=1:n
                   Nii = N * (ii-1);
                   UU(1:L,ii)   = p(1) * U( SIMP(1:L,1)+Nii, 1 ) ...
                                + p(2) * U( SIMP(1:L,2)+Nii, 1 ) ...
                                + p(3) * U( SIMP(1:L,3)+Nii, 1 );
                   UUX(1:L,ii)  = px(1:L,1) .* U( SIMP(1:L,1)+Nii, 1 ) ...
                                + px(1:L,2) .* U( SIMP(1:L,2)+Nii, 1 ) ...
                                + px(1:L,3) .* U( SIMP(1:L,3)+Nii, 1 );
                   UUY(1:L,ii)  = py(1:L,1) .* U( SIMP(1:L,1)+Nii, 1 ) ...
                                + py(1:L,2) .* U( SIMP(1:L,2)+Nii, 1 ) ...
                                + py(1:L,3) .* U( SIMP(1:L,3)+Nii, 1 );
                   UUg(1:L,ii)  = p(1) * Ug( SIMP(1:L,1)+Nii, 1 ) ...
                                + p(2) * Ug( SIMP(1:L,2)+Nii, 1 ) ...
                                + p(3) * Ug( SIMP(1:L,3)+Nii, 1 );
                   UUgX(1:L,ii) = px(1:L,1) .* Ug( SIMP(1:L,1)+Nii, 1 ) ...
                                + px(1:L,2) .* Ug( SIMP(1:L,2)+Nii, 1 ) ...
                                + px(1:L,3) .* Ug( SIMP(1:L,3)+Nii, 1 );
                   UUgY(1:L,ii) = py(1:L,1) .* Ug( SIMP(1:L,1)+Nii, 1 ) ...
                                + py(1:L,2) .* Ug( SIMP(1:L,2)+Nii, 1 ) ...
                                + py(1:L,3) .* Ug( SIMP(1:L,3)+Nii, 1 );
               end

               %%% neumann/edge term (do the load vector while here...)
               if ((assKey == 2) | (assKey == 3))
                  if (s==1)        
                     evalKey = 1;
                     thetam = fu_v(evalKey,parm,L, xm,ym,ftype_edge, ...
                        (UU+UUg), (UUX+UUgX), (UUY+UUgY), ...
                        p(r), px(1:L,r),py(1:L,r));
                     for ii=1:n
                        fr(1:L,ii) = fr(1:L,ii) + WEM(1:L,1) .* thetam(1:L,ii);
                     end
                  end
               else
                  %%% do nothing
               end;

               %%% robin/volume term (evaluate bilinear form at quad point)
               if (assKey == 2)
                  evalKey = 1;
                  thetam = dfu_wv(evalKey,parm,L, xm,ym,ftype_edge, ...
                               (UU+UUg), (UUX+UUgX), (UUY+UUgY), ...
                               p(s),px(1:L,s),py(1:L,s), ...
                               p(r),px(1:L,r),py(1:L,r));
                  for ii=1:n2
                     ars(1:L,ii) = ars(1:L,ii) + WEM(1:L,1) .* thetam(1:L,ii);
                  end
               else
                  %%% do nothing
               end;
            end %%% the "edge=1:T" loop
         end %%% the "m=1:QE" loop

         %%% Add contribution to global stiffness matrix.
         i(1:L,1)  = SIMP(1:L,r);
         j(1:L,1)  = SIMP(1:L,s);
         ik(1:L,1) = VERT( SIMP(1:L,r), 6 );
         jk(1:L,1) = VERT( SIMP(1:L,s), 6 );
         for l=1:L
            if ( (~vdiri(ik(l))) & (~vdiri(jk(l))) )
               idx=0;
               for ii=1:n
                  Nii = N * (ii-1);
                  for jj=1:n
                     Njj = N * (jj-1);
                     idx=idx+1; 
                     A(i(l)+Nii,j(l)+Njj) = A(i(l)+Nii,j(l)+Njj) + ars(l,idx);
                  end 
               end 
            else
               if ( vdiri(ik(l)) )
                  for ii=1:n
                     Nii = N * (ii-1);
                     A(i(l)+Nii,i(l)+Nii) = 1.0;
                  end
               end
               if ( vdiri(jk(l)) )
                  for jj=1:n
                     Njj = N * (jj-1);
                     A(j(l)+Njj,j(l)+Njj) = 1.0;
                  end
               end
            end
         end
      end %%% the "s=1:T" loop

      %%% Add contribution to global load vector
      i(1:L,1)  = SIMP(1:L,r);
      ik(1:L,1) = VERT( SIMP(1:L,r), 6 );
      for l=1:L
         if ( ~vdiri(ik(l)) )
            for ii=1:n
               Nii = N * (ii-1);
               F(i(l)+Nii) = F(i(l)+Nii) - fr(l,ii);
            end
         else
            for ii=1:n
               Nii = N * (ii-1);
               F(i(l)+Nii) = 0.;
            end
         end
      end
   end %%% the "r=1:T" loop

