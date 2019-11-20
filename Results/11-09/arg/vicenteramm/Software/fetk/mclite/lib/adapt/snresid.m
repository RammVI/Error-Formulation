function [indic] = snresid(VERT,SIMP,U,Ug)
%RESID  Calculate the strong nonlinear residual indicator in all elements
%
% Usage: [indic] = snresid(VERT,SIMP,U,Ug);
%
% Input:
%
%    U      = current solution
%    Ug     = dirichlet function
%    VERT   = <see "read.m" for description of the datastructure>
%    SIMP   = <see "read.m" for description of the datastructure>
%
% Output:
%
%    indic  = the value of the indicator for each element
%
% Notes:
%
%    This routine computes the strong nonlinear residual indicator
%    for all elements in the mesh.
%
% Author:   Michael Holst
% rcsid="$Id: snresid.m,v 1.1.1.1 2007/04/27 08:28:06 hrg Exp $"

%%% first recover various problem dimensions.
[N,eight]     = size(VERT);
[L,seventeen] = size(SIMP);
T=3;

%%% determine how many unknowns per spatial point.
[Nn,one] = size(U);
n  = Nn / N;
n2 = n * n;

%%% some parameters 
VSMALL = 1.0e-7;
parm = zeros(10,1);

%%% get master element info (use high-order quadrature!).
deg = 2;
[INTPT,  W,  PHI,  PHIX,  PHIY,   ...
 INTPTE0,WE0,PHIE0,PHIXE0,PHIYE0, ...
 INTPTE1,WE1,PHIE1,PHIXE1,PHIYE1, ...
 INTPTE2,WE2,PHIE2,PHIXE2,PHIYE2] = mastr('P1',1);
[Q,two]  = size(INTPT);
[QE,two] = size(INTPTE0);

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

%%% THE BIG LOOP OVER THE ELEMENTS
indic = zeros(L,1);
for el=1:L

    %%% initialize the error estimate for this element
    errEst = 0.0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% compute the indicator contribution from this VOLUME
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% initialize contribution from this volume
    errEstS = 0.0;

    %%% Cycle thru quad points
    face = -1;
    for m=1:Q

        %%% quadrature weights and transformation jacobians
        WM = W(m) * D(el);

        %%% get quad pts by mapping master Gauss pts to our element
        xm=f00(el)*INTPT(m,1)+f01(el)*INTPT(m,2)+b0(el);
        ym=f10(el)*INTPT(m,1)+f11(el)*INTPT(m,2)+b1(el);

        %%% evaluate function at the current quadrature point m
        p  = zeros(T,1);
        px = zeros(T,1);
        py = zeros(T,1);
        for q=1:T
            p(q) = PHI(q,m);
            px(q) = g00(el) * PHIX(q,m) + g10(el) * PHIY(q,m);
            py(q) = g01(el) * PHIX(q,m) + g11(el) * PHIY(q,m);
        end

        %%% current solution and diri func at quad points in all simps
        UU   = zeros(n,1);
        UUX  = zeros(n,1);
        UUY  = zeros(n,1);
        UUg  = zeros(n,1);
        UUgX = zeros(n,1);
        UUgY = zeros(n,1);
        for i=1:n
            Nii = N * (i-1);
            UU(i)   = p(1) * U( SIMP(el,1)+Nii, 1 ) ...
                    + p(2) * U( SIMP(el,2)+Nii, 1 ) ...
                    + p(3) * U( SIMP(el,3)+Nii, 1 );
            UUX(i)  = px(1) .* U( SIMP(el,1)+Nii, 1 ) ...
                    + px(2) .* U( SIMP(el,2)+Nii, 1 ) ...
                    + px(3) .* U( SIMP(el,3)+Nii, 1 );
            UUY(i)  = py(1) .* U( SIMP(el,1)+Nii, 1 ) ...
                    + py(2) .* U( SIMP(el,2)+Nii, 1 ) ...
                    + py(3) .* U( SIMP(el,3)+Nii, 1 );
            UUg(i)  = p(1) * Ug( SIMP(el,1)+Nii, 1 ) ...
                    + p(2) * Ug( SIMP(el,2)+Nii, 1 ) ...
                    + p(3) * Ug( SIMP(el,3)+Nii, 1 );
            UUgX(i) = px(1) .* Ug( SIMP(el,1)+Nii, 1 ) ...
                    + px(2) .* Ug( SIMP(el,2)+Nii, 1 ) ...
                    + px(3) .* Ug( SIMP(el,3)+Nii, 1 );
            UUgY(i) = py(1) .* Ug( SIMP(el,1)+Nii, 1 ) ...
                    + py(2) .* Ug( SIMP(el,2)+Nii, 1 ) ...
                    + py(3) .* Ug( SIMP(el,3)+Nii, 1 );
        end

        %%% get strong residual in element:  b(u)^i - a(u)^{iq}_{~;q}
        evalKey = 0;
        type    = 0;
        nvec    = [ 0 , 0 ]; %% dummy unit normal...
        F       = fu( evalKey,parm, xm,ym,type,nvec, ...
                      UU+UUg, UUX+UUgX, UUY+UUgY );
        nval = 0.0;
        for i=1:n
            tmp = F(i);
            nval = nval + tmp*tmp;
        end
        errEstS = errEstS + ( WM * nval );

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% compute the indicator contribution from all FACES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% initialize contribution from this face
    errEstFi = 0.0;
    errEstFb = 0.0;

    %%% Cycle thru faces
    for face=1:T

        %%% ftype mask indicating which faces are interior faces
        ftype_face = ftype(el,face);
        maskit = vint(ftype_face);

        %%% get face-nabor across face, his local face number, and trans
        [el_N,face_N] = getnabor(VERT,SIMP,el,face);

        %%% deal with boundary case
        if (el_N == 0)
            assertw( (face_N == 0), 'snresid 1' );
            assertw( (maskit == 0), 'snresid 2' );
            el_N   = 1;
            face_N = 1;
        end

        %%% calculate unit normal vector to this face
        if (face==1) v1=2; v2=3; end;
        if (face==2) v1=3; v2=1; end;
        if (face==3) v1=2; v2=1; end;
        g1 = SIMP(el,v1);
        g2 = SIMP(el,v2);
        w1 = zeros(2,1);
        w2 = zeros(2,1);
        w1(1) = VERT(g1,1);
        w1(2) = VERT(g1,2);
        w2(1) = VERT(g2,1);
        w2(2) = VERT(g2,2);
        w3 = w1 - w2;
        len = sqrt( w3' * w3 );
        assertw( (len ~= 0), 'snresid 3' );
        nvec = zeros(1,2);
        nvec(1) = -w3(2) / len;
        nvec(2) =  w3(1) / len;

        %%% use the appropriate quadrature info for particular face
        if (face==1)
            INTE  = INTPTE0;
            WE    = WE0;
            PHIE  = PHIE0;  
            PHIXE = PHIXE0;  
            PHIYE = PHIYE0;  
        elseif (face==2)
            INTE = INTPTE1;
            WE    = WE1;    
            PHIE  = PHIE1;  
            PHIXE = PHIXE1;  
            PHIYE = PHIYE1;  
        elseif (face==3)
            INTE  = INTPTE2;
            WE    = WE2;    
            PHIE  = PHIE2;  
            PHIXE = PHIXE2;  
            PHIYE = PHIYE2;  
        end             

        %%% use the appropriate quadrature info for particular face
        if (face_N==1)
            INTE_N  = INTPTE0;
            WE_N    = WE0;
            PHIE_N  = PHIE0;  
            PHIXE_N = PHIXE0;  
            PHIYE_N = PHIYE0;  
        elseif (face_N==2)
            INTE_N  = INTPTE1;
            WE_N    = WE1;    
            PHIE_N  = PHIE1;  
            PHIXE_N = PHIXE1;  
            PHIYE_N = PHIYE1;  
        elseif (face_N==3)
            INTE_N  = INTPTE2;
            WE_N    = WE2;    
            PHIE_N  = PHIE2;  
            PHIXE_N = PHIXE2;  
            PHIYE_N = PHIYE2;  
        end             

        %%% Cycle thru quad points (diff face integs poss)
        for m=1:QE

            %%% the face nabor's "numbering" of quad pts is reversed!
            m_N = QE - m + 1;

            %%% quadrature weights and transformation jacobians
            WEM   = WE(m)   * DE(el,face);
            WEM_N = WE(m_N) * DE(el,face_N);

            %%% Get quad pts by mapping master Gauss pts to our element
            xm   = f00(el)*INTE(m,1) ...
                 + f01(el)*INTE(m,2)       + b0(el);
            ym   = f10(el)*INTE(m,1) ...
                 + f11(el)*INTE(m,2)       + b1(el);
            xm_N = f00(el_N)*INTE_N(m_N,1) ...
                 + f01(el_N)*INTE_N(m_N,2) + b0(el_N);
            ym_N = f10(el_N)*INTE_N(m_N,1) ...
                 + f11(el_N)*INTE_N(m_N,2) + b1(el_N);

            %%% TEST: see if we have the quad points lined up!
            for i=1:2
                assertw( (maskit*abs(xm - xm_N)<VSMALL), 'snresid 4' );
                assertw( (maskit*abs(ym - ym_N)<VSMALL), 'snresid 5' );
            end

            %%% evaluate function at the current quadrature point m
            p    = zeros(T,1);
            px   = zeros(T,1);
            py   = zeros(T,1);
            p_N  = zeros(T,1);
            px_N = zeros(T,1);
            py_N = zeros(T,1);
            for q=1:T
                p(q)    = PHIE(q,m);
                px(q)   = g00(el) * PHIXE(q,m) + g10(el) * PHIYE(q,m);
                py(q)   = g01(el) * PHIXE(q,m) + g11(el) * PHIYE(q,m);

                p_N(q)  = PHIE_N(q,m_N);
                px_N(q) = g00(el) * PHIXE_N(q,m) + g10(el) * PHIYE_N(q,m);
                py_N(q) = g01(el) * PHIXE_N(q,m) + g11(el) * PHIYE_N(q,m);
            end

            %%% current solution and diri func at quad points in all simps
            UU     = zeros(n,1);
            UUX    = zeros(n,1);                                            
            UUY    = zeros(n,1);
            UUg    = zeros(n,1);
            UUgX   = zeros(n,1);
            UUgY   = zeros(n,1);
            UU_N   = zeros(n,1);
            UUX_N  = zeros(n,1);                                            
            UUY_N  = zeros(n,1);
            UUg_N  = zeros(n,1);
            UUgX_N = zeros(n,1);
            UUgY_N = zeros(n,1);
            for ii=1:n
                Nii = N * (ii-1);

                UU(ii)   = p(1) * U( SIMP(el,1)+Nii, 1 ) ...
                         + p(2) * U( SIMP(el,2)+Nii, 1 ) ...
                         + p(3) * U( SIMP(el,3)+Nii, 1 );
                UUX(ii)  = px(1) .* U( SIMP(el,1)+Nii, 1 ) ...
                         + px(2) .* U( SIMP(el,2)+Nii, 1 ) ...
                         + px(3) .* U( SIMP(el,3)+Nii, 1 );
                UUY(ii)  = py(1) .* U( SIMP(el,1)+Nii, 1 ) ...
                         + py(2) .* U( SIMP(el,2)+Nii, 1 ) ...
                         + py(3) .* U( SIMP(el,3)+Nii, 1 );   
                UUg(ii)  = p(1) * Ug( SIMP(el,1)+Nii, 1 ) ...     
                         + p(2) * Ug( SIMP(el,2)+Nii, 1 ) ...  
                         + p(3) * Ug( SIMP(el,3)+Nii, 1 );   
                UUgX(ii) = px(1) .* Ug( SIMP(el,1)+Nii, 1 ) ...
                         + px(2) .* Ug( SIMP(el,2)+Nii, 1 ) ...
                         + px(3) .* Ug( SIMP(el,3)+Nii, 1 );   
                UUgY(ii) = py(1) .* Ug( SIMP(el,1)+Nii, 1 ) ...
                         + py(2) .* Ug( SIMP(el,2)+Nii, 1 ) ...
                         + py(3) .* Ug( SIMP(el,3)+Nii, 1 );   

                UU_N(ii)   = p(1) * U( SIMP(el_N,1)+Nii, 1 ) ...
                           + p(2) * U( SIMP(el_N,2)+Nii, 1 ) ...
                           + p(3) * U( SIMP(el_N,3)+Nii, 1 );
                UUX_N(ii)  = px(1) .* U( SIMP(el_N,1)+Nii, 1 ) ...
                           + px(2) .* U( SIMP(el_N,2)+Nii, 1 ) ...
                           + px(3) .* U( SIMP(el_N,3)+Nii, 1 );
                UUY_N(ii)  = py(1) .* U( SIMP(el_N,1)+Nii, 1 ) ...
                           + py(2) .* U( SIMP(el_N,2)+Nii, 1 ) ...
                           + py(3) .* U( SIMP(el_N,3)+Nii, 1 );   
                UUg_N(ii)  = p(1) * Ug( SIMP(el_N,1)+Nii, 1 ) ...     
                           + p(2) * Ug( SIMP(el_N,2)+Nii, 1 ) ...  
                           + p(3) * Ug( SIMP(el_N,3)+Nii, 1 );   
                UUgX_N(ii) = px(1) .* Ug( SIMP(el_N,1)+Nii, 1 ) ...
                           + px(2) .* Ug( SIMP(el_N,2)+Nii, 1 ) ...
                           + px(3) .* Ug( SIMP(el_N,3)+Nii, 1 );   
                UUgY_N(ii) = py(1) .* Ug( SIMP(el_N,1)+Nii, 1 ) ...
                           + py(2) .* Ug( SIMP(el_N,2)+Nii, 1 ) ...
                           + py(3) .* Ug( SIMP(el_N,3)+Nii, 1 );   
            end

            %%% get side1 of jump on interior face:  a(u1)^{iq} n_q
            evalKey = 2;
            type    = 0;
            F       = fu( evalKey,parm,xm,ym,type,nvec, ...
                          UU+UUg,UUX+UUgX,UUY+UUgY );

            %%% get side2 of jump on interior face:  a(u2)^{iq} n_q
            evalKey = 2;
            type    = 0;
            F_N     = fu( evalKey,parm,xm,ym,type,nvec, ...
                          UU_N+UUg_N,UUX_N+UUgX_N,UUY_N+UUgY_N );

            %%% jump term  [a(u1)^{iq} - a(u2)^{iq}] n_q
            nval = 0.0;
            for i=1:n
                tmp = F(i) - F_N(i);
                nval = nval + tmp*tmp;
            end
            errEstFi = errEstFi + ( maskit * WEM * nval );

            %%% if neumann boundary face, penalize using neumann residual
            %%% neumann condition residual:  c(u)^i + a(u)^{iq} n_q
            evalKey = 1;
            type    = ftype_face;
            F       = fu( evalKey,parm,xm,ym,type,nvec, ...
                          UU+UUg,UUX+UUgX,UUY+UUgY );
            nval = 0.0;
            for i=1:n
                tmp = F(i);
                nval = nval + tmp*tmp;
            end;
            errEstFb = errEstFb + ( (1-maskit) * WEM * nval );

        end
    end

   %%% take square root as final result over the element
   errEst = sqrt( errEstS  * D(el) * D(el) ...
                + errEstFi * 0.5 * DE(el,face) ...
                + errEstFb * DE(el,face) ...
                );
   indic(el) = errEst;
   %%% fprintf('element=<%d>   indicator=<%g>\n', el, indic(el));
end

return;

