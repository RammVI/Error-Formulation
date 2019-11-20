function [VERT_F,SIMP_F,EDGE_F,QUE_F] = closure(VERT,SIMP,EDGE,QUE,whichq);
%CLOSURE  One-shot closure of a red-refined mesh via quad-/tri-/bi-section
%
% Usage: [VERT_F,SIMP_F,EDGE_F,QUE_F] = closure(VERT,SIMP,EDGE,QUE,whichq);
%
% Input:
%
%    whichq = which que we are currently using
%    VERT   = <see "read.m" for description of the datastructure>
%    SIMP   = <see "read.m" for description of the datastructure>
%    EDGE   = <see "read.m" for description of the datastructure>
%    QUE    = list of simplices to be refined
% 
% Output:
%    
%    VERT_F = list of vertices in the refined mesh
%    SIMP_F = conformally refined set of simplices
%    EDGE_F = list of edges created during refinement
%    QUE_F  = new list of simplices to be refined (EMPTY!)
%
% Author:   Michael Holst
% rcsid="$Id: closure.m,v 1.1.1.1 2007/04/27 08:28:05 hrg Exp $"

%%% recover array dimensions
[N,eight]     = size(VERT);
[L,seventeen] = size(SIMP); 
[K,eight]     = size(EDGE);
[P,one]       = size(QUE);
T=3;
    
%%% the refinement loop
KK = K; 
NN = N;
LL = L;
if (whichq==1)
    swapq = 2;
else
    swapq = 1;
end 
assertw( ((whichq==1)&(swapq==2)) | ((whichq==2)&(swapq==1)), 'closure 1' );
    
%%% loop over the elements
NEWQUE = [];
for refGuy = 1:P

    %%% get my name and generation number
    l=QUE(refGuy);
    mark = SIMP(l,18);
    type = SIMP(l,19);
    gen  = SIMP(l,20);

    %%% child mark, type, generation
    cmark = mark;
    ctype = 0;
    cgen  = gen+1;

    %%% check que tag status
    refq = SIMP(l,whichq+15);
    assertw( (refq ~= 0), 'closure 2' );

    %%% clear our que tags
    SIMP(l,whichq+15) = 0;
    SIMP(l,swapq+15)  = 0;

    %%% get our three nabor simplex numbers before killing the rings
    nabors    = zeros(3,1);
    [nabors(1),cf1] = getnabor(VERT,SIMP,l,1);  %%% nabor opposite vertex 1
    [nabors(2),cf2] = getnabor(VERT,SIMP,l,2);  %%% nabor opposite vertex 2
    [nabors(3),cf3] = getnabor(VERT,SIMP,l,3);  %%% nabor opposite vertex 3

    %%% remove this simplex from rings (will reuse as first child)
    [VERT,SIMP] = delsring(VERT,SIMP,l);

    %%% face/vertex info for the T vertices/faces forming this simplex
    vn = zeros(T);      %%% vertex numbers of this simp
    ve = zeros(T);      %%% first edge of each vertices' edge list
    v  = zeros(T,3);    %%% the x,y,z coordinates of the vertices
    ft = zeros(3,1);    %%% the 3 integer face types
    nb = zeros(3,1);    %%% the 3 nabors
    for i=1:T
        vn(i)  = SIMP( l, i );
        ft(i)  = SIMP( l, i+4 );
        ve(i)  = VERT( vn(i), 8 );
        v(i,:) = VERT( vn(i), 1:3 );
        nb(i)  = nabors( i );
    end

    %%% get existing vertices
    vn_NEW = zeros(3,3);
    for i=1:T
        for j=i+1:T
            edg = ve(i);
            found = 0;
            while ((found == 0) & (edg ~= 0))
                ev = EDGE(edg,1:3);
                if ( ((ev(1) == vn(i)) & (ev(2) == vn(j)) ) ...
                   | ((ev(1) == vn(j)) & (ev(2) == vn(i)) ) )
                    found = 1;
                else
                    ii=-1;
                    if ( ev(1) == vn(i) ) ii = 1; end;
                    if ( ev(2) == vn(i) ) ii = 2; end;
                    assertw( (ii >= 0), 'closure 3' );
                    edg = EDGE( edg, ii+6);
                end
            end

            if (found)
                vn_NEW(i,j) = ev(3);
            end
        end
    end
    c = vn_NEW;

    %%% The three bisection cases
    if ( c(1,2) & ~c(1,3) & ~c(2,3))

        %%% should only bisect a quadrasection child!
        assertw( (type>=1), 'closure 4' );

        %%% child simplex 1; create/add to rings (reuse parent structure)
        %%% (bisection child noted via the "0" in snew slot 19)
        %%% (generation of child noted in snew slot 20)
        LC=zeros(4,1);
        LC(1)=l;
        snew = [ vn(1),vn_NEW(1,2),vn(3),0, 0,ft(2),ft(3),0, ...
                 LC(1),0,0, 0,0,0,0, 0,0,cmark,ctype,cgen ];
        SIMP(LC(1),:) = snew; 
        [VERT,SIMP] = addsring(VERT,SIMP,LC(1));

        %%% child simplex 2; create/add to rings (new simplex location)
        %%% (bisection child noted via the "0" in snew slot 19)
        %%% (generation of child noted in snew slot 20)
        LL=LL+1;
        LC(2)=LL;
        snew = [ vn(2),vn(3),vn_NEW(1,2),0, 0,ft(3),ft(1),0, ...
                 LC(2),0,0, 0,0,0,0, 0,0,cmark,ctype,cgen ]; 
        SIMP = [ SIMP ; snew ];
        [VERT,SIMP] = addsring(VERT,SIMP,LC(2));

    elseif ( ~c(1,2) & c(1,3) & ~c(2,3))

        %%% should only bisect a quadrasection child!
        assertw( (type>=1), 'closure 5' );

        %%% child simplex 1; create/add to rings (reuse parent structure)
        %%% (bisection child noted via the "0" in snew slot 19)
        %%% (generation of child noted in snew slot 20)
        LC=zeros(4,1);
        LC(1)=l;
        snew = [ vn(1),vn(2),vn_NEW(1,3),0, 0,ft(2),ft(3),0, ...
                 LC(1),0,0, 0,0,0,0, 0,0,cmark,ctype,cgen ];
        SIMP(LC(1),:) = snew; 
        [VERT,SIMP] = addsring(VERT,SIMP,LC(1));

        %%% child simplex 2; create/add to rings (new simplex location)
        %%% (bisection child noted via the "0" in snew slot 19)
        %%% (generation of child noted in snew slot 20)
        LL=LL+1;
        LC(2)=LL;
        snew = [ vn(3),vn_NEW(1,3),vn(2),0, 0,ft(1),ft(2),0, ...
                 LC(2),0,0, 0,0,0,0, 0,0,cmark,ctype,cgen ]; 
        SIMP = [ SIMP ; snew ];
        [VERT,SIMP] = addsring(VERT,SIMP,LC(2));

    elseif ( ~c(1,2) & ~c(1,3) & c(2,3))

        %%% should only bisect a quadrasection child!
        assertw( (type>=1), 'closure 6' );

        %%% child simplex 1; create/add to rings (reuse parent structure)
        %%% (bisection child noted via the "0" in snew slot 19)
        %%% (generation of child noted in snew slot 20)
        LC=zeros(4,1);
        LC(1)=l;
        snew = [ vn(2),vn_NEW(2,3),vn(1),0, 0,ft(3),ft(1),0, ...
                 LC(1),0,0, 0,0,0,0, 0,0,cmark,ctype,cgen ];
        SIMP(LC(1),:) = snew; 
        [VERT,SIMP] = addsring(VERT,SIMP,LC(1));

        %%% child simplex 2; create/add to rings (new simplex location)
        %%% (bisection child noted via the "0" in snew slot 19)
        %%% (generation of child noted in snew slot 20)
        LL=LL+1;
        LC(2)=LL;
        snew = [ vn(3),vn(1),vn_NEW(2,3),0, 0,ft(1),ft(2),0, ...
                 LC(2),0,0, 0,0,0,0, 0,0,cmark,ctype,cgen ]; 
        SIMP = [ SIMP ; snew ];
        [VERT,SIMP] = addsring(VERT,SIMP,LC(2));

    %%% The three trisection cases
    elseif ( c(1,2) & c(1,3) & ~c(2,3))

        %%% should only bisect a quadrasection child!
        assertw( (type>=1), 'closure 7' );

        %%% child simplex 1; create/add to rings (reuse parent structure)
        %%% (tri-section child noted via the "0" in snew slot 19)
        %%% (generation of child noted in snew slot 20)
        LC=zeros(4,1);
        LC(1)=l;
        snew = [ vn(1),vn_NEW(1,2),vn_NEW(1,3),0, 0,ft(2),ft(3),0, ...
                 LC(1),0,0, 0,0,0,0, 0,0,cmark,ctype,cgen ];
        SIMP(LC(1),:) = snew; 
        [VERT,SIMP] = addsring(VERT,SIMP,LC(1));

        %%% child simplex 2; create/add to rings (new simplex location)
        %%% (tri-section child noted via the "0" in snew slot 19)
        %%% (generation of child noted in snew slot 20)
        LL=LL+1;
        LC(2)=LL;
        snew = [ vn(3),vn_NEW(1,3),vn(2),0, 0,ft(1),ft(2),0, ...
                 LC(2),0,0, 0,0,0,0, 0,0,cmark,ctype,cgen ]; 
        SIMP = [ SIMP ; snew ];
        [VERT,SIMP] = addsring(VERT,SIMP,LC(2));

        %%% child simplex 3; create/add to rings (new simplex location)
        %%% (tri-section child noted via the "0" in snew slot 19)
        %%% (generation of child noted in snew slot 20)
        LL=LL+1;
        LC(3)=LL;
        snew = [ vn(2),vn_NEW(1,3),vn_NEW(1,2),0, 0,ft(3),0,0, ...
                 LC(3),0,0, 0,0,0,0, 0,0,cmark,ctype,cgen ]; 
        SIMP = [ SIMP ; snew ];
        [VERT,SIMP] = addsring(VERT,SIMP,LC(3));

    elseif ( c(1,2) & ~c(1,3) & c(2,3))

        %%% should only bisect a quadrasection child!
        assertw( (type>=1), 'closure 8' );

        %%% child simplex 1; create/add to rings (reuse parent structure)
        %%% (tri-section child noted via the "0" in snew slot 19)
        %%% (generation of child noted in snew slot 20)
        LC=zeros(4,1);
        LC(1)=l;
        snew = [ vn(1),vn_NEW(1,2),vn_NEW(2,3),0, 0,0,ft(3),0, ...
                 LC(1),0,0, 0,0,0,0, 0,0,cmark,ctype,cgen ];
        SIMP(LC(1),:) = snew; 
        [VERT,SIMP] = addsring(VERT,SIMP,LC(1));

        %%% child simplex 2; create/add to rings (new simplex location)
        %%% (tri-section child noted via the "0" in snew slot 19)
        %%% (generation of child noted in snew slot 20)
        LL=LL+1;
        LC(2)=LL;
        snew = [ vn(3),vn(1),vn_NEW(2,3),0, 0,ft(1),ft(2),0, ...
                 LC(2),0,0, 0,0,0,0, 0,0,cmark,ctype,cgen ]; 
        SIMP = [ SIMP ; snew ];
        [VERT,SIMP] = addsring(VERT,SIMP,LC(2));

        %%% child simplex 3; create/add to rings (new simplex location)
        %%% (tri-section child noted via the "0" in snew slot 19)
        %%% (generation of child noted in snew slot 20)
        LL=LL+1;
        LC(3)=LL;
        snew = [ vn(2),vn_NEW(2,3),vn_NEW(1,2),0, 0,ft(3),ft(1),0, ...
                 LC(3),0,0, 0,0,0,0, 0,0,cmark,ctype,cgen ]; 
        SIMP = [ SIMP ; snew ];
        [VERT,SIMP] = addsring(VERT,SIMP,LC(3));

    elseif ( ~c(1,2) & c(1,3) & c(2,3))

        %%% should only bisect a quadrasection child!
        assertw( (type>=1), 'closure 9' );

        %%% child simplex 1; create/add to rings (reuse parent structure)
        %%% (tri-section child noted via the "0" in snew slot 19)
        %%% (generation of child noted in snew slot 20)
        LC=zeros(4,1);
        LC(1)=l;
        snew = [ vn(1),vn_NEW(2,3),vn_NEW(1,3),0, 0,ft(2),0,0, ...
                 LC(1),0,0, 0,0,0,0, 0,0,cmark,ctype,cgen ];
        SIMP(LC(1),:) = snew; 
        [VERT,SIMP] = addsring(VERT,SIMP,LC(1));

        %%% child simplex 2; create/add to rings (new simplex location)
        %%% (tri-section child noted via the "0" in snew slot 19)
        %%% (generation of child noted in snew slot 20)
        LL=LL+1;
        LC(2)=LL;
        snew = [ vn(2),vn_NEW(2,3),vn(1),0, 0,ft(3),ft(1),0, ...
                 LC(2),0,0, 0,0,0,0, 0,0,cmark,ctype,cgen ]; 
        SIMP = [ SIMP ; snew ];
        [VERT,SIMP] = addsring(VERT,SIMP,LC(2));

        %%% child simplex 3; create/add to rings (new simplex location)
        %%% (tri-section child noted via the "0" in snew slot 19)
        %%% (generation of child noted in snew slot 20)
        LL=LL+1;
        LC(3)=LL;
        snew = [ vn(3),vn_NEW(1,3),vn_NEW(2,3),0, 0,ft(1),ft(2),0, ...
                 LC(3),0,0, 0,0,0,0, 0,0,cmark,ctype,cgen ]; 
        SIMP = [ SIMP ; snew ];
        [VERT,SIMP] = addsring(VERT,SIMP,LC(3));

    %%% The one quadrasection case
    elseif ( c(1,2) & c(1,3) & c(2,3))

        %%% should only bisect a quadrasection child!
        assertw( (type>=1), 'closure 10' );

        %%% child simplex 1; create/add to rings (reuse parent structure)
        %%% (non-interior quad child noted via the "1" in snew slot 19)
        %%% (generation of child noted in snew slot 20)
        LC=zeros(4,1);
        LC(1)=l;
        ctype=1;
        snew = [ vn(1),vn_NEW(1,2),vn_NEW(1,3),0, 0,ft(2),ft(3),0, ...
                 LC(1),0,0, 0,0,0,0, 0,0,cmark,ctype,cgen ];
        SIMP(LC(1),:) = snew; 
        [VERT,SIMP] = addsring(VERT,SIMP,LC(1));

        %%% child simplex 2; create/add to rings (new simplex location)
        %%% (non-interior quad child noted via the "1" in snew slot 19)
        %%% (generation of child noted in snew slot 20)
        LL=LL+1;
        LC(2)=LL;
        ctype=1;
        snew = [ vn(2),vn_NEW(2,3),vn_NEW(1,2),0, 0,ft(3),ft(1),0, ...
                 LC(2),0,0, 0,0,0,0, 0,0,cmark,ctype,cgen ]; 
        SIMP = [ SIMP ; snew ];
        [VERT,SIMP] = addsring(VERT,SIMP,LC(2));

        %%% child simplex 3; create/add to rings (new simplex location)
        %%% (non-interior quad child noted via the "1" in snew slot 19)
        %%% (generation of child noted in snew slot 20)
        LL=LL+1;
        LC(3)=LL;
        ctype=1;
        snew = [ vn(3),vn_NEW(1,3),vn_NEW(2,3),0, 0,ft(1),ft(2),0, ...
                 LC(3),0,0, 0,0,0,0, 0,0,cmark,ctype,cgen ];
        SIMP = [ SIMP ; snew ];
        [VERT,SIMP] = addsring(VERT,SIMP,LC(3));
    
        %%% child simplex 4; create/add to rings (new simplex location)
        %%% (interior quad child noted via the "2" in snew slot 19)
        %%% (generation of child noted in snew slot 20)
        LL=LL+1;
        LC(4)=LL;
        ctype=2;
        snew = [ vn_NEW(2,3),vn_NEW(1,3),vn_NEW(1,2),0, 0,0,0,0, ...
                 LC(4),0,0, 0,0,0,0, 0,0,cmark,ctype,cgen ];
        SIMP = [ SIMP ; snew ];
        [VERT,SIMP] = addsring(VERT,SIMP,LC(4));

    %%% The error case
    else
        assertw( 0, 'closure 11' );
    end

end

%%% check new sizes
[N,six]  = size(VERT);
[L,nine] = size(SIMP);
assertw( (L == LL), 'closure 12' );
assertw( (N == NN), 'closure 13' );

%%% return the results
VERT_F = VERT;
SIMP_F = SIMP;
EDGE_F = EDGE;
QUE_F  = NEWQUE;

