function [VERT_F,SIMP_F,EDGE_F,QUE_F] = refinbis(VERT,SIMP,EDGE,QUE, ...
                                                 whichq,closekey);
%REFINBIS  Refine the finite element mesh using BISECTION
%
% Usage: [VERT_F,SIMP_F,EDGE_F,QUE_F] = refinbis(VERT,SIMP,EDGE,QUE, ...
%                                                whichq,closekey);
%
% Input:
%
%    VERT     = <see "read.m" for description of the datastructure>
%    SIMP     = <see "read.m" for description of the datastructure>
%    EDGE     = <see "read.m" for description of the datastructure>
%    QUE      = list of simplices to be refined
%    whichq   = which que we are currently using
%    closekey = do we worry about closure for the children
% 
% Output:
%    
%    VERT_F = list of vertices in the refined mesh
%    SIMP_F = conformally refined set of simplices
%    EDGE_F = list of edges created during refinement
%    QUE_F  = new list of simplices to be refined
%
% A useful mapping:
%
%    mapV[3][3] ==> The vertex permutation map; this mapping allows
%                   us to map any of the three refinement
%                   scenarios onto the canonical one.
%                   Each row represents the appropriate vertex
%                   number mappings for each of the three cases.
%
% Author:   Michael Holst
% rcsid="$Id: refinbis.m,v 1.1.1.1 2007/04/27 08:28:06 hrg Exp $"

%%% this mapping allows us to map the three refinement cases to one case
mapV  = [ 1,  2,  3;
          3,  1,  2;
          2,  3,  1 ];

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
assertw( ((whichq==1)&(swapq==2)) | ((whichq==2)&(swapq==1)), 'refinbis 5' );

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
    assertw( (refq ~= 0), 'refinbis 1' );

    %%% clear our que tags
    SIMP(l,whichq+15) = 0;
    SIMP(l,swapq+15)  = 0;

    %%% get our three nabor simplex numbers before killing the rings
    nabors = zeros(3,1);
    [nabors(1),cf1] = getnabor(VERT,SIMP,l,1);  %%% nabor opposite vertex 1
    [nabors(2),cf2] = getnabor(VERT,SIMP,l,2);  %%% nabor opposite vertex 2
    [nabors(3),cf3] = getnabor(VERT,SIMP,l,3);  %%% nabor opposite vertex 3

    %%% remove this simplex from rings (will reuse as first child)
    [VERT,SIMP] = delsring(VERT,SIMP,l);

    %%% which is longest edge
    canon = longedge(VERT,SIMP,l);

    %%% face/vertex info for the T vertices/faces forming this simplex
    vl = zeros(T);      %%% the canonical mapping array
    vn = zeros(T);      %%% vertex numbers of this simp
    ve = zeros(T);      %%% first edge of each vertices' edge list
    v  = zeros(T,3);    %%% the x,y,z coordinates of the vertices
    ft = zeros(3,1);    %%% the 3 integer face types
    nb = zeros(3,1);    %%% the 3 nabors
    for i=1:T
        vl(i)  = mapV(canon,i);
        vn(i)  = SIMP( l, vl(i) );
        ft(i)  = SIMP( l, vl(i)+4 );
        ve(i)  = VERT( vn(i), 8 );
        v(i,:) = VERT( vn(i), 1:3 );
        nb(i)  = nabors( vl(i) );
    end

    %%% get existing vertices or make new ones as needed
    edg = ve(1);
    found = 0;
    while ((found == 0) & (edg ~= 0))
        ev = EDGE(edg,1:3);
        if ( ((ev(1) == vn(1)) & (ev(2) == vn(2)) ) ...
            | ((ev(1) == vn(2)) & (ev(2) == vn(1)) ) )
            found = 1;
        else
            ii=-1;
            if ( ev(1) == vn(1) ) ii = 1; end;
            if ( ev(2) == vn(1) ) ii = 2; end;
            assertw( (ii >= 0), 'refinbis 2' );
            edg = EDGE( edg, ii+6);
        end
    end

    if (found == 0)
        %%% determine type for new vertex
        vt_NEW = ft(3);

        %%% create the new vertex
        NN=NN+1;
        vn_NEW = NN;
        v_NEW = edgsplit( v(1,:), v(2,:) );
        vnew = [ v_NEW, NN,0,vt_NEW, 0,0 ];
        VERT = [ VERT ; vnew ];

        %%% create the marking edge and add to rings
        KK=KK+1;
        enew = [ vn(1),vn(2),vn_NEW, KK,0,vt_NEW, 0,0 ];
        EDGE = [ EDGE ; enew ];
        [VERT,EDGE] = addering(VERT,EDGE,KK);

        %%% we created edge, so have to notify nonconforming nabor
        %%% if nabor on either queue already, we don't do anything
        tnab = nb(3);
        if (tnab > 0)
            if ( ~((SIMP(tnab,whichq+15)>0) | (SIMP(tnab,swapq+15)>0)) )
                SIMP(tnab,swapq+15) = 1;
                NEWQUE = [ NEWQUE ; tnab ];
            end
        end
    else
        vn_NEW = ev(3);
    end

    %%% child simplex 1; create/add to rings (reuse parent structure)
    %%% (bisection child noted via the "0" in snew slot 18)
    %%% (generation of child noted in snew slot 20)
    LC=zeros(2,1);
    LC(1)=l;
    snew = [ vn(1),vn_NEW,vn(3),0, 0,ft(2),ft(3),0, ...
             LC(1),0,0, 0,0,0,0, 0,0,cmark,ctype,cgen ];
    SIMP(LC(1),:) = snew; 
    [VERT,SIMP] = addsring(VERT,SIMP,LC(1));

    %%% child simplex 2; create/add to rings (new simplex location)
    %%% (bisection child noted via the "0" in snew slot 18)
    %%% (generation of child noted in snew slot 20)
    LL=LL+1;
    LC(2)=LL;
    snew = [ vn(2),vn(3),vn_NEW,0, 0,ft(3),ft(1),0, ...
             LC(2),0,0, 0,0,0,0, 0,0,cmark,ctype,cgen ];
    SIMP = [ SIMP ; snew ];
    [VERT,SIMP] = addsring(VERT,SIMP,LC(2));

    %%% see if children are nonconforming (only do for recursive closure)
    if (closekey == 2)
        edgexx = zeros(3,1);
        for child=1:2
            vv        = SIMP(LC(child),1:3);
            edgexx(1) = getedge(VERT,EDGE,vv(2),vv(3));  %%% edge opposite v1
            edgexx(2) = getedge(VERT,EDGE,vv(1),vv(3));  %%% edge opposite v2
            edgexx(3) = getedge(VERT,EDGE,vv(1),vv(2));  %%% edge opposite v3
            if (norm(edgexx) > 0)
                 SIMP(LC(child),swapq+15) = 1;
                 NEWQUE = [ NEWQUE ; LC(child) ];
            end
        end
    end
end

%%% check new sizes
[N,six]  = size(VERT);
[L,nine] = size(SIMP);
assertw( (L == LL), 'refinbis 3' );
assertw( (N == NN), 'refinbis 4' );

%%% return the results
VERT_F = VERT;
SIMP_F = SIMP;
EDGE_F = EDGE;
QUE_F  = NEWQUE;

