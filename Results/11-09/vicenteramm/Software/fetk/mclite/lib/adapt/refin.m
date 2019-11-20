function [VERT_F,SIMP_F,P_CF,ONER] = refin(VERT,SIMP,QUE,refkey,closekey);
%REFIN  Refine the finite element mesh
%
% Usage: [VERT_F,SIMP_F,P_CF,ONER] = refin(VERT,SIMP,QUE,refkey,closekey);
%
%    Refine the finite element mesh and produce a linear 
%    prolongation operator which relates functions on the two meshes.
%    We perform a conformity loop after the initial refinement loop,
%    and continue until a conforming mesh is obtained.
%
% Input:
%
%    refkey   = which refinement stategy to use
%    closekey = which closeure stategy to use
%    VERT     = <see "read.m" for description of the datastructure>
%    SIMP     = <see "read.m" for description of the datastructure>
%    QUE      = list of marked simplices to be refined
% 
% Output:
%    
%    VERT_F = list of vertices in the refined mesh
%    SIMP_F = conformally refined set of simplices
%    P_CF   = prolongation matrix relating previous and current meshes
%    ONER   = list of new vertices and their one-ring edge nabor vertices
%
% Author:   Michael Holst
% rcsid="$Id: refin.m,v 1.1.1.1 2007/04/27 08:28:06 hrg Exp $"

%%% sizes
EDGE = [];
[numRc,whatever] = size(VERT);

%%% build simplex rings (assume done already!)
%%% [VERT,SIMP] = bldsring(VERT,SIMP);

%%% while we have some guys to refine, do it
[numQ,one] = size(QUE);
whichq = 1;
if (numQ > 0)

    %%% initial refinement of marked simplices
    if (refkey == 0)
        %%% bisection (bisect marked guys by longest edge)
        [VERT,SIMP,EDGE,QUE] = refinbis(VERT,SIMP,EDGE,QUE,whichq,closekey);
    elseif (refkey == 1)
        %%% quadrasection-pure (quadrasect marked guys)
        [VERT,SIMP,EDGE,QUE] = refinqud(VERT,SIMP,EDGE,QUE,whichq,closekey);
    elseif (refkey == 2)
        %%% quadrasection-patch (quadrasect marked + 1-rings)
        [VERT,SIMP,EDGE,QUE] = addoneri(VERT,SIMP,EDGE,QUE,whichq);
        [VERT,SIMP,EDGE,QUE] = refinqud(VERT,SIMP,EDGE,QUE,whichq,closekey);
    else
        assertw( 0, 'refin 1' );
    end

    %%% queue swap
    [numQ,one] = size(QUE);
    if (whichq == 1)
         whichq = 2;
    else
         whichq = 1;
    end

    %%% closure
    if (closekey == 0)
        %%% no-op (i.e., return a pure red-refined mesh)
    elseif (closekey == 1)
        %%% green (pure green closure; quad-/tri-/bi-section)
        [VERT,SIMP,EDGE,QUE] = closure(VERT,SIMP,EDGE,QUE,whichq);
    elseif (closekey == 2)
        %%% bisection (RECURSIVE bisection by longest edge)
        iter=0;
        while (numQ > 0)
            iter=iter+1;
            [VERT,SIMP,EDGE,QUE] = refinbis(VERT,SIMP,EDGE,QUE,whichq,closekey);
            %%% queue swap
            [numQ,one] = size(QUE);
            if (whichq == 1)
                 whichq = 2;
            else
                 whichq = 1;
            end
        end
    else
        assertw( 0, 'refin 2' );
    end
end

%%% build the prolongation operator (and initialize 1-ring list)
[numR,whatever] = size(VERT);
[numE,whatever] = size(EDGE);
assertw( (numE == (numR-numRc)), 'refin 3' );
PRO = speye(numRc,numRc);
for z = 1:numRc
    if (vdiri(VERT(z,6)))
        PRO(z,:) =  zeros(1,numRc);
    end;
end;
RING = [];
for k=1:numE
    v0    = EDGE(k,1);
    v1    = EDGE(k,2);
    v01   = EDGE(k,3);
    v01tp = EDGE(k,6);
    assertw( (v01 > numRc), 'refin 4' );
    thisRow = zeros(1,numRc);
    if ( ~vdiri(v01tp) )
        RING  = [ RING ; v01 ];
        if ( (v0 <= numRc) & (v1 <= numRc) )
            thisRow(v0) = 0.5;
            thisRow(v1) = 0.5;
        else
            assertw( 0, 'refin 5' );
        end
    end
    PRO = [ PRO ; thisRow ];
end

%%% build the 1-ring vertex list around newly introduced vertices
[RING_len,one] = size(RING);
for k=1:RING_len

    %%% first grab vertex id of one of the new vertices
    me = RING(k);

    %%% now go through simplex ring and add all vertices to 1-ring list
    ll = VERT( me, 7 );
    while (ll ~= 0)

        %%% grab three vertices and vertex types of 1-ring neighbors
        locv = SIMP(ll,1:3);
        typv = VERT(locv,6);

        %%% add all non-dirichlet guys to list (we deal with duplicates later)
        locvv = [];
        for i=1:3
            if ( ~vdiri(typv(i)) )
                locvv = [ locvv , locv(i) ];
            end;
        end;
        RING = [ RING ; locvv' ];

        %%% go to the next simplex in the simplex ring
        lfind = 0;
        if (locv(1)==me) lfind = 1; end;
        if (locv(2)==me) lfind = 2; end;
        if (locv(3)==me) lfind = 3; end;
        assertw( (lfind>0), 'addoneri 3' );
        ll = SIMP(ll,lfind+11);
    end
end;    

%%% condense the 1-ring vertex list to just the unique vertices
RING3 = [];
[LL,one] = size(RING);
if (LL > 0)
    RING2 = sort(RING);
    RING3 = [ RING3 ; RING2(1) ];
    jj = 1;
    for ii=2:LL
        if (RING2(ii) ~= RING3(jj))
            jj = jj + 1;
            RING3 = [ RING3 ; RING2(ii) ];
        end
    end
end

%%% destroy edges before returning
EDGE = [];

%%% kill simplex rings (assume done already!)
%%% [VERT,SIMP] = kilring(VERT,SIMP);

%%% return the results
VERT_F = VERT;
SIMP_F = SIMP;
P_CF   = PRO;
ONER   = RING3;

