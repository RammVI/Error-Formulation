function [VERT_F,SIMP_F,EDGE_F,QUE_F] = addoneri(VERT,SIMP,EDGE,QUE,whichq);
%ADDONERI  Add vertex 1-rings of simplices to the refinement queue
%
% Usage: [VERT_F,SIMP_F,EDGE_F,QUE_F] = addoneri(VERT,SIMP,EDGE,QUE,whichq);
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
% rcsid="$Id: addoneri.m,v 1.1.1.1 2007/04/27 08:28:05 hrg Exp $"

%%% recover array dimensions
[N,eight]     = size(VERT);
[L,seventeen] = size(SIMP); 
[K,eight]     = size(EDGE);
[P,one]       = size(QUE);

%%% loop over the elements
NEWQUE = QUE;
for refGuy = 1:P

    %%% get my name
    l=QUE(refGuy);
    gen = SIMP(l,20);

    %%% check que tag status
    refq = SIMP(l,whichq+15);
    assertw( (refq ~= 0), 'addoneri 1' );

    %%% the vertices making up this simplex
    v1 = SIMP(l,1);
    v2 = SIMP(l,2);
    v3 = SIMP(l,3);

    %%% vertex 1: add simplices in the 1-rings to refinement queue
    vv = v1;
    ll = VERT( vv, 7 );
    while (ll ~= 0)
        %%% first do a level check!
        gg = SIMP(ll,20);
        assertw( (gg == gen), 'addoneri 2');
        if (gg ~= gen)
           fprintf('gg=<%d>  gen=<%d>\n',gg,gen);
        end
        %%% add guy to the queue
        if (SIMP(ll,whichq+15) == 0)
            SIMP(ll,whichq+15) = 1;
            NEWQUE = [ NEWQUE ; ll ];
        end
        %%% go to the next simplex in the simplex ring
        locv = SIMP(ll,1:3);
        lfind = 0;
        if (locv(1)==vv) lfind = 1; end;
        if (locv(2)==vv) lfind = 2; end;
        if (locv(3)==vv) lfind = 3; end;
        assertw( (lfind>0), 'addoneri 3' );
        ll = SIMP(ll,lfind+11);
    end

    %%% vertex 2: add simplices in the 1-rings to refinement queue
    vv = v2;
    ll = VERT( vv, 7 );
    while (ll ~= 0)
        %%% first do a level check!
        gg = SIMP(ll,20);
        assertw( (gg == gen), 'addoneri 4');
        if (gg ~= gen)
           fprintf('gg=<%d>  gen=<%d>\n',gg,gen);
        end
        %%% add guy to the queue
        if (SIMP(ll,whichq+15) == 0)
            SIMP(ll,whichq+15) = 1;
            NEWQUE = [ NEWQUE ; ll ];
        end
        %%% go to the next simplex in the simplex ring
        locv = SIMP(ll,1:3);
        lfind = 0;
        if (locv(1)==vv) lfind = 1; end;
        if (locv(2)==vv) lfind = 2; end;
        if (locv(3)==vv) lfind = 3; end;
        assertw( (lfind>0), 'addoneri 5' );
        ll = SIMP(ll,lfind+11);
    end


    %%% vertex 3: add simplices in the 1-rings to refinement queue
    vv = v3;
    ll = VERT( vv, 7 );
    while (ll ~= 0)
        %%% first do a level check!
        gg = SIMP(ll,20);
        assertw( (gg == gen), 'addoneri 6');
        if (gg ~= gen)
           fprintf('gg=<%d>  gen=<%d>\n',gg,gen);
        end
        %%% add guy to the queue
        if (SIMP(ll,whichq+15) == 0)
            SIMP(ll,whichq+15) = 1;
            NEWQUE = [ NEWQUE ; ll ];
        end
        %%% go to the next simplex in the simplex ring
        locv = SIMP(ll,1:3);
        lfind = 0;
        if (locv(1)==vv) lfind = 1; end;
        if (locv(2)==vv) lfind = 2; end;
        if (locv(3)==vv) lfind = 3; end;
        assertw( (lfind>0), 'addoneri 7' );
        ll = SIMP(ll,lfind+11);
    end

end

%%% return the results
VERT_F = VERT;
SIMP_F = SIMP;
EDGE_F = EDGE;
QUE_F  = NEWQUE;

