function [SIMP_F,QUE_F] = mark(VERT,SIMP,U,Ug,markkey,level,estkey,tol);
%MARK  Mark some simplices for refinement
%
% Usage: [SIMP_F,QUE_F] = mark(VERT,SIMP,U,Ug,markkey,level,estkey,tol);
%
% Input:
%
%    VERT    = <see "read.m" for description of the datastructure>
%    SIMP    = <see "read.m" for description of the datastructure>
%    markkey = which marking restrictions are enabled
%    level   = mesh level for refinement (if marking is restricted)
%    estkey  = which error estimator to use to use
%    tol     = indicator tolerance (mark simplex if indicator exceeds tol)
%
% Output:
%
%    SIMP_F = marked set of simplices
%    QUE_F  = list of the simplex numbers that were marked for refinement
%
% Notes:
%
%    markkey ==> 0=no restriction on marking for refinement
%                1=mark ONLY elements where were created as the
%                  INTERIOR child in a quadrasection in the
%                  previous level.
%                2=mark ONLY elements which were created as
%                  ANY child in a quadrasection in the
%                  next previous level.
%
%    estkey  ==> 0=mark all
%                1=geometric marking
%                2=a posteriori indicator via nonlinear strong residuals
%

%
% Author:   Michael Holst
% rcsid="$Id: mark.m,v 1.1.1.1 2007/04/27 08:28:06 hrg Exp $"

%%% get some sizes
[N,eight]     = size(VERT);
[L,seventeen] = size(SIMP);

%%% build the simplex rings (assume done already!)
%%% [VERT,SIMP] = bldsring(VERT,SIMP);

%%% we always place the elements into the first q
whichq = 1;

%%% initialize the queue of marked simplices
QUE = [];

%%% the uniform refinement indicator
if (estkey==0)
    fprintf('..uniform marking..');
    for l=1:L
        mark = SIMP(l,18);
        % mark = 1;
        type = SIMP(l,19);
        gen  = SIMP(l,20);
        if (level==1)
            SIMP(l,whichq+15) = 1;
        elseif (markkey==0)
            SIMP(l,whichq+15) = 1;
        elseif (markkey==1)
            if (mark & (gen==level) & (type==2))
                SIMP(l,whichq+15) = 1;
            end
        elseif (markkey==2)
            if (mark & (gen==level) & (type>=1))
                SIMP(l,whichq+15) = 1;
            end
        end
    end
end;

%%% the geometric test indicator
if (estkey==1)
    fprintf('..geometric marking..');
    rad = .25;
    for l=1:L
        d1 = norm( VERT( SIMP(l,1), 1:3 ) );
        d2 = norm( VERT( SIMP(l,2), 1:3 ) );
        d3 = norm( VERT( SIMP(l,3), 1:3 ) );
        if ( ~ ((d1 < rad) & (d2 < rad) & (d3 < rad)) ...
           & ~ ((d1 > rad) & (d2 > rad) & (d3 > rad)) )
            markit = 1;
        else
            markit = 0;
        end
        if (markit)
            mark = SIMP(l,18);
            % mark = 1;
            type = SIMP(l,19);
            gen  = SIMP(l,20);
            if (level==1)
                SIMP(l,whichq+15) = 1;
            elseif (markkey==0)
                SIMP(l,whichq+15) = 1;
            elseif (markkey==1)
                if (mark & (gen==level) & (type==2))
                    SIMP(l,whichq+15) = 1;
                end
            elseif (markkey==2)
                if (mark & (gen==level) & (type>=1))
                    SIMP(l,whichq+15) = 1;
                end
            end
        end
    end
end

%%% the strong residual indicator
if (estkey==2)
    fprintf('..error marking..');
    indicator = snresid(VERT,SIMP,U,Ug);
    failed = find( indicator > tol );
    [len,one] = size(failed);
    for ll=1:len
        l=failed(ll);
        mark = SIMP(l,18);
        % mark = 1;
        type = SIMP(l,19);
        gen  = SIMP(l,20);
        if (level==1)
            SIMP(l,whichq+15) = 1;
        elseif (markkey==0)
            SIMP(l,whichq+15) = 1;
        elseif (markkey==1)
            if (mark & (gen==level) & (type==2))
                SIMP(l,whichq+15) = 1;
            end
        elseif (markkey==2)
            if (mark & (gen==level) & (type>=1))
                SIMP(l,whichq+15) = 1;
            end
        end
    end
end

%%% build a queue of marked simplices
QUE = find( SIMP(:,whichq+15) == 1 );
[P,one] = size(QUE);

%%% set new legal marking subset
SIMP(:,18) = zeros(L,1);
SIMP(QUE,18) = ones(P,1);

%%% kill the simplex rings (assume done already!)
%%% [VERT,SIMP] = kilring(VERT,SIMP);

%%% return the marked simplices
SIMP_F = SIMP;
QUE_F  = QUE;

