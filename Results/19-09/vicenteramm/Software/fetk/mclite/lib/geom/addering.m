function [VERT_F,EDGE_F] = addering(VERT,EDGE,which);
%ADDERING  Add an edge to the edge rings associated with each vertex
%
% Usage: [VERT_F,EDGE_F] = addering(VERT,EDGE,which);
%
% Input:
%
%    which  = the edge number to add to rings
%    VERT   = <see "read.m" for description of the datastructure>
%    EDGE   = <see "read.m" for description of the datastructure>
%    
% Output:
%
%    VERT_F = list of vertices in the mesh with the new ring info
%    EDGE_F = list of edges in the mesh with the new ring info
%
% Author:   Michael Holst
% rcsid="$Id: addering.m,v 1.1.1.1 2007/04/27 08:28:06 hrg Exp $"

%%% add the "which" edge to the two rings

   for i=1:2
      oldFirstE = VERT( EDGE(which,i) , 8 );
      VERT( EDGE(which,i) , 8 ) = which;
      EDGE(which,i+6) = oldFirstE;
   end

%%% return the results

   VERT_F = VERT;
   EDGE_F = EDGE;

