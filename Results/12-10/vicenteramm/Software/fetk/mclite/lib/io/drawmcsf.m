function h = drawmcsf(vert,simp,varargin)
%DRAWMCSF  Draw a 2D mesh from MCSF (MC-Simplex-format) file
%
%   If VERT and SIMP have been defined by executing a 2D MCSF file,
%   DRAWMCSF(VERT,SIMP) draws the mesh.
%
%   DRAWMCSF(VERT,SIMP,C) uses C to specify the coloring.  The
%   default is to use the Z coordinate, i.e., color by height.
%   DRAWMCSF(VERT,SIMP,0) will use a constant color.
%
%   DRAWMCSF(...,'param','value','param','value',...) allows
%   patch param/value pairs to be used when creating the patch
%   object.  E.g., DRAWMCSF(VERT,SIMP,'FACECOLOR','INTERP').
%
%   A handle may be specified as an output argument for later
%   handle graphics manipulation: H = DRAWMCSF(...)
%
%   It may be useful to use "axis equal" and/or "view(2)" (for flat
%   meshes) and/or "rotate3d" (for surfaces) with the displayed mesh.
%
%   See also PATCH, TRISURF.
%
% Author:   Doug Arnold
% rcsid="$Id: drawmcsf.m,v 1.1.1.1 2007/04/27 08:28:07 hrg Exp $"

if nargin < 2
  error('Requires at least 2 input arguments')
end
ss = size(simp); sv = size(vert);
if (ss(2) ~= 9)
  error('simp must by M x 9')
end
if (sv(2) ~= 5)
  error('vert must be N x 5')
end

if nargout > 0
  h = trisurf(simp(:,7:9)+1,vert(:,3),vert(:,4),vert(:,5),varargin{:});
else
  trisurf(simp(:,7:9)+1,vert(:,3),vert(:,4),vert(:,5),varargin{:});
end

