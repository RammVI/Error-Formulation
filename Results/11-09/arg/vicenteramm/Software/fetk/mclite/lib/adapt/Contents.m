% Error estimation and simplex subdivision library
%
% This library implements high-level simplex mesh adaptivity algorithms
% and (a posteriori) error estimation algorithms.
%
% The available routines are as follows:
%
%      ADDONERI      add vertex 1-rings of simplices to refinement queue
%      CLOSURE       one-shot closure of a red-refined mesh via quad/tri/bisect
%      MARK          tag simplices for refinement (e.g., error estim.)
%      REFIN         refine a set of marked simplices until conformity
%      REFINBIS      one pass through simplices, bisecting marked ones
%      REFINQUD      one pass through simplices, quadrasect marked ones
%      SNRESID       strong nonlinear residual error indicator
%
% Author:   Michael Holst
% rcsid="$Id: Contents.m,v 1.1.1.1 2007/04/27 08:28:05 hrg Exp $"
