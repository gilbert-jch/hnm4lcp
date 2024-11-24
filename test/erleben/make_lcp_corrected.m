function [x, q] = make_lcp_corrected(M,F)

% [x, q] = MAKE_LCP(M,F): Make LCP
%
% INPUT:
%
%   M -- The coefficient matrix in the LCP
%   F -- The fraction of zero-values in the x-solution.
%
% OUTPUT:
%
%   x -- A solution for the LCP problem.
%   q -- The right hand side vector
%
% Copyright 2011, Kenny Erleben, DIKU

%--- Get number of variables ----------------------------------------------

n = size(M,1);

%--- Generate a random LCP solution ---------------------------------------

x    = rand(n,1);
I    = find(x<F);
x(I) = 0;

%--- Generate a right-hand-side vector that is ----------------------------
%--- consistent with random solution           ----------------------------

y    = zeros(n,1);
y(I) = rand(length(I),1);	% hence y(I) > 0, while x(I) = 0;
q    = y-M*x;			% hence M*x+q = y ≥ 0
end
