function [M,q,x,x1] = bg2012(n,verb)

%%
% [M,q,x,x1] = bg2012(n,verb)
%
% Return data for the BenGharbia-Gilbert problem with n variables (n
% muts be ≥ 3). On return, M is a sparse P-matrix. Data are taken from
%
%   I. Ben Gharbia, J.Ch. Gilbert (2012), "Nonconvergence of the plain
%   Newton-min algorithm for linear complementarity problems with a
%   P-matrix, Mathematical Programming, 134, 349-364,
%   http://dx.doi.org/10.1007/s10107-010-0439-6.
%
% If verb is false, no printing.

  M  = [];	% matrix M of the problem
  q  = [];	% vector q of the problem
  x  = [];	% solution to the problem
  x1 = [];	% suggested initial point

  if isempty(n) | (n < 3)
    if verb; fprintf('\n\n### bg2012: the dimension n = %i must be ≥ 3\n\n',n); end
    return
  end

% Introductory message

  cl = clock;
  dstr = sprintf('%i-%i-%i, %i:%i:%i',cl(3),cl(2),cl(1),cl(4),cl(5),fix(cl(6)));
  if verb
    fprintf('\nFathi problem (%s)',dstr);
    fprintf('\n. n = %i\n',n);
  end

% Compute M

  M = speye(n);

  if fix(n/2)*2 ~= n	% n is odd
    alpha  = 2;
    a      = alpha*ones(n,1);
    M      = spdiags(a,-1,M);
    M(1,n) = alpha;
  else			% n is even
    alpha  = 4/3;
    a      = alpha*ones(n,1);
    M      = spdiags(a,-1,M);
    M(1,n) = alpha;
    beta   = 0.5;
    b      = beta*ones(n,1);
    M      = spdiags(b,-2,M);
    M      = spdiags(b,n-2,M);
  end

% Set q

  q = ones(n,1);

% Set x (the solution)

  x = zeros(n,1);

% Set x1 (initial point)

  x1    = zeros(n,1);
  x1(1) = -1;

end
