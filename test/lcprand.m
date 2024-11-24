function [M,q,x,x1] = lcprand(n,na,ni,verb)

%%
% [M,q,x,x1] = lcprand(n,na,ni,verb)
%
% Return data for a random linear complementarity problem 0 <= x _|_
% (M*x+q) >= 0, with a solution x and y=M*x+q satisfying
%
%    n  = length(x),
%    na = length(find(x < y)),
%    ni = length(find(x > y)).
%
% Then, there are ne=n-na-ni indices i such that x(i)=y(i)=0. The vector
% x1 on output is a possible initial point for an iterative method.
%
% If verb is zero or false, the code works without printing anything.

  M  = [];
  q  = [];
  x  = [];
  x1 = [];

% Introductory message

  if verb
    cl = clock;
    dstr = sprintf('%i-%i-%i, %i:%i:%i',cl(3),cl(2),cl(1),cl(4),cl(5),fix(cl(6)));
    fprintf('\nRandom LCP problem (%s)\n',dstr);
  end

  ne = n-na-ni;
  if verb; fprintf('\nn = %i, na = %i, ne = %i, ni = %i\n',n,na,ne,ni); end
  if (n <= 0) | (na < 0) | (ne < 0) | (ni < 0)
    if verb; fprintf('Inconsistent dimensions\n\n'); end
    return
  end

% Compute M

  rng('default');	% initialize the random generator to its default setting

  A = 10*rand(n,n)-5;
  B = 10*rand(n,n)-5;
  B = 0.5*(B-B');	% B is skew-symmetric
  d = 0.3*rand(n,1);
  M = A'*A+B+diag(d);

% Set x (the solution)

  x            = zeros(n,1);
  x(na+ne+1:n) = rand(ni,1);

% Compute q

  p       = M*x;
  pmax    = norm(p(1:na),Inf);
  y       = zeros(n,1);
  y(1:na) = pmax*rand(na,1);
  q       = y-p;	% Hence y = M*x+q

% Set initial x

  x1 = 100*rand(n,1)-50;

end
