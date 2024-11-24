function [M,q,x,x1] = fathi(n,p,verb)

%%
% [M,q,x,x1] = fathi(n,p,verb)
%
% Return data for the p-th linear complementarity problem of Fathi with
% n variables. Data taken from
%
% Y. Fathy (1979). "Computational complexity of LCPs associated with
% positive definite symmetric matrices". Mathematical Programming, 17,
% 335–344.
%
% If verb ~= 0, there are some printings.

  M  = [];
  q  = [];
  x  = [];
  x1 = [];

  if isempty(p) | (p < 1) | (p > 3)
    if verb; fprintf('\n\n### fathi: unknown problem %i\n\n',p); end
    return
  end

% Introductory message

  cl = clock;
  dstr = sprintf('%i-%i-%i, %i:%i:%i',cl(3),cl(2),cl(1),cl(4),cl(5),fix(cl(6)));
  if verb
    fprintf('\nFathi problem %i (%s)',p,dstr);
    fprintf('\n. n = %i\n',n);
  end

% Compute M

  M = 2*ones(n,n);
  M = tril(M);
  M = M - diag(ones(n,1));
  M = M*M';

  switch p
  case 1
  case 2
    E = ones(n,n);
    M = M+4*E;
  case 3
    E = ones(n,n);
    E(:,1) = 0;
    M = M-4*E;
  end

% Set x (the solution)

  x    = zeros(n,1);
  x(1) = 1;

% Set q

  q = -ones(n,1);

% Set initial x

  x1 = zeros(n,1);

end
