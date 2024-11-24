function [M,q,x,x1] = murty(n,verb)

%%
% [M,q,x,x1] = murty(n,verb)
%
% Return data for the linear complementarity problem of Murty, which
% originated in
%
%   Katta G. Murty (1978). "Computational complexity of complementarity
%   pivot methods". Mathematical Programming Study, 7, 61-73.
%
% If verb ~= 0, there are some printings.

% Introductory message

  cl = clock;
  dstr = sprintf('%i-%i-%i, %i:%i:%i',cl(3),cl(2),cl(1),cl(4),cl(5),fix(cl(6)));
  if verb
    fprintf('\nMurty problem (%s)',dstr);
    fprintf('\n. n = %i\n',n);
  end

% Compute M

  M = 2*ones(n,n);
  M = tril(M);
  M = M - diag(ones(n,1));

% Set x (the solution)

  x = zeros(n,1);
  x(1) = 1;

% Set q

  q = -ones(n,1);

% Set intial x

  x1 = zeros(n,1);

end
