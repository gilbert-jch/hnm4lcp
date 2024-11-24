function [M,q,x,x1] = csizmadia(n,var,verb)

%%
% [M,q,x,x1] = csizmadia(n,var,verb)
%
% Return data for linear complementarity problem of Csizmadia with n
% variables. The matrix M has 1 on its diagonal, 0 on its upper
% triangular part and -1 on its lower triangular part. It is therefore
% an M-matrix. The vector q is e-M*e, where e is the vector of all ones.
% The problem is known to be difficult to solve for some interior point
% algorithms because the matrix has an exponential handicap.
%
% If verb = true, there are some printings.
%
% References:
% - E. de Klerk and M. E.-Nagy (2011), "On the complexity of computing
%   the handicap of a sufficient matrix", Mathematical Programming,
%   129:2, 383-402.
% - Z. Darvay, T. Illés, J. Povh, P.R. Rigó (2020), "Predictor-corrector
%   interior-point algorithm for sufficient linear complementarity
%   problems based on a new search direction. SIAM Journal on
%   Optimization, 30:3, 2628–2658.

  M  = [];
  q  = [];
  x  = [];
  x1 = [];

% Introductory message

  cl = clock;
  dstr = sprintf('%i-%i-%i, %i:%i:%i',cl(3),cl(2),cl(1),cl(4),cl(5),fix(cl(6)));
  if verb
    fprintf('\nCsizmadia problem');
    fprintf('\n. variant = %i',var);
    fprintf('\n. n       = %i\n',n);
  end

% Compute M

  M = 2*eye(n)-tril(ones(n));

% Set q, x and x1

  e  = ones(n,1);

  if var == 0	% values set by Z. Darvay, T. Illés, J. Povh, P.R. Rigó (2020)

    q  = e - M*e;		% q = (0, 1, 2, ..., n-1)
    x  = zeros(n,1);
    x1 = e;

  elseif var == 1		% personal values to avoid q >= 0, x1 = e

    x  = rem([1:n]',2);		% x     = (1, 0, 1, 0, 1, ...)
    q  = e - x - M*x;		% M*x+q = (0, 1, 0, 1, 0, ...)
    x1 = e;

  elseif var == 2		% personal values to avoid q >= 0, x1 = 0

    x  = rem([1:n]',2);		% x     = (1, 0, 1, 0, 1, ...)
    q  = e - x - M*x;		% M*x+q = (0, 1, 0, 1, 0, ...)
    x1 = zeros(n,1);

  end

end
