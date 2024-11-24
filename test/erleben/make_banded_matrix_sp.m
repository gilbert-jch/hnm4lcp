function A = make_banded_matrix(B,N)
% A = MAKE_BANDED_MATRIX(B,N): Make random banded diagonal symmetric matrix
%
% INPUT:
%
%   B -- Number of off-diagonal bands
%   N -- Number of variables
%
% OUTPUT:
%
%   A -- A banded symmetric positive definite matrix
%
% Copyright 2011, Kenny Erleben, DIKU
% Jean Charles Gilbert, 2024

  % reset random generator

  rng(N)

  % Set A symmetric, B banded, with R(b) on its diagonal b

  A  = speye(N,N);
  R  = -rand(B,1);
  N2 = N*N;
  for b = 1:B
    A(1+b:N+1:N*(N-b)) = R(b);
    A(N*b+1:N+1:N2)    = R(b);
  end

  % modify A in a PSD matrix

  eigmin = min(eig(A));					% correct but expensive
% eigmin = eigs(A,1,'smallestreal');			% correct, but do not converges for n ≥ 175776
% eigmin = min(-eigs(-A));				% incorrect (takes the largest of -A in absolute valeu)
% eigmin = -eigs(-A,1)					% smallest eigenvalue of A (wrong value)
% eigmin = -eigs(-A,1,'largestabs','Tolerance',1.e-6)	% smallest eigenvalue of A (wrong value)
  if eigmin <= 0
    A = A + (0.5-eigmin)*speye(N,N);
  end

%>   [L,d,e] = cholmod_sp(A,1.e-10,1.e+10);
%>   A = A+diag(e);
