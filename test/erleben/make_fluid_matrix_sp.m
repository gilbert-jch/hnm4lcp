function [A] = make_fluid_matrix_sp(G)
%%
% [A] = make_fluid_matrix_sp(G)
% 
% Make matrix corresponding to a Poisson equation. Sparse version.
%
% INPUT:
%
%   G -- Grid cell size.
%
% OUTPUT:
%
%   A -- The LCP coefficient matrix.
%
% Copyright 2011, Kenny Erleben, DIKU
% Jean Charles Gilbert, 2024

%--- Generate a random fluid matrix for 3D fluid problem on a -------------
%--- grid of size GxGxG ---------------------------------------------------

  B = 3;      % Number of off-diagonal bands
  N = G*G*G;  % Number of variables
  A = make_banded_matrix_sp(B,N);
