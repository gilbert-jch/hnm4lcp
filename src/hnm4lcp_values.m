function [values] = hnm4lcp_values ()

%%
% [values] = hnm4lcp_values ()
%
% Set the constant values used in NMAR.

% Version 1.0, November 2024.
%
% HAL Id: hal-04799965
%
% If you found this piece of software useful, please cite the paper:
% Jean-Pierre Dussault, Mathieu Frappier, Jean Charles Gilbert (2024).
% 'Polyhedral {Newton}-min algorithms for complementarity problems',
% Mathematical Programming (in revision). See also
% http://hal.inria.fr/hal-02306526/en.

%-----------------------------------------------------------------------
%
% Authors:
% - Jean-Pierre Dussault (Univ. of Sherbrooke, Canada),
% - Jean Charles Gilbert (INRIA, France) & (Univ. of Sherbrooke, Canada)
%
% Copyright 2024, INRIA (France) and Univ. of Sherbrooke (Canada).
%
% Hnm4lcp is distributed under the terms of the Q Public License version
% 1.0.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the Q Public
% License version 1.0 for more details.
%
% You should have received a copy of the Q Public License version 1.0
% along with this program.  If not, see
% <https://doc.qt.io/archives/3.3/license.html>.
%
%-----------------------------------------------------------------------

% Diagnosis values

  values.success                      = int8(  0);	% solution found
  values.fail_on_argument             = int8(  1);	% an argument is wrong
  values.stop_on_itermax              = int8(  2);	% exit on the max number of iterations reached
  values.stop_on_dxmin                = int8(  3);	% exit on dxmin
  values.fail_with_quadprog           = int8(  4);	% failure with Quadprog
  values.fail_with_quadprog_iter      = int8(  5);	% failure with Quadprog (too many iterations)
  values.fail_with_quadprog_infeas    = int8(  6);	% failure with Quadprog (infeasible system)
  values.fail_on_nonnegative_slope    = int8(  7);	% exit on nonnegative slope of the least-square function
  values.stop_on_iterlsmax            = int8(  8);	% exit on too many linesearch trials
  values.stop_on_alphamin             = int8(  9);	% linesearch blocked on minimal stepsize (rounding error suspected)
  values.fail_with_index_computation  = int8( 10);	% exit due to a failure in the computation of index sets (hence this computation needs to be revised)
  values.fail_on_degeneracy           = int8( 11);	% a principal submatrix of M is singular
  values.fail_no_index_set_change     = int8( 12);	% too many times identical consecutuve index sets

% Running values

  values.fail_on_recursivity_before   = int8(101);	% recursivity fails before calling nmvar (no decrease in the number of complementarity constraints)
  values.fail_on_recursivity_after    = int8(102);	% recursivity fails after calling nmvar
  values.time_direction               = int8(103);	% CPU time of the computation of the directions
  values.success_tol                  = int8(104);	% stop on the given tolerance on theta
  values.success_same_polyhedron      = int8(105);	% stop on a displacement with unit stepsize in the same polyhedron
  values.hybrid                       = 0;		% hybrid direction

% Other values

  values.dline = '------------------------------------------------------------------------------------------';
  values.eline = '==========================================================================================';

% Return

return
