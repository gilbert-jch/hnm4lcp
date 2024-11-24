function [M,q,x,y,options,info] = hnm4lcp_prelim (M,q,x,options,values,info);

%%
% [M,q,x,y,options,info] = hnm4lcp_prelim (M,q,x,options,values,info);
%
% This function realizes the following preliminary tasks:
% - set options (default values if absent, numeric values for
%   lexical options)
% - check the given options
% - get the dimensions
% - initial printings

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

  y = [];

% Set fout and verb

  if ~isstruct(options)
    warning('hnm4lcp:OptionsNotStruct','argument ''options'' is expected to be a structure\n');
    info.flag = values.fail_on_argument;
    return
  end
  if isfield(options,'fout')
    if isempty(options.fout) | (options.fout <= 0) | (fix(options.fout)-options.fout)
      if options.verb, fprintf('\n### hnm4lcp_prelim: options.fout = %g not recongnized\n',options.fout); end
      info.flag = values.fail_on_argument;
      return
    end
  end
  if ~isfield(options,'fout')
    options.fout = 1;			% default output channel
  end
  fout = options.fout;

  if isfield(options,'verb')
    if isempty(options.verb)
      fprintf(fout,'\n### hnm4lcp_prelim: options.verb must be nonnegative\n');
      info.flag = values.fail_on_argument;
      return
    end
  else
    options.verb = 1;			% default verbosity level
  end
  verb = options.verb;			% to reduce time access

  n = length(q);
  if n <= 0
    fprintf(fout,'\n### hnm4lcp_prelim: the length of q must be positive (it is equal to %i)\n',n);
    info.flag = values.fail_on_argument;
    return
  end

% Check options

  % Check options fields

  possible_options = {'abstol', ...
                      'direction', ...
                      'dxmin', ...
                      'dymin', ...
                      'eta_descent', ...
                      'eta_inexact', ...
                      'fout', ...
                      'inexact', ...
                      'itermax', ...
                      'interpolation_safeguard', ...
                      'linesearch', ...
                      'linesearch_itermax', ...
                      'linesearch_m', ...
                      'qpmaxdim', ...
                      'qpmaxiter', ...
                      'qptol', ...
                      'reltol', ...
                      'scaling', ...
                      'stop_norm', ...
                      'tau', ...
                      'verb'};

  actual_options = fieldnames(options);
  for i=1:length(actual_options)
    if ~ismember(actual_options(i,:),possible_options)
      fprintf(fout,'\n### hnm4lcp_prelim: field ''%s'' of structure ''options'' is not recognized\n',char(actual_options(i)));
      info.flag = values.fail_on_argument;
      return
    end
  end

  if isfield(options,'itermax')
    if isempty(options.itermax) | (options.itermax <= 0)
      if verb; fprintf(fout,'\n### hnm4lcp_prelim: options.itermax must be > 0\n'); end
      info.flag = values.fail_on_argument;
      return
    end
  else
    options.itermax = 1000;			% default maximal number of iterations
  end

  if isfield(options,'dxmin')
    if isempty(options.dxmin) | (options.dxmin <= 0)
      if verb; fprintf(fout,'\n### hnm4lcp_prelim: options.dxmin must be > 0\n'); end
      info.flag = values.fail_on_argument;
      return
    end
  else
    options.dxmin = 1.e-10;			% default resolution in x
  end

  if isfield(options,'dymin')
    if isempty(options.dymin) | any(options.dymin <= 0)
      if verb
        fprintf(fout,'\n### hnm4lcp_prelim: options.dymin must be a positive vector of size %i\n',n);
      end
      info.flag = values.fail_on_argument;
      return
    end
  else
    options.dymin = 1.e-10*ones(n,1);		% default resolution in y
  end

  if isfield(options,'abstol') & isfield(options,'reltol')
    if verb; fprintf(fout,'\n### hnm4lcp_prelim: at most one of the tolerances options.abstol and options.reltol must be specified\n'); end
    info.flag = values.fail_on_argument;
    return
  end
  if isfield(options,'stop_norm')
    if ~ismember(options.stop_norm,[2,Inf])
      if verb; fprintf(fout,'\n### hnm4lcp_prelim: options.stop_norm can only be 2 or Inf.\n'); end
      info.flag = values.fail_on_argument;
      return
    end
  else
    options.stop_norm = 2;			% default l2 norm used in the stopping criterion
  end

  if isfield(options,'abstol')
    if isempty(options.abstol) | (options.abstol <= 0)
      if verb; fprintf(fout,'\n### hnm4lcp_prelim: options.abstol must be > 0\n'); end
      info.flag = values.fail_on_argument;
      return
    end
  elseif isfield(options,'reltol')
    if isempty(options.reltol) | (options.reltol <= 0) | (options.reltol > 1)
      if verb; fprintf(fout,'\n### hnm4lcp_prelim: options.reltol must be in (0,1]\n'); end
      info.flag = values.fail_on_argument;
      return
    end
    if isempty(x)
      if verb; fprintf(fout,'\n### hnm4lcp_prelim: x0 must be nonempty when the tolerance options.reltol is used\n'); end
      info.flag = values.fail_on_argument;
      return
    end
  else
    options.abstol = 1.e-8;			% default required relative precision
  end

  % scaling

  if isfield(options,'scaling')
    if isempty(options.scaling) | (fix(options.scaling) ~= options.scaling) | (options.scaling < 0) | (options.scaling > 4)
      if verb; fprintf(fout,'\n### hnm4lcp_prelim: options.scaling must be an integer between 0 and 3\n'); end
      info.flag = values.fail_on_argument;
      return
    end
  else
    options.scaling = 2;			% default diagonal scaling
  end

  % direction

  if isfield(options,'direction')
    if isempty(options.direction) | (fix(options.direction) ~= options.direction) | (options.direction < 0) | (options.direction > 1)
      if verb; fprintf(fout,'\n\n### nmvar: options.direction must be 0 or 1\n'); end
      info.flag = values.fail_on_argument;
      return
    end
  else
    options.direction = values.hybrid;		% default: use the hybrid direction
  end

  if isfield(options,'inexact') & ~isempty(options.inexact)
    if ~islogical(options.inexact)
      if verb; fprintf(fout,'\n\n### nmvar: options.inexact must be a logical\n'); end
      info.flag = values.fail_on_argument;
      return
    end
  else
    options.inexact = true;			% default: check whether the computed d is however an inexact PNM direction
  end

  if isfield(options,'projection')
    if isempty(options.projection) | (fix(options.projection) ~= options.projection) | (options.projection < 0) | (options.projection > 1)
      if verb; fprintf(fout,'\n\n### nmvar: options.projection must be 0 or 1\n'); end
      info.flag = values.fail_on_argument;
      return
    end
  else
    options.projection = 0;		% default: use the hybrid projection
  end

  % tau, qpmaxdim, qptol, eta_descent, eta_inexact

  if isfield(options,'tau') & ~isempty(options.tau)
    if options.tau < 0
      if verb; fprintf(fout,'\n### hnm4lcp_prelim: options.tau must be >= 0, possibly infinite\n'); end
      info.flag = values.fail_on_argument;
      return
    end
  else
    options.tau = 100*options.dxmin;
  end

  if isfield(options,'qpmaxdim')
    if isempty(options.qpmaxdim)
      options.qpmaxdim = Inf;
    elseif options.qpmaxdim < 0
      if verb; fprintf(fout,'\n### hnm4lcp_prelim: options.qpmaxdim must be >= 0, possibly infinite\n'); end
      info.flag = values.fail_on_argument;
      return
    elseif fix(options.qpmaxdim) ~= options.qpmaxdim
      options.qpmaxdim = max(1,fix(options.qpmaxdim));
      if verb; fprintf(fout,'\n### hnm4lcp_prelim: options.qpmaxdim must be integer, reset to %i\n',options.qpmaxdim); end
    end
  else
    options.qpmaxdim = Inf;			% defaul: no limitation of the number of inequalities of the QPs
  end

  if isfield(options,'qpmaxiter')
    if isempty(options.qpmaxiter)
      options.qpmaxiter = 100;
    elseif (options.qpmaxiter <= 0) | (fix(options.qpmaxiter) ~= options.qpmaxiter)
      if verb; fprintf(fout,'\n### hnm4lcp_prelim: options.qpmaxiter must be a positive integer, possibly infinite\n'); end
      info.flag = values.fail_on_argument;
      return
    end
  else
    options.qpmaxiter = 100;			% default max iterations for the QP solver
  end

  if isfield(options,'qptol') & ~isempty(options.qptol)
    if options.qptol <= 0
      if verb; fprintf(fout,'\n### hnm4lcp_prelim: options.qptol must be > 0\n'); end
      info.flag = values.fail_on_argument;
      return
    end
  else
    options.qptol = 1.e-7;
  end

  if isfield(options,'eta_descent')
    if isempty(options.eta_descent) | (options.eta_descent < 0) | (options.eta_descent >= 1)
      if verb; fprintf(fout,'\n### hnm4lcp_prelim: options.eta_descent must be in [0,1)\n'); end
      info.flag = values.fail_on_argument;
      return
    end
  else
    options.eta_descent = 0.9;			% eta in formula (2.19) in the paper
  end

  if isfield(options,'eta_inexact')
    if isempty(options.eta_inexact) | (options.eta_inexact < 0) | (options.eta_inexact >= 1)
      if verb; fprintf(fout,'\n### hnm4lcp_prelim: options.eta_inexact must be in [0,1)\n'); end
      info.flag = values.fail_on_argument;
      return
    end
  else
    options.eta_inexact = 0.9;			% eta in formula (2.22) in the paper
  end

  % Type of linesearch

  if isfield(options,'linesearch')
    if isempty(options.linesearch) | (fix(options.linesearch) ~= options.linesearch) | (options.linesearch < 0) | (options.linesearch > 2)
      if verb; fprintf(fout,'\n### hnm4lcp_prelim: ''options.linesearch'' must be in [0:2]\n'); end
      info.flag = values.fail_on_argument;
      return
    end
  else
    options.linesearch = 1;			% default Armijo linesearch with backtracking by stepsize halving
  end

  if isfield(options,'linesearch_itermax')
    if isempty(options.linesearch_itermax) || (fix(options.linesearch_itermax) ~= options.linesearch_itermax) || (options.linesearch_itermax < 1)
      if verb; fprintf(fout,'\n### hnm4lcp_prelim: options.linesearch_itermax (= %g) must be an integer >= 1\n',options.linesearch_itermax); end
      info.flag = values.fail_on_argument;
      return
    end
  else
    options.linesearch_itermax = 100;		% Defaul memory depth of the nonmonotone LS (hence 1 means monotone LS)
  end

  if isfield(options,'interpolation_safeguard')
    if isempty(options.interpolation_safeguard) | options.interpolation_safeguard <= 0 | options.interpolation_safeguard > 0.5
      if verb; fprintf(fout,'\n### hnm4lcp_prelim: options.interpolation_safeguard (= %g) must be in the interval (0,0.5]\n',options.interpolation_safeguard); end
      info.flag = values.fail_on_argument;
      return
    end
  else
    options.interpolation_safeguard = 0.1;	% Default interpolation safeguard
  end

  if isfield(options,'linesearch_m')
    if isempty(options.linesearch_m)
      options.linesearch_m = 1;
    elseif (fix(options.linesearch_m) ~= options.linesearch_m) | options.linesearch_m < 1
      if verb; fprintf(fout,'\n### hnm4lcp_prelim: options.linesearch_m (= %g) must be an integer >= 1\n',options.linesearch_m); end
      info.flag = values.fail_on_argument;
      return
    end
  else
    options.linesearch_m = 1;			% Defaul memory depth of the nonmonotone LS (hence 1 means monotone LS)
  end

% Initial x, y and options.abstol

  x1_is_empty = false;
  if isempty(x)
    x1_is_empty = true;
    x = zeros(n,1);
  end
  y = M*x+q;
  if ~isfield(options,'abstol')	% then, by the logic above, options.reltol is set
    options.abstol = options.reltol*abs(min(x,y));
  end

% Consistency of the parameters

  if options.projection == 1			% project the computed NM direction
    if options.direction ~= values.hybrid	% no NM direction is computed
      options.projection = 0;			% project 0 on the polyhedron, instead of a NM direction
      if verb
        fprintf(fout,'\n### hnm4lcp_prelim: options.projection == %i, while options.direction ~= 0, hence options.projection reset to 0\n', ...
          options.projection);
      end
    end
  end

  if options.tau < options.dymin
    if verb
      fprintf(fout,'\n### hnm4lcp_prelim: options.tau (=%11.5e) cannot be < options.dymin (=%11.5e), reset to %11.5e\n', ...
        options.tau,options.dymin,options.dymin);
    end
  end

% Print

  if verb
    fprintf(fout,'\n%s',values.eline);
    fprintf (fout,'\nHNM4LCP (LCP solver, version 1.0, November 2024)',n);
    fprintf (fout,'\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^',n);
    fprintf (fout,'\nOutput channel  = %i',fout);
    fprintf (fout,'\nVerbosity level = %i',verb);
    fprintf (fout,'\n\nProblem');
    if issparse(M)
      fprintf (fout,'\n. Sparse matrix M');
    else
      fprintf (fout,'\n. Dense matrix M');
    end
    fprintf (fout,'\n. Number of variables          = %6i',n);
    fprintf (fout,'\n. Resolution in x (dxmin)      =      %8.2e',options.dxmin);
    fprintf (fout,'\n. Resolution in y (dymin)      =      %8.2e(mean), %8.2e(var)',mean(options.dymin),var(options.dymin));
    if x1_is_empty
      fprintf (fout,'\n. x1 initialized to zero');
    else
      fprintf (fout,'\n. x1 taken on entry');
    end

    fprintf (fout,'\n\nMethod');
    if options.direction == values.hybrid
      fprintf (fout,'\n. Hybrid Newton-min direction');
    else
      fprintf (fout,'\n. Secure polyhedral Newton-min direction');
    end
    if options.projection == 0
      fprintf (fout,'\n. PNM direction obtained by projecting 0');
    else
      fprintf (fout,'\n. PNM direction obtained by projecting an NM direction');
    end
    switch options.linesearch
      case 0
        fprintf (fout,'\n. Unit stepsize');
      case 1
        fprintf (fout,'\n. Armijo linesearch, backtracking by stepsize halving');
      case 2
        fprintf (fout,'\n. Armijo linesearch, backtracking by stepsize interpolation (safeguard = %7.1e)',options.interpolation_safeguard);
      case 3
        fprintf (fout,'\n. Exact linesearch');
    end
    if options.linesearch_m > 1
      fprintf (fout,', nonmonotone memory depth = %i',options.linesearch_m);
    end
    fprintf (fout,'\n. Max stepsize trials per linesearch         = %6i',options.linesearch_itermax);
    fprintf (fout,'\n. Sufficient descent parameter (eta_descent) =      %8.2e',options.eta_descent);
    if options.inexact == true
      fprintf (fout,'\n. Inexactness parameter (eta_inexact)        =      %8.2e',options.eta_inexact);
    end
    fprintf (fout,'\n. Kink proximity (tau)                       = ');
    if options.tau == Inf
      fprintf (fout,'   Inf');
    else
      fprintf (fout,'     %8.2e',options.tau);
    end
    fprintf (fout,'\n. QP max dimension (qpmaxdim)                = %6i',options.qpmaxdim);
    fprintf (fout,'\n. QP max iter (qpmaxiter)                    = %6i',options.qpmaxiter);
    fprintf (fout,'\n. Max iterations (itermax)                   = %6i',options.itermax);
    fprintf (fout,'\n. Stopping norm (stop_norm)                  = %6i',options.stop_norm);
    fprintf (fout,'\n. Max iterations (itermax)                   = %6i',options.itermax);
    if isfield(options,'abstol')
      fprintf (fout,'\n. Absolute tolerance (abstol)                =      %8.2e',options.abstol);
    else
      fprintf (fout,'\n. Relative tolerance (reltol)                =      %8.2e',options.reltol);
    end

  end

  return
