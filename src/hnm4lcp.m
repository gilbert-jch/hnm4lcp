function [x,info] = hnm4lcp(M,q,x0,options)

%%
% [x] = hnm4lcp(M,q,x0,options)
%
% Hnm4lcp solves a standard linear complementarity problem (LCP), which
% is a problem of the form
%
%    0 <= x _|_ (M*x+q) >= 0,
%
% where M is a real matrix of size n, while x and q are vectors of
% length n. This means that a vector x is sought such that the
% components of x and M*x+q are nonnegative and x'*(M*x+q) = 0. The
% n-vector x0 is an optional initial vector for the iterative method and
% options is an optional structure of options (see its description
% below).
%
% The matrix M is assumed to be nondegenerate (i.e., M(I,I) is
% nonsingular for all subset I of [1:n]) and no provision is made to
% deal with cases where this assumption is not verified. In case of
% degeneracy is encountered, HNM4LCP fails and mentions the problem (see
% info.flag below).
%
% HNM4LCP implements the HNM algirthms described in the paper
%
%    Jean-Pierre Dussault, Mathieu Frappier, Jean Charles Gilbert
%    (2024), 'Polyhedral Newton-min algorithms for complementarity
%    problems', Mathematical Programming, in revision [hal-02306526],
%
% referred to as the DFG paper below. The globalization of the
% algorithms is made by using the least-squares merit function, which is
% the map Theta defined at x in Rn by
%
%    Theta(x) = 0.5*norm(H(x))^2,
%
% where
%
%    H(x) = min(x,M*x+q),
%
% componentwise.
%
% HNM4LCP is a research code, with many options that are not necessary
% easy to understand. If an option is not clear, just ignore it and its
% default value will be set. Here is the list of options and their
% meaning.
%
% . options.abstol (not compatible with 'options.reltol'): an iterate x
%     is considered to be a solution when the norm of H(x) is less than
%     options.abstol; the norm is the l2 (Euclidean) norm if
%     options.stop_norm == 2 and the max norm otherwise;
% . options.direction (default 0): specifies the type of search
%     direction the algorithm must use:
%     = 0: hybrid Newton-min direction (algorithm HNM or 3.8 in the DFG
%          paper, recommended),
%     = 1: (obsolete) secure polyhedral Newton-min direction (algorithm
%          PNM or 2.11 in the DFG paper, less efficient),
% . options.dxmin (default 1.e-10): positive number giving the
%     resolution in x; if the increment in x to improve the solution
%     approximation is less than options.dxmin in the sup-norm, the
%     algorithm considers that it cannot do better and stops;
% . options.dymin (default 1.e-10*ones(n,1)): scalar or vector of
%     positive numbers giving the resolution in y = M*x+q; if abs(y-yy)
%     < options.dymin, the algorithm considers that y and yy are the
%     same vectors; if options.dymin is set as a scalar, it has the same
%     effect as options.dymin*ones(n,1); when dymin is a vector of non
%     indentical components, it introduces a kind of preconditioning;
%     this value is used for determininng index sets (hence has a major
%     impact, in particular for determining the NM direction), for
%     preventing options.tau from being smaller than dymin, for checking
%     some inexact direction inequalities and for printings;
% . options.eta_descent (default 0.9): parameter in [0,1) that specifies
%     the parameter eta in the sufficient descent condition (2.19) of
%     the DFG paper; it is used in the hybrid algorithm to test whether
%     the NM direction has the sufficient descent property (2.19);
% . options.eta_inexact (default 0.9): parameter in [0,1) that specifies
%     the parameter eta in the inexactness system (2.22) of the DFG
%     paper; it is used when the QP solver declares that its computed
%     direction d is not a solution (usually because of rounding errors
%     or too demanding stopping criteria); set eta_inexact large to
%     accept that approximate QP solution more often;
% . options.fout: output channel for the outputs (0 for the screen, fout
%     = fopen('res.txt','w') for writing on the file 'res.txt');
% . options.hybrid: if true, use the hybrid algorithm, i.e., try first
%     the NM direction;
% . options.inexact (default true): specify whether a direction computed
%     by the QP solver, but that does not satisfy the requested
%     precision, must be accepted when it is an inexact PNM direction;
% . options.itermax: maximum number of iterations allowed;
% . options.interpolation_safeguard (default 1.e-2): a value in (0,0.5];
%     when options.linesearch == 2, the stepsize found by quadratic
%     interpolation is safeguarded in the interval [isg*alpha,
%     (1-isg)*alpha], where alpha is the current stepsize and isg =
%     options.interpolation_safeguard;
% . options.linesearch (obsolete, must be set to 1):
%     = 0: unit stepsize;
%     = 1: with Armijo linesearch and backtracking by stepsize halving;
%     = 2: with Armijo linesearch and backtracking by stepsize
%          interpolation;
%     = 3: with exact linesearch (disactivated);
% . options.linesearch_m (default = 1): for nonmonotone LS, set
%     linesearch_m > 1, which then specifies the depth of the previous
%     merit function vallues to memorize;
% . options.linesearch_itermax (default = 100): maximum number of
%     stepsize trials per iteration during the linesearch; this number
%     may be large when the problem is ill-conditioned;
% . options.qpmaxdim: positive integer used to limit the number of
%     variables of the QP to solve; the dimension of the QP is never
%     smaller than the number of indices i such that x(i) = y(i) < 0; if
%     this number is larger than qpmaxdim, the dimension of the QP is
%     limited to qpmaxdim; the selected indices are those in {i:
%     |x(i)-y(i)| < tau, x(i)<0, y(i)<0} having the smallest value of
%     |x(i)-y(i)|; default value is Inf (i.e., the option is inactive);
% . options.qpmaxiter (defaut 100): max number of iterations of the QP
%     solver;
% . options.qptol (default 1.e-7): tolerance used to set the required
%     precision of the QP solver;
% . options.reltol (not compatible with 'options.abstol'): an iterate x
%     is considered to be a solution when the value Theta(x) is less
%     than options.reltol*Theta(x0);
% . options.scaling
% . info.linesearches: number of iterations with linesearch;
% . options.stop_norm specifies the norm used in the stopping criterion
%     (see options.abstol); this one is the Euclidean norm if
%     options.stop_norm == 2 and the max norm otherwise;
% . options.tau (default 100*options.dxmin): real number measuring the
%     `negative kink' proximity; this value is crucial in the
%     determination of the search directions;
% . options.verb: verbosity level of the output
%     = 0: no output; works silently,
%     = 1: one line per iteration,
%     = 2: one paragraph per iteration, including linesearch
%          information,
%     = 3: more output, requiring extra-computation, not essemtial for
%          finding a solution.
%
% On return
% . x: solution to the LCP,
% . info.etime: elapsed time in hnm4lcp
% . info.etime_ws: elapsed time in hnm4lcp, descardong the time spent in
%   sclaing the problem
% . info.flag: diagonis
%     =  0: solution found,
%     =  1: an input argument is wrong,
%     =  2: max number of iterations reached,
%     =  3: stop on dxmin,
%     =  4: failure of the Harker-Pang technique (obsolete),
%     =  5: failure with Quadprog (too many iterations),
%     =  6: failure with Quadprog (infeasible system),
%     =  7: positive slope of the least-square function,
%     =  8: too many linesearch trials,
%     =  9: linesearch blocked at minimal stepsize (rounding error
%           suspected),
%     = 10: index set computation failure;
%     = 11: M is degenerate (i.e., a principal submatrix of M is
%           singular);
%     = 12: too many times identical consecutuve index sets;
% . info.iter: number of iterations;
% . info.qpsolve: number of QP solve;
% . info.qpnvar: sum of the dimensions of the solved QPs.
% . info.stepsize: median stepsize

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

%------------------------
% Initialization
%------------------------

  tic_init = tic;	% Launch the timer

  [values] = hnm4lcp_values();	% Set the constant values used in HNM4LCP

  no_index_set_change_ctr =  0;	% counter of the warning "no change in index sets"
  no_index_set_change_max =  5;	% max number of warning "no change in index sets"

  info.flag   = values.success;
  info.qpnvar = 0;

  %------------------------
  % Check input arguments
  %------------------------

  if nargin < 4, options = struct();
    if nargin < 3, x0 = [];
      if nargin < 2
        fprintf('\n### hnm4lcp: the first 2 arguments are required\n\n');
        x = x0;
        info.flag = values.fail_on_argument;
        info.timer(values.cput_total) = toc(tic_init);
        return
      end
    end
  end

  q = q(:);
  x = x0(:);
  infos = [];

  %------------------------
  % Decode options
  %------------------------

  if ~exist('options') | isempty(options); options = struct(); end

  [M,q,x,y,options,info] = hnm4lcp_prelim(M,q,x,options,values,info);
  if info.flag; return; end

  n = length(q);

  fout = options.fout;
  verb = options.verb;

  %------------------------
  % Initialization
  %------------------------

  info.nb_unit_stepsizes = 0;
  ls_m = options.linesearch_m;	% memory depth of the nonmonotone LS (hence 1 means monotone LS)
  ls_p = 0;			% pointer in ls_v and ls_i (increased by 1 before memorization)
  ls_v = -Inf*ones(ls_m,1);	% hence the max can always be taken on ls_v(1:ls_m), even though nothing has been memorized in ls_v(ls_p+1:ls_m)
  ls_i = zeros(ls_m,1);		% iteration number corresponding toe the max value in ls_v

  %------------------------
  % Scaling
  %------------------------

  tic_scaling = tic;	% Launch the timer for scaling

  if options.scaling == 0

    if verb; fprintf (fout,'\n. No scaling'); end

  elseif options.scaling == 1
   
    sigma = 1/norm(M,'fro');
    M = sigma*M;
    q = sigma*q;
    if verb
      fprintf (fout,'\n\nScaling M and q');
      fprintf (fout,'\n. Scalar scaling factor = %8.2e',sigma);
    end
    y = M*x+q;

  elseif options.scaling == 2

    fprintf (fout,'\n\nScaling M and q');
    [M,q,dl] = hnm4lcp_lscale (M,q,options);
    y = M*x+q;

  elseif options.scaling == 3

    fprintf (fout,'\n\nScaling M and q');

    [M,q,dl,dr] = hnm4lcp_biscale (M,q,2,1.e-1,options);

    x = x./dr;
    y = M*x+q;

  elseif options.scaling == 4

    S = zeros(n,1);
    for i = 1:n
      s = 1/norm(M(i,:));
      M(i,:) = s*M(i,:);
      S(i) = s;
    end
    q = S.*q;
    sigma  = sum(S)/n;
    sigma2 = norm(q);
    if sigma2 ~= 0
      q = q/sigma2;
    end
    fprintf (fout,'\n. Diagonal scaling of M and q by      %8.2e (average value)',sigma);
    if sigma2 ~= 0; fprintf (fout,'\n. Additional scaling of q by          %8.2e',1/sigma2); end
    y = M*x+q;

  end

  scaling_time = toc(tic_scaling);

  %------------------------
  % Parameters
  %------------------------

  dxmin = options.dxmin;
  dymin = options.dymin;

  ls_itermax      = options.linesearch_itermax;	% max iter in LS
  ls_step_factor  = 0.5;			% stepsize factor in backtraking without interpolation
  ls_step_factor9 = ls_step_factor^9;		% stepsize after 9 backtraking
  ls_isg          = options.interpolation_safeguard;	% interpolation safeguard (must be in (0,0.5]), the new alpha will be in [ls_isg*t, (1-ls_isg)*t]
  ls_omega        = 1.e-4;			% fraction of the slope in the Armijo condition
  itblocksize     = 500;			% block size for reallocation of stepsizes
  stepsizes       = zeros(itblocksize,1);

  qp_options = optimoptions ('Quadprog', ...
                             'Display'             , 'off', ...
                             'Diagnostics'         , 'off', ...
                             'MaxIter'             , options.qpmaxiter, ...
                             'OptimalityTolerance' , options.qptol, ...
                             'Algorithm'           , 'interior-point-convex');
%                            'Algorithm'           , 'trust-region-reflective');
%                            'Algorithm'           , 'active-set');
%                            'Algorithm'           , 'interior-point-convex');

%------------------------
% Loop
%------------------------

  info.qpsolve = 0;
  iter         = 0;

  Anm = [];
  Inm = [];
  A   = [];
  I   = [];
  Eyp = [];
  Emt = [];

  tol_merit = options.abstol*options.abstol;	% options.abstol has been defined in hnm4lcp_prelim

  alpha = [];

  % It is assumed that the following values are known at the begining of the loop: x, y = M*x+q

  while true

    %------------------------
    % Iteration counter
    %------------------------

    if iter >= options.itermax
      info.flag = values.stop_on_itermax;
      break
    end
    iter = iter+1;

    %------------------------
    % Print
    %------------------------

    H = min(x,y);
    merit = 0.5*(H'*H);

    if verb
      if verb == 1
        if mod(iter,50) == 1
          fprintf (fout,'\n%s\niter  |min(x,M*x+q)|   |A|   |I| |Emt|  QP',values.dline);
          if options.linesearch; fprintf (fout,'  stepsize'); end
          fprintf (fout,'\n%s',values.dline);
        end
        fprintf (fout,'\n%4i  %14.8e',iter,sqrt(2*merit));
      elseif verb > 1
        fprintf (fout,'\n%s',values.dline);
        fprintf (fout,'\niter = %0i, merit = %12.5e, |min(x,M*x+q)| = %11.5e(l2) ou %11.5e(max)',iter,merit,sqrt(2*merit),max(abs(H)));
        if verb > 3
          fprintf('\nx ='); fprintf(' %12.5e',x);
          fprintf('\ny ='); fprintf(' %12.5e',y);
        end
      end
    end

    %------------------------
    % Stopping criterion: on the norm of min(x,M*x+q)
    %------------------------

    if options.stop_norm == 2
      if sqrt(2*merit) <= options.abstol	% convergence on the l2-norm
        iter = iter-1;
        info.flag = values.success;
        success = values.success_tol;
        break
      end
    else
      if max(abs(H)) <= options.abstol		% convergence on the max-norm
        iter = iter-1;
        info.flag = values.success;
        success = values.success_tol;
        break
      end
    end

    %------------------------
    % Index sets
    %------------------------

    absxmy = abs(x-y);

    % index sets for computing the slope along the search direction and detection of their change

    if options.linesearch

      Esl = find(absxmy <= dymin);
      Asl = setdiff( find(x < y), Esl );
      Isl = setdiff( find(x > y), Esl );

      if verb > 1; fprintf (fout,'\nIndex sets for slope: |A| = %i, |I| = %i, |E| = %i',length(Asl),length(Isl),length(Esl)); end

    end

    % index sets for NM

    if options.direction == values.hybrid

      Anm_old = Anm;
      Inm_old = Inm;

      Anm = find(x <= y+dymin);
      Inm = setdiff([1:n],Anm);

      if verb > 1; fprintf (fout,'\nIndex sets for NM:    |A| = %i, |I| = %i',length(Anm),length(Inm)); end

      if (length(setdiff(Anm,Anm_old)) == 0) && (length(setdiff(Anm_old,Anm)) == 0) && (length(setdiff(Inm,Inm_old)) == 0) && (length(setdiff(Inm_old,Inm)) == 0)
        no_index_set_change_ctr = no_index_set_change_ctr+1;
        if verb > 1; fprintf (fout,'   ### (hnm4lcp) Warning %i: no change',no_index_set_change_ctr); end
        if no_index_set_change_ctr > no_index_set_change_max
          info.flag = values.fail_no_index_set_change;
          if verb; fprintf(fout,'\n\n### Too many times identical consecutuve index sets\n'); end
          break
        end
      else
        no_index_set_change_ctr = 0;
      end

    end

    % index sets for QP

    Emt = find( (absxmy < max(dymin,options.tau)) & (x < 0) & (y < 0) );
    if length(Emt) > options.qpmaxdim
      [~,I] = sort(absxmy(Emt));
      Emt = Emt(I(1:options.qpmaxdim));
    end
    A = setdiff( find(x<=y+dymin), Emt);
    I = setdiff( find(x> y+dymin), Emt);

    if verb == 1
      fprintf (fout,' %5i %5i %4i ',length(A),length(I),length(Emt));
    elseif verb == 2
      fprintf (fout,'\nIndex sets for QP:    |A| = %i, |I| = %i, |Emt| = %i ',length(A),length(I),length(Emt));
    elseif verb > 2
      fprintf (fout,'\nIndex sets for QP:');
      if length(A)
        fprintf (fout,'\n. A =');
        fprintf (fout,' %0i',A);
      end
      if length(I)
        fprintf (fout,'\n. I =');
        fprintf (fout,' %0i',I);
      end
      if length(Emt)
        fprintf (fout,'\n. Emtl  =');
        fprintf (fout,' %0i',Emt);
      end
    end

    if length(A)+length(I)+length(Emt) ~= n
      info.flag = values.fail_with_index_computation;
      if verb == 1
        fprintf(fout,'\n\n### Failure in the computation of index sets\n');
      elseif verb > 1
        fprintf(fout,'\nFailure with the computation of index sets:\n. x   ='); fprintf(fout,'  %12.5e',x);
        fprintf(fout,'\n. y   ='); fprintf(fout,'  %12.5e',y);
        fprintf(fout,'\n. A   ='); fprintf(fout,'  %0i',A);
        fprintf(fout,'\n. I   ='); fprintf(fout,'  %0i',I);
        fprintf(fout,'\n. Emt ='); fprintf(fout,'  %0i',Emt);
      end
      break
    end

    %------------------------
    % Linesearch preparation (useful in hnm4lcp_direction)
    %------------------------

    if options.linesearch

      % Nonmonotine linesearch material (also appropriate for monotone linesearch): ls_p is a pointer that cycles in the table
      % of theta-values ls_v

      ls_p = ls_p+1;
      if ls_p > ls_m
        ls_p = 1;
      end

      ls_v(ls_p) = merit;
      ls_i(ls_p) = iter;

      [ls_vmax,imax] = max(ls_v);
      ls_imax        = ls_i(imax(1));

      % Compute the slope of theta at the current point and armijo_slope (the negative term of the Armijo condition)

      armijo_slope = -2*ls_omega*(1-options.eta_descent)*merit;	

    end

    %------------------------
    % Direction
    %------------------------

    [d,info] = hnm4lcp_direction (M,q,x,y,Anm,Inm,A,I,Emt,H,merit,ls_vmax,armijo_slope,qp_options,options,values,info);

    if info.take_unit_stepsize	% one can take a unit stepsize along the computed direction
      if(iter > 1) & (mod(iter,itblocksize) == 1)
        stepsizes(end+itblocksize) = 0;	% reallocate 'itblocksize' elements in stepsizes
      end
      stepsizes(iter) = 1;
      info.nb_unit_stepsizes = info.nb_unit_stepsizes+1;
      if verb == 1
        fprintf(fout,'      %8.2e',1);
      end
      x = x+d;
      y = M*x+q;
      continue;			% take the next iteration
    end

    if isempty(d)
      if info.flag
        if verb == 1
          fprintf (fout,'\n\n### Quadprog failed to find a solution to the QP');
          if info.flag == values.fail_with_quadprog_iter; fprintf (fout,' (too many iterations)'); end
          if info.flag == values.fail_with_quadprog_infeas; fprintf (fout,' (infeasible problem)'); end
          fprintf (fout,'\n');
        end
      end
      break
    end

    if verb == 1
      if info.aqpissolved
        if info.qpsave
          fprintf (fout,'  *');
        else
          fprintf (fout,'   ');
        end
        fprintf (fout,'1');
      else
        fprintf (fout,'    ');
      end
    end

    %------------------------
    % Linesearch
    %------------------------

    if options.linesearch

      % find the smallest (alphamin) and largest (alphamax) stepsizes in (0,1) leading to a kink of H (may not be a kink of
      % theta), if any (alphamin = 0 otherwise); alphamax is not really used actually

      ymx  = y-x;
      Md   = M*d;
      dmMd = d-Md;
      Ind  = find(dmMd~=0);		% indices in [1:n]
      alphak = ymx(Ind)./dmMd(Ind);
      [alphak,ik] = sort(alphak(find((alphak*norm(d)>dxmin)&(alphak<=1))));	% stepsises giving a kink of H (may not be a kink of theta)
      alphamin = 0;
      alphamax = 0;
      if length(alphak) > 0
        alphamin = alphak(1);
        alphamax = alphak(end);
      end

      % slope theta'(x;d)

      slope = x(Asl)'*d(Asl) + x(Esl)'*min(d(Esl),Md(Esl)) + y(Isl)'*Md(Isl);

      if slope >= 0
        info.flag  = values.fail_on_nonnegative_slope;
        break
      end

      % set the initial stepsize

      if (no_index_set_change_ctr == no_index_set_change_max) && (alphamin > 0)
        alpha = alphamin;
        alpha_set_to_alphamin = true;
      else
        alpha = 1;
        alpha_set_to_alphamin = false;
      end

      % preliminaries

      dn = norm(d,Inf);
      xp = x+alpha*d;
      yp = M*xp+q;
      iter_ls = 0;
      if verb > 1
        fprintf (fout,'\n\nLinesearch')
        fprintf (fout,'\n  min/max stepsizes = %11.5e/%11.5e, |d| = %11.5e, slope = %12.5e',alphamin,alphamax,norm(d),slope);
        if slope >= 0; fprintf (fout,'\n  ### warning: positive computed slope, try anyway'); end
        fprintf (fout,'\n  reference value      %11.5e (iteration %i)',ls_vmax,ls_imax);
        fprintf (fout,'\n  merit value          %11.5e',merit);
        if (no_index_set_change_ctr == no_index_set_change_max) && (alphamin > 0)
          fprintf(fout,'\n  too many times identical consecutuve index sets => move to next kink');
        end
        fprintf (fout,'\n  iter    stepsize     refer decr    slope estim');
      end

      if (no_index_set_change_ctr == no_index_set_change_max) && (alphamin > 0)
        if (iter > 1) & (mod(iter,itblocksize) == 1)
          stepsizes = [stepsizes;zeros(itblocksize,1)];	% reallocate 'itblocksize' elements in stepsizes
        end
        Hp = min(xp,yp);
        meritp = 0.5*(Hp'*Hp);
        if verb == 1; fprintf (fout,'  %8.2e',alpha); end
        if verb > 1; fprintf (fout,'\n  %4i  %12.6e  %12.5e  %12.5e',1,alpha,meritp-ls_vmax,(meritp-merit)/alpha); end
        stepsizes(iter) = alpha;
        x = xp;
        y = yp;
	continue
      end

      % linesearch loop (it is supposed that the next trial point (xp,yp) is known, when exiting the loop the new point is (x,y))

      alphamindn_dxmin = alphamin*dn-dxmin;

      while true

        if alpha*dn < dxmin
          info.flag  = values.stop_on_dxmin;
          break;
        end

        iter_ls = iter_ls+1;
        if iter_ls > ls_itermax
          info.flag = values.stop_on_iterlsmax;
          if verb > 1; fprintf (fout,'\n\n### too many linesearch trials (%i)\n',ls_itermax); end
          break
        end

        Hp = min(xp,yp);
        meritp = 0.5*(Hp'*Hp);

        if verb > 1; fprintf (fout,'\n  %4i  %12.6e  %12.5e  %12.5e',iter_ls,alpha,meritp-ls_vmax,(meritp-merit)/alpha); end

        % accept the stepsie alpha
        % - either if the Armijo condition is satisfied
        % - or ...

        exit_ls = false;

        if meritp <= ls_vmax + alpha*armijo_slope
          exit_ls = true;
        elseif (alphamin < 1) && (alpha < alphamin) && (meritp <= (1+1.e-1)*ls_vmax)
          exit_ls = true;
          if verb > 1; fprintf (fout,'\n  ### rounding errors: stepsize set to its min value %12.6e, with merit function limited to %11.5e',alphamin,(1+1.e-1)*ls_vmax); end
          alpha = alphamin;
          xp = x + alpha*d;
          yp = M*xp+q;
          Hp = min(xp,yp);
          meritp = 0.5*Hp'*Hp;
        end

        if exit_ls
          if (iter > 1) & (mod(iter,itblocksize) == 1)
            stepsizes = [stepsizes;zeros(itblocksize,1)];	% reallocate 'itblocksize' elements in stepsizes
          end
          stepsizes(iter) = alpha;
          x = xp;
          y = yp;
          if verb == 1; fprintf (fout,'  %8.2e',alpha); end
          break;
        end

        % new stepsize alpha

        alpha0 = alpha;

        if options.linesearch == 1	% stepsize halving

          if (meritp == Inf) && (alpha > alphamin) && (alphamindn_dxmin > 0)	% after 10 trials, set alpha to alphamin
            alpha = alphamin;
            alpha_set_to_alphamin = true;
          else
            alpha = alpha*ls_step_factor;
          end

        elseif options.linesearch == 2	% quadratic interpolation (the quadratic model is t -> a*t^2 + b*t + c)

          a = (meritp-merit-slope*alpha)/alpha^2;
          alpham = -slope/(2*a);	% minimizer of the quadratic model
          if alpham < ls_isg*alpha
            alpha = ls_isg*alpha;
          elseif alpham > (1.-ls_isg)*alpha
            alpha = (1-ls_isg)*alpha;
          else
            alpha = alpham;
          end

        end

        if (alpha0 > alphamin) && (alpha < alphamin) && ~alpha_set_to_alphamin
          alpha = alphamin;
          alpha_set_to_alphamin = true;
        end

        xp = x+alpha*d;
        yp = M*xp+q;

      end	% of linesearch loop

      % count the number of consecutive unit stepsizes (alpha == 1 when iter_ls == 1)

      if alpha == 1
        info.nb_unit_stepsizes = info.nb_unit_stepsizes+1;
      else
        info.nb_unit_stepsizes = 0;
      end

      if (info.flag == values.stop_on_alphamin) | (info.flag == values.stop_on_iterlsmax) | (info.flag == values.stop_on_dxmin)
        break;
      end

    else	% unit stepsize

      alpha = 1;
      x = x+d;
      y = M*x+q;

    end

  end		% of iteration loop

%------------------------
% Various tunings
%------------------------

  info.iter  = iter;

  % descaling of x

  if options.scaling == 3
    x = x.*dr;
  end

  % Mean/median stepsize

  stepsizes     = stepsizes(1:info.iter);	% purge stepsizes from its final zero values
  stepsizes     = stepsizes(stepsizes~=1);	% only get the stepsizes < 1
  log2stepsizes = log2(stepsizes);

  info.linesearches        = length(stepsizes);
  info.stepsize_mean       = [];
  info.stepsize_log2mean   = [];
  info.stepsize_median     = [];
  info.stepsize_log2median = [];
  if info.linesearches
    info.stepsize_mean       = mean(stepsizes);
    info.stepsize_log2mean   = 2^mean(log2stepsizes);
    info.stepsize_median     = median(stepsizes);
    info.stepsize_log2median = 2^median(log2stepsizes);
  end

%------------------------
% Final printings
%------------------------

  if verb < 0
    if info.flag == values.fail_on_recursivity_before
      fprintf (fout,'\n\n### Recursivity cannot be carried on since |E| = n = %i\n',n);
    end
  end

  if verb > 0

    fprintf(fout,'\n%s',values.dline);

    if verb
      switch info.flag
        case values.success
          fprintf (fout,'\nSolution found');
          if success == values.success_tol
            fprintf (fout,' (merit function smaller than tolerance)');
          else
            fprintf (fout,' (unit stepsize in the same polyhedron)');
          end
          fprintf (fout,'\n%s\nNumber of iterations = %0i\n|min(x,M*x+q)|_2 = %11.5e',values.dline,info.iter,sqrt(2*merit));
        case values.stop_on_itermax
          fprintf (fout,'\n### Stopping on the maximum number of iterations allowed\n%s',values.dline);
        case values.stop_on_dxmin
          fprintf (fout,'\n### Too small change in x (= %8.2e < dxmin = %8.2e, alpha = %8.2e)\n%s',alpha*dn,dxmin,alpha,values.dline);
        case values.fail_with_quadprog
          % print already done
        case values.fail_on_recursivity_after
          fprintf(fout,'\n### Recursive call to hnm4lcp failled, flag on output = %i\n%s',hnm4lcp_info.flag,values.dline);
        case values.fail_on_nonnegative_slope
          fprintf (fout,'\n### Nonnegative slope (%11.5e) of the least-square function\n%s',slope,values.dline);
        case values.fail_on_degeneracy
          fprintf('\n### hnm4lcp: M is degenerate\n\n');
      end
    end

    if info.flag == values.success

      % Print

      if verb

        E = find(abs(x-y) < dymin);
        A = setdiff( find(x < y), E);
        I = setdiff( find(x > y), E);
        fprintf (fout,'\n|A| = %0i, |I| = %0i, |E| = %0i',length(A),length(I),length(E));

        fprintf(fout,'\nNumber of QPs = %i',info.qpsolve);
        if info.qpsolve
          fprintf(fout,'\nMean dimension of the QPs = %.1f',info.qpnvar/info.qpsolve);
        end
        if options.linesearch
          fprintf(fout,'\nNumber of iterations with linesearch = %i',info.linesearches);
          fprintf(fout,'\nStepsizes: mean = %11.5e (log2mean = %11.5e), median = %11.5e (log2median = %11.5e)', ...
            info.stepsize_mean,info.stepsize_log2mean,info.stepsize_median,info.stepsize_log2median);
          fprintf(fout,'\nNumber of final unit stepsizes = %0i',info.nb_unit_stepsizes);
        end
        info.etime    = toc(tic_init);
        fprintf(fout,'\nCPU time = %g sec',info.etime);
        if options.scaling
          info.etime_ws = info.etime - scaling_time;
          fprintf(fout,' (%g sec without scaling)',info.etime_ws);
        end

      end

    end

    fprintf(fout,'\n%s\n',values.eline);

  end

  return
