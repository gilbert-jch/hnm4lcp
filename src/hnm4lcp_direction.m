function [d,info] = hnm4lcp_direction (M,q,x,y,Anm,Inm,A,I,Emt,H,merit,max_merit,armijo_slope,qp_options,options,values,info)

%%
% [d,info] = hnm4lcp_direction (M,q,x,y,A,I,Emt,H,merit,max_merit,armijo_slope,qp_options,options,values,info)
%
% Computes a descent direction d at x for the PNM of HNM algorithms,
% which solve the linear complementarity problem 0 <= x _|_ y >= 0, with
% y = M*x+q.
%
% The direction d computed by this procedure may be a solution to tthe
% linear system
%
%   x(AE) + d(AE)    = 0
%   y(I)  + M(I,:)*d = 0
%
% where I = find(x>y) and AE = find(x<=y) or the minimum l2-norm
% direction satisfying
%
%   x(A) + d(A)          = 0,
%   y(I) + M(I,:)*d      = 0,
%   x(Emt) + d(Emt)     >= 0,
%   y(Emt) + M(Emt,:)*d >= 0,
%
% where the index sets A, I, and Emt are given on entry. Call first the
% NM direction and next the QP direction.
%
% The logic for choosing one or the orther is the following:
% - if length(Emt) == 0, the NM solution is computed,
% - otherwise
%   = if options.direction == values.hybrid, the NM direction is
%     first computed
%     + if this one is a sufficient descent direction for the
%       least-square merit function, it is adopted,
%     + otherwise the QP direction is computed,
%   = otherwise the QP direction is computed.
%
% It is assumed that the matrix M is nondegenerate, meaning that all its
% principal submatrices are nonsingular.
%
% If the direction cannot be computed, d is empty on output and the
% reason why this occurs is given in info.flag.
%   = values.fail_with_quadprog: other failures of quadprog than the two
%     below;
%   = values.fail_with_quadprog_infeas: there is no d satisfying the
%     conditions above, according to quadprog;
%   = values.fail_with_quadprog_iter: convergence does not occur with
%     the prescribed number of iterations.

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

  n    = length(x);
  nA   = length(A);
  nI   = length(I);
  nEmt = length(Emt);

  dymin = options.dymin;
  etad  = options.eta_descent;
  etai  = options.eta_inexact;
  tau   = options.tau;

  verb = options.verb;
  fout = options.fout;

  info.aqpissolved = false;
  info.qpsave      = false;

%------------------------
% For the hybrid algorithm, compute the NM direction and check that it
% is convenient 
%------------------------

  info.take_unit_stepsize = false;

  if options.direction == values.hybrid			% see whether solving a single linear system is enough

    if verb > 1; fprintf (fout,'\nDirection computation: hybrid direction\n. the NM direction is computed'); end

    % compute the NM direction

    nInm = length(Inm);

    if nInm
      if nInm < n
        d  = zeros(n,1);
        d(Anm) = -x(Anm);					% d(A) = -x(A)
        [L,U,P] = lu(M(Inm,Inm));				% then P*MII = L*U
        try
          d(Inm) = U\(L\(P*(M(Inm,Anm)*x(Anm)-y(Inm))));	% MII*d(I) = MIA*x(A)-y(I)
        catch
          info.flag = values.fail_on_degeneracy;
          return
        end
      else
        [L,U,P] = lu(M);				% then P*M = L*U
        try
          d = -(U\(L\(P*y)));				% M*d = -y
        catch
          info.flag = values.fail_on_degeneracy;
          return
        end
      end
    else
      d = -x;
    end

    % check whether the unit stepsize is acccepted along the NM direction

    xp = x+d;
    Hp = min(xp,M*xp+q);
    meritp = 0.5*(Hp'*Hp);
    if meritp <= max_merit + armijo_slope
      info.take_unit_stepsize = true;
      if verb > 1
        fprintf (fout,' and accepted with unit stepsize (|d| = %12.5e)',norm(d));
        fprintf (fout,'\n  reference merit value = %12.5e',max_merit);
        fprintf (fout,'\n  current merit value   = %12.5e',merit);
        fprintf (fout,'\n  new merit value       = %12.5e',meritp);
      end
      return
    end

    % check whether the direction is a sufficient descent direction, in the sense that it verifies condition (2.19b) in the DFG
    % paper

%   secure_direction = true;

    rho = zeros(nEmt,1);
    rho = max((x(Emt)+d(Emt))./x(Emt),(y(Emt)+M(Emt,:)*d)./y(Emt));
    if (rho.*H(Emt))'*H(Emt) <= etad*merit
      if verb > 1
        fprintf (fout,' and used in linesearch (|d| = %12.5e)',norm(d))
      end
      return
    end

    if verb > 1
      fprintf (fout,' but not a sufficient descent direction (|d| = %12.5e)',norm(d));
    end

  end

%------------------------
% Compute the PNM direction
%------------------------

  if nEmt	% a QP is solved in order to force the satisfaction of inequality constraints

    info.aqpissolved = true;
    info.qpsolve     = info.qpsolve+1;	% total number of QP solves
    info.qpnvar      = info.qpnvar+nEmt;	% sum of the number of variables of the QP's

    if verb > 1
      if options.direction ~= values.hybrid fprintf (fout,'\nDirection computation:'); end
      fprintf (fout,'\n. a QP to solve with %i variable',nEmt);
      if nEmt > 1; fprintf(fout,'s'); end
    end

    d = zeros(n,1);

    if nI
      [L,U,P]   = lu(M(I,I));			% P*MII = L*U
      try
        MIImMIEmt = U\(L\(P*M(I,Emt)));		% MII*MIImMIEmt = M(I,Emt);
      catch
        info.flag = values.fail_on_degeneracy;
        return
      end
      Ms        = M(Emt,I)*MIImMIEmt - M(Emt,Emt);	% negative Schur complement or Jacobian of the inequality constraints
      Hq        = eye(nEmt)+MIImMIEmt'*MIImMIEmt;	% Hessian of the objective of the quadratic problem
      if nA
        v = U\(L\(P*(y(I)-M(I,A)*x(A))));	% MII*v = y(I)-MIA*x(A)
        gq = MIImMIEmt'*v;			% gradient of the objective of the quadratic problem
        [d(Emt),~,qp_info] = quadprog(Hq,gq,Ms,y(Emt)-M(Emt,A)*x(A)-M(Emt,I)*v,[],[],-x(Emt),[],[],qp_options);
        d(A) = -x(A);
      else
        v = U\(L\(P*y(I)));			% MII*v = y(I)
        gq = MIImMIEmt'*v;			% gradient of the objective of the quadratic problem
        [d(Emt),~,qp_info] = quadprog(Hq,gq,Ms,y(Emt)-M(Emt,I)*v,[],[],-x(Emt),[],[],qp_options);
      end
      d(I) = -v-MIImMIEmt*d(Emt);
    else
      if nA
        [d(Emt),~,qp_info] = quadprog(eye(nEmt),zeros(nEmt,1),-M(Emt,Emt),y(Emt)-M(Emt,A)*x(A),[],[],-x(Emt),[],[],qp_options);
        d(A) = -x(A);
      else
        [d,~,qp_info] = quadprog(eye(n),zeros(n,1),-M,y,[],[],-x,[],[],qp_options);
      end
    end

    if qp_info ~= 1
      info.flag = values.fail_with_quadprog;
      if verb > 1; fprintf(fout,': failure ('); end
      switch qp_info
        case 0
          if verb > 1; fprintf (fout,'too many iterations)'); end
          info.flag = values.fail_with_quadprog_iter;
        case -2
          if verb > 1; fprintf (fout,'infeasible problem)'); fprintf (fout,', |d| = %12.5e',norm(d)); end
          info.flag = values.fail_with_quadprog_infeas;
        case -3
          if verb > 1; fprintf (fout,'unbounded problem)'); end
        otherwise
          if verb > 1; fprintf (fout,'exit flag = %0i)',qp_info); end
      end

      % Rerun quadprog with more output for clarification

      if verb > 1
        fprintf('\n. the QP is run again with more output to clarify what happened\n');
        qp_options = optimoptions ('Quadprog', ...
                                   'Display'             , 'iter', ...
                                   'Diagnostics'         , 'off', ...
                                   'MaxIter'             , options.qpmaxiter, ...
                                   'OptimalityTolerance' , 1.e-7, ...
                                   'Algorithm'           , 'interior-point-convex');
        if nI
          [L,U,P]   = lu(M(I,I));			% P*MII = L*U
          try
            MIImMIEmt = U\(L\(P*M(I,Emt)));		% MIImMIEmt = MII\M(I,Emt);
          catch
            info.flag = values.fail_on_degeneracy;
            return
          end
          Ms        = M(Emt,I)*MIImMIEmt - M(Emt,Emt);% negative Schur complement or Jacobian of the inequality constraints
          Hq        = eye(nEmt)+MIImMIEmt'*MIImMIEmt;	% Hessian of the objective of the quadratic problem
          if nA
            v = U\(L\(P*(y(I)-M(I,A)*x(A))));		% v = MII\(y(I)-MIA*x(A))
            gq = MIImMIEmt'*v;			% gradient of the objective of the quadratic problem
            [d(Emt),~,qp_info] = quadprog(Hq,gq,Ms,y(Emt)-M(Emt,A)*x(A)-M(Emt,I)*v,[],[],-x(Emt),[],[],qp_options);
            d(A) = -x(A);
          else
            v = U\(L\(P*y(I)));			% v = MII\y(I)
            gq = MIImMIEmt'*v;			% gradient of the objective
            [d(Emt),~,qp_info] = quadprog(Hq,gq,Ms,y(Emt)-M(Emt,I)*v,[],[],-x(Emt),[],[],qp_options);
          end
          d(I) = -v-MIImMIEmt*d(Emt);
        else
          if nA
            [d(Emt),~,qp_info] = quadprog(eye(nEmt),zeros(nEmt,1),-M(Emt,Emt),y(Emt)-M(Emt,A)*x(A),[],[],-x(Emt),[],[],qp_options);
            d(A) = -x(A);
          else
            [d,~,qp_info] = quadprog(eye(n),zeros(n,1),-M,y,[],[],-x,[],[],qp_options);
          end
        end
        fprintf('qp_info = %i, |d| = %12.5e',qp_info,norm(d));
      end

      % Nevertheless, accept the d computed by quadprog (for some erroe flags) if it is an inexact PNM direction

      if options.inexact & ismember(qp_info,[0,-2])

        IS = find( (0 <= x) & (x <= y) );
        if ~isempty(IS) & any((1-etai)*x(IS)+d(IS) > dymin)
          if verb > 1
            fprintf('\n\n. d is not accepted as an inexact PNM direction');
            I = IS(find((1-etai)*x(IS)+d(IS) > dymin));
            i = I(1);
            fprintf('\n  i = %i, x(i) = %12.5e is >= 0, y(i) = %12.5e is >= x(i) but',i,x(i),y(i));
            fprintf('\n  (1-eta_inexact)*x(i)+d(i) = (1-%g)*%12.5e+%12.5e = %12.5e should be <= dymin',etai,x(i),d(i),(1-etai)*x(i)+d(i));
          end
          d = [];
          return
        end

        IS = find( (x < 0) & (x - tau < y) );
        if ~isempty(IS) & any((1-etai)*x(IS)+d(IS) < -dymin)
          if verb > 1
            fprintf('\n\n. d is not accepted as an inexact PNM direction');
            I = IS(find((1-etai)*x(IS)+d(IS) < -dymin));
            i = I(1);
            fprintf('\n  i = %i, x(i) = %12.5e is < 0, y(i) = %12.5e is > x(i)-tau = %12.5 but',i,x(i),y(i),x(i)-tau);
            fprintf('\n  (1-eta_inexact)*x(i)+d(i) = (1-%g)*%12.5e+%12.5e = %12.5e should be >= -dymin',etai,x(i),d(i),(1-etai)*x(i)+d(i));
          end
          d = [];
          return
        end

        IS = find( (0 <= y) & (y < x) );
        if ~isempty(IS) & any((1-etai)*y(IS)+M(IS,:)*d > dymin)
          if verb > 1
            fprintf('\n\n. d is not accepted as an inexact PNM direction');
            I = IS(find((1-etai)*y(IS)+M(IS,:)*d > dymin));
            i = I(1);
            fprintf('\n  i = %i, x(i) = %12.5e is > y(i), y(i) = %12.5e is >= 0 but',i,x(i),y(i));
            fprintf('\n  (1-eta_inexact)*y(i)+M(i,:)*d = (1-%g)*%12.5e+%12.5e = %12.5e should be <= dymin',etai,y(i),M(i,:)*d,(1-etai)*y(i)+M(i,:)*d);
          end
          d = [];
          return
        end

        IS = find( (y < 0) & (y - tau < x) );
        if ~isempty(IS) & any((1-etai)*y(IS)+M(IS,:)*d < -dymin)
          if verb > 1
            fprintf('\n\n. d is not accepted as an inexact PNM direction');
            I = IS(find((1-etai)*y(IS)+M(IS,:)*d < -dymin));
            i = I(1);
            fprintf('\n  i = %i, x(i) = %12.5e is > y(i)-tau = %12.5e, y(i) = %12.5e is < 0 but',i,x(i),y(i)-tau,y(i));
            fprintf('\n  (1-eta_inexact)*y(i)+M(i,:)*d = (1-%g)*%12.5e+%12.5e = %12.5e should be >= -dymin',etai,y(i),M(i,:)*d,(1-etai)*y(i)+M(i,:)*d);
          end
          d = [];
          return
        end

        info.qpsave = true;

        if verb > 1; fprintf('\n\n. d is accepted as an inexact PNM direction'); end

        return

      else

        d = [];
        return

      end

    end

  else		% solve a single linear system

    if nI
      if verb > 1; fprintf (fout,'\n. a single linear system to solve'); end
      if nA
        d = zeros(n,1);
        d(A) = -x(A);				% d(A) = -x(A)
        [L,U,P]=lu(M(I,I));			% then P*MII = L*U
        try
          d(I) = U\(L\(P*(M(I,A)*x(A)-y(I))));	% MII*d(I) = MIA*x(A)-y(I)
        catch
          info.flag = values.fail_on_degeneracy;
          return
        end
      else
        [L,U,P]=lu(M);		% then P*M = L*U
        try
          d = -(U\(L\(P*y)));	% M*d = -y
        catch
          info.flag = values.fail_on_degeneracy;
          return
        end
      end
    else
      if verb > 1; fprintf (fout,' no linear system to solve (d = -x)'); end
      d = -x;
    end

  end

end
