function [] = certify_table(table_flag, hnm4lcp_flag, pathlcp_flag, lcpsolve_flag)

%%
% [] = certify_table(table_flag, hnm4lcp_flag, pathlcp_flag, lcpsolve_flag)
%
% This function prints on the standard output the table 4.table_flag of
% the DFG paper (reference in the README file). The run may take much
% time (several hours). To run the function, just enter
%
%    certify_table(table_flag, hnm4lcp_flag, pathlcp_flag, lcpsolve_flag)
%
% in the Matlab window, where
%
% - 'table_flag' denotes the number of the desired table in section 4 of
%   the paper,
% - 'hnm4lcp_flag' must be set to 'true' if the results on the Hnm4lcp
%   solver are desired,
% - 'pathlcp_flag' must be set to 'true' if the results on the Pathlcp
%   solver are desired,
% - 'lcpsolve_flag' must be set to 'true' if the results on the LCPsolve
%   solver are desired.
%
% To get the results of the last two solvers, these must have been
% installed before with the names 'pathlcp' and 'lcpsolve',
% respectively. Pathlcp and LCPsolve can be downloaded at the addresses
%
%    https://pages.cs.wisc.edu/~ferris/path.html
%    http://github.com/erleben/num4lcp
%
% The solvers are run twice on the first problem in the list, so that
% the computing time evaluated by tic-toc is only correct at the second
% run.

%-----------------------------------------------------------------------
% Add paths
%-----------------------------------------------------------------------

  addpath ../src			% for hnm4lcp
  addpath ~/codes/ferris/path		% PATH solver (Ferris et al.)
  addpath ~/codes/lcpsolve		% for LCPsolve
  addpath erleben			% for Erleben's test problems
  addpath diphasic			% for diphasic test problems

%-----------------------------------------------------------------------
% Check input arguments
%-----------------------------------------------------------------------

  if nargin < 4, lcpsolve_flag = false;
    if nargin < 3, pathlcp_flag = false;
      if nargin < 2, hnm4lcp_flag = false;
        if nargin < 1
          fprintf('\n### (certify_table): the first argument is required\n\n');
          return
        end
      end
    end
  end

%-----------------------------------------------------------------------
% Lists of solvers
%-----------------------------------------------------------------------

  solver_list = [];
  if hnm4lcp_flag
    if exist('hnm4lcp') ~= 2
      fprintf('\n### (certify_table): the results of ''hnm4lcp'' are desired but the solver has not been found, solver ignored\n\n');
    else
      solver_list = [solver_list, "hnm4lcp"];
    end
  end
  if pathlcp_flag
    if exist('pathlcp') ~= 2
      fprintf('\n### (certify_table): the results of ''pathlcp'' are desired but the solver has not been found, solver ignored\n\n');
    else
      solver_list = [solver_list, "pathlcp"];
    end
  end
  if lcpsolve_flag
    if exist('lcpsolve') ~= 2
      fprintf('\n### (certify_table): the results of ''lcpsolve'' are desired but the solver has not been found, solver ignored\n\n');
    else
      solver_list = [solver_list, "lcpsolve"];
    end
  end
  if isempty(solver_list)
    fprintf('\n### (certify_table): no solver required\n');
  end

  % default setting of hnm4lcp

  dxmin = 1.e-12;
  dymin = 1.e-8;

%-----------------------------------------------------------------------
% Lists of problems
%-----------------------------------------------------------------------

  if table_flag == 1	% Murty (1978)

    problem_list = ["murty-512", ...
                    "murty-512", ...
                    "murty-1024", ...
                    "murty-2048", ...
                    "murty-4096", ...
                    "murty-8192"];


  elseif table_flag == 2	% Fathy (1979)

    problem_list = ["fathi-512-1", ...
                    "fathi-512-1", ...
                    "fathi-1024-1", ...
                    "fathi-2048-1", ...
                    "fathi-4096-1", ...
                    "fathi-8192-1"];

  elseif table_flag == 3	% Fathy (1979)

    dymin = 1.e-12;
    problem_list = ["fathi-512-1", ...
                    "fathi-512-1", ...
                    "fathi-1024-1", ...
                    "fathi-2048-1", ...
                    "fathi-4096-1", ...
                    "fathi-8192-1"];

  elseif table_flag == 4	% Ben Gharbia, Gilbert (2012)

    problem_list = ["bg2012-8192", ...
                    "bg2012-8192"];


  elseif table_flag == 5	% Csizmadia-a

    problem_list = ["csizmadia-v0-n8192", ...
                    "csizmadia-v0-n8192"];

  elseif table_flag == 6	% Csizmadia-b

    problem_list = ["csizmadia-v1-n128", ...
                    "csizmadia-v1-n128", ...
                    "csizmadia-v1-n256", ...
                    "csizmadia-v1-n512", ...
                    "csizmadia-v1-n1024"];

  elseif table_flag == 7	% Lcprand

    problem_list = ["lcprand-c=n0512=na130=ni130=tol1.e-10", ...
                    "lcprand-c=n0512=na130=ni130=tol1.e-10", ...
                    "lcprand-d=n1024=na250=ni250=tol1.e-10", ...
                    "lcprand-e=n2048=na400=ni400=tol1.e-09", ...
                    "lcprand-f=n4096=na700=ni700=tol1.e-09", ...
                    "lcprand-g=n8192=na1000=ni1000=tol1.e-08"];


  elseif table_flag == 8	% Contact problem from https://github.com/erleben/num4lcp/blob/master/matlab/make_contact_matrix.m

    problem_list = ["erleben-contact=k86=s5", ...
                    "erleben-contact=k86=s5", ...
                    "erleben-contact=k171=s5", ...
                    "erleben-contact=k342=s5", ...
                    "erleben-contact=k683=s5", ...
                    "erleben-contact=k1366=s5"];

  elseif table_flag == 9	% Fluid problem from https://github.com/erleben/num4lcp/blob/master/matlab/make_contact_matrix.m

    problem_list = ["erleben-fluid=g8=s5", ...
                    "erleben-fluid=g8=s5", ...
                    "erleben-fluid=g11=s5", ...
                    "erleben-fluid=g13=s5", ...
                    "erleben-fluid=g16=s5", ...
                    "erleben-fluid=g21=s5", ...
                    "erleben-fluid=g26=s5", ...
                    "erleben-fluid=g32=s5", ...
                    "erleben-fluid=g41=s5", ...
                    "erleben-fluid=g51=s5", ...
                    "erleben-fluid=g64=s5", ...
                    "erleben-fluid=g81=s5"];

  elseif table_flag == 10	% diphasic problem

    problem_list = ["diphasic-101-a", ...
                    "diphasic-101-a", ...
                    "diphasic-101-b", ...
                    "diphasic-101-c", ...
                    "diphasic-101-d", ...
                    "diphasic-201-a", ...
                    "diphasic-201-b", ...
                    "diphasic-501-a", ...
                    "diphasic-501-b"];

  else

    fprintf('\n### (certify_table): there is no table 4.%i in the DFG paper\n\n',table_flag);
    return

  end

%-----------------------------------------------------------------------
% Solver options
%-------------------------------------------------------------------------------------------------

  hnm4lcp_options.itermax      = Inf;
  hnm4lcp_options.direction    = 0;		% 0: hybrid, 1: secured PNM
  hnm4lcp_options.tau          = 1.e-7;		% >=0: parameter tau measuring the kink proximity; -1: automatic determination of tau
  hnm4lcp_options.qpmaxdim     = 10;		% desirable bound on the QP sizes
  hnm4lcp_options.linesearch   = 1;		% 0: unit stepsize, 1: stesize halfing; 2: quadratic interpolation
  hnm4lcp_options.interpolation_safeguard = 1.e-1;	% interpolation safeguard
  hnm4lcp_options.linesearch_itermax = 1100;	% max linesearch stepsize trial per iteration (may be high for ill-conditioned problems)
  hnm4lcp_options.linesearch_m = 10;		% 1: monotone LS, >1: nonmonotone LS on the last linesearch_m values of the merit function
  hnm4lcp_options.dxmin        = dxmin;		% resolution in x
  hnm4lcp_options.dymin        = dymin;		% resolution in y

  hnm4lcp_options.scaling      = 0;	% (incide pnm4lpc) 0: no scaling, 1: scalar scaling; 2: diagonal scaling; 3: biscaling
  hnm4lcp_options.fout         = 1;
  hnm4lcp_options.verb         = 0;	% 0: silent, 1: one line per iteration, 2: standard

% Select the output channels and their level of verbosity

  fout  = 1;	% fopen('res','w');
  verb  = 0;	% channel 1 (0: silent, 1: error messages, 2: standard)

  fout2 = 1;	% fopen('res','w');
  verb2 = 0;	% channel 2 (0: silent, 1: sign vector, 2: + directions, 3: + intermediate, 4: + verification)

%-------------------------------------------------------------------------------------------------
% Run each solver in the given list on each test problem in the given list
%-------------------------------------------------------------------------------------------------

  np = length(problem_list);	% number of problems
  ns = length(solver_list);	% number of solvers

  % dash line

  dline = strcat('--------------------------');
  for js = 1:ns
    solver = solver_list(js);
    if strcmp(solver,"hnm4lcp")
      dline = strcat(dline,'-------------------------------------------------');
    elseif strcmp(solver,"pathlcp")
      dline = strcat(dline,'-----------');
    elseif strcmp(solver,"lcpsolve")
      dline = strcat(dline,'-----------');
    end
  end

  % head of the table

  fprintf('\nTable 4.%i of the DFG paper',table_flag);
  fprintf('\n%s',dline);
  fprintf('\n|        Problem         |');
  for js = 1:ns
    solver  = solver_list(js);
    slength = strlength(solver);
    if strcmp(solver,"hnm4lcp")
      fprintf('%s%s%s|',repmat(' ',1,ceil((48-slength)/2)),solver,repmat(' ',1,floor((48-slength)/2)));
    elseif strcmp(solver,"pathlcp")
      fprintf('%s%s%s|',repmat(' ',1,ceil((10-slength)/2)),solver,repmat(' ',1,floor((10-slength)/2)));
    elseif strcmp(solver,"lcpsolve")
      fprintf('%s%s%s|',repmat(' ',1,ceil((10-slength)/2)),solver,repmat(' ',1,floor((10-slength)/2)));
    end
  end

  fprintf('\n|     Name            n  |');
  for js = 1:ns
    solver =  solver_list(js);
    if strcmp(solver,"hnm4lcp")
      fprintf('   iter  qp   |qp|   %%ls   alpha   ac1s    sec  |');
    elseif strcmp(solver,"pathlcp")
      fprintf('    sec   |');
    elseif strcmp(solver,"lcpsolve")
      fprintf('    sec   |');
    end
  end
  fprintf('\n%s',dline);

%-----------------------------------------------------------------------
% Big loop on the problems and solvers
%-----------------------------------------------------------------------

  for ip = 1:np

    %-------------------------------------------------------------------
    % Select one problem and get its data
    %-------------------------------------------------------------------

    problem = problem_list{ip};		% with {}, problem is a string vector, not a string arrwy

    abstol = 1.e-10;	% default precision

    if strcmp(problem(1:5),'murty')

      problemname = 'murty';

      dash = regexp(problem,'-');
      n    = str2num(problem(dash(1)+1:end));
      fprintf('\n| murty%s  %6i |',repmat(' ',14-strlength(problemname),1),n);

      [M,q,xsol,x1] = murty(n,0);

    elseif strcmp(problem(1:5),'fathi')

      problemname = 'fathi';

      dash = regexp(problem,'-');
      n = str2num(problem(dash(1)+1:dash(2)-1));
      m = str2num(problem(dash(2)+1:end));
      fprintf('\n| fathi%s  %6i |',repmat(' ',14-strlength(problemname),1),n);

      [M,q,xsol,x1] = fathi(n,m,0);

    elseif strcmp(problem(1:6),'bg2012')

      problemname = 'bg2012';

      dash = regexp(problem,'-');
      n    = str2num(problem(dash(1)+1:end));
      fprintf('\n| bg2012%s  %6i |',repmat(' ',14-strlength(problemname),1),n);

      [M,q,xsol,x1] = bg2012(n,0);
      M = full(M);

    elseif strcmp(problem(1:9),'csizmadia')

      problemname = problem(1:9);

      dash  = regexp(problem,'-');
      ndash = length(dash);
      for i = 1:ndash
        i1 = dash(i)+1;
        if i == ndash
          i2 = length(problem);
        else
          i2 = dash(i+1)-1;
        end
        if strcmp(problem(i1:i1),'v')
          var = str2num(problem(i1+1:i2));
          if var == 0
            problemname = strcat(problemname,'-a');
          elseif var == 1
            problemname = strcat(problemname,'-b');
            linesearch_itermax = 1100;
            abstol = 1.e-12;
          elseif var == 2
            problemname = strcat(problemname,'-c');
          else
            if verb; fprintf('\n| %s: unknown problem\n\n',problemname);
            end
          end
        elseif strcmp(problem(i1:i1),'n')
          n = str2num(problem(i1+1:i2));
        end
      end

      if (var == 1)
        if n == 128
          abstol = 1.e-15;
        elseif n == 512
          abstol = 1.e-14;
        end
      end

      fprintf('\n| %s%s  %6i |',problemname,repmat(' ',14-strlength(problemname),1),n);

      [M,q,xsol,x1] = csizmadia(n,var,0);

    elseif strcmp(problem(1:7),'lcprand')

      problemname = problem(1:7);

      equal = regexp(problem,'=');
      nequal = length(equal);
      for i = 1:nequal
        i1 = equal(i)+1;
        if i == nequal
          i2 = length(problem);
        else
          i2 = equal(i+1)-1;
        end
        if strcmp(problem(i1:i1+1),'na')
          na = str2num(problem(i1+2:i2));
        elseif strcmp(problem(i1:i1+1),'ni')
          ni = str2num(problem(i1+2:i2));
        elseif strcmp(problem(i1:i1),'n')
          n = str2num(problem(i1+1:i2));
        elseif strcmp(problem(i1:i1+2),'tol')
          abstol = str2num(problem(i1+3:i2));
        end
      end

      fprintf('\n| %s%s  %6i |',problemname,repmat(' ',14-strlength(problemname),1),n);

      [M,q,xsol,x1] = lcprand(n,na,ni,0);

    elseif (length(problem)>14) && strcmp(problem(1:15),'erleben-contact')

      problemname = problem(9:15);

      equal = regexp(problem,'=');
      nequal = length(equal);
      for i = 1:nequal
        i1 = equal(i)+1;
        if i == nequal
          i2 = length(problem);
        else
          i2 = equal(i+1)-1;
        end
        if strcmp(problem(i1:i1),'k')
          k = str2num(problem(i1+1:i2));
          n = k*6;
        elseif strcmp(problem(i1:i1),'s')
          s = 0.1*str2num(problem(i1+1:i2));
        end
      end

      fprintf('\n| %s%s  %6i |',problemname,repmat(' ',14-strlength(problemname),1),n);

      % get M
      M        = make_contact_matrix(k);

      % force M to have a positive diagonal

      diagM    = diag(M);
      diagMmin = min(diagM(diagM>0));
      M        = M + 1.e-8*diagMmin*eye(n);

      % get q and the solution

      [xsol,q] = make_lcp_corrected(M,s);		% xsol is the solution with sparcity s in (0,1)
      x1       = [];

    elseif strcmp(problem(1:13),'erleben-fluid')

      problemname = problem(9:13);

      equal = regexp(problem,'=');
      nequal = length(equal);
      for i = 1:nequal
        i1 = equal(i)+1;
        if i == nequal
          i2 = length(problem);
        else
          i2 = equal(i+1)-1;
        end
        if strcmp(problem(i1:i1),'g')
          g = str2num(problem(i1+1:i2));
          n = g*g*g;
        elseif strcmp(problem(i1:i1),'s')
          s = 0.1*str2num(problem(i1+1:i2));
        end
      end

      fprintf('\n| %s%s  %6i |',problemname,repmat(' ',14-strlength(problemname),1),n);

      % get M

%     M = make_fluid_matrix(g);		% dense version
      M = make_fluid_matrix_sp(g);	% sparse version

      % get q and the solution
      [xsol,q] = make_lcp_corrected(M,s);		% xsol is the solution with sparcity s in (0,1)
      x1       = [];

    elseif strcmp(problem(1:8),'diphasic')

      problemname = problem;

      % get M and q

      [M,q] = diphasic(strcat(problem,'.tdt'),0,0);
      xsol  = [];
      x1    = [];
      n     = length(q);

      fprintf('\n| %s%s  %6i |',problemname,repmat(' ',14-strlength(problemname),1),n);

    else

      fprintf('\n| %s%s       | unkwon problem',problemname,repmat(' ',13-strlength(problem),1));
      continue

    end

    %-------------------------------------------------------------------
    % Run the solvers
    %-------------------------------------------------------------------

    for js = 1:ns

      solver = solver_list(js);

      %-----------------------------------------------------------------
      if strcmp(solver,"hnm4lcp")
      %-----------------------------------------------------------------

        hnm4lcp_options.abstol = abstol;

        tic
        [x,info] = hnm4lcp(M,q,x1,hnm4lcp_options);

        if ~isempty(xsol)
          if (info.flag == 0) && ~isempty(x)
            if (norm(x-xsol) <= 1.e-8) || (norm(x-xsol)/norm(xsol) <= 1.e-8)
              fprintf('%7i %4i',info.iter,info.qpsolve);
              if info.qpsolve > 0
                fprintf(' %5.1f ',info.qpnvar/info.qpsolve);
              else
                fprintf('   --  ');
              end
              if info.linesearches == 0
                fprintf(' %5.2f    --   ',0);
              else
                fprintf(' %5.2f  %7.1e',info.linesearches*100/info.iter,info.stepsize_log2median);
              end
              fprintf(' %3i %8.2f |',info.nb_unit_stepsizes,toc);
            else
              fprintf('                      (13)                      |');
            end
          else
            fprintf('                      (%i)                       |',info.flag);
          end
        else
          fprintf('?%6i %4i',info.iter,info.qpsolve);
          if info.qpsolve > 0
            fprintf(' %5.1f ',info.qpnvar/info.qpsolve);
          else
            fprintf('   --  ');
          end
          if info.linesearches == 0
            fprintf(' %5.2f    --   ',0);
          else
            fprintf(' %5.2f  %7.1e',info.linesearches*100/info.iter,info.stepsize_log2median);
          end
          fprintf(' %3i %8.2f |',info.nb_unit_stepsizes,toc);
        end

      %-----------------------------------------------------------------
      elseif strcmp(solver,"pathlcp")
      %-----------------------------------------------------------------

        try
          tic
          [x,mu] = pathlcp(M,q);
          if isempty(x)
            fprintf('   (3)    |');
          elseif ~isempty(xsol)
            if (norm(x-xsol) <= 100*abstol) || (norm(x-xsol)/norm(xsol) <= 100*abstol)
              fprintf(' %8.2f |',toc);
            else
              fprintf('   (1)    |');
            end
          else
              fprintf('?%8.2f |',toc);
          end
        catch ME
          ME_message = ME.message;
          if strcmp(ME_message(1:19),'Solution not finite')
            fprintf('   (2)    |');
          elseif strcmp(ME_message(1:27),'Path fails to solve problem')
            fprintf('   (4)    |');
          else
            fprintf('   (?)    |');
          end
        end

      %-----------------------------------------------------------------
      elseif strcmp(solver,"lcpsolve")
      %-----------------------------------------------------------------

        tic
        [y,x,retcode] = LCPsolve(M,q);
        if ~isempty(xsol)
          if ~isempty(x) && ( (norm(x-xsol) <= abstol) || (norm(x-xsol)/norm(xsol) <= abstol) )
            fprintf(' %8.2f |',toc);
          else
            fprintf('    (%i)   |',retcode(1));
          end
        else
          fprintf('?%8.2f |',toc);
        end

      %-----------------------------------------------------------------

      end

    end

  end

  fprintf('\n%s\n\n',dline);

  return
