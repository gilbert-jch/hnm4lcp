% Matlab setting

  clc		% clear command window
  clf
  close all	% close all figures
  clear all
  format long
  format compact

% Modify the Matlab search path

  addpath ../src/nmvar					% for nmvar
  addpath ../src/hnm4lcp				% for hnm4lcp
  addpath ../data					% for data
% addpath ~/codes/compecon2011_64_20110718/CEtools	% for lcplemke from CompEcon (Fackler and Miranda)
  addpath ~/codes/ferris/path				% PATH solver (Ferris et al.), version of Oct. 2024 for 'arm64' architecture
  addpath ~/codes/lcpsolve				% for LCPsolve
  addpath erleben					% for Erleben's test problems
  addpath diphasic					% for diphasic test problems

% Set the solver to run data

  solver.hnm4lcp  = true;
  solver.pathlcp  = false;
  solver.lcpsolve = false;

  solver.nmvar    = false;

%---------------------------------------
% Select data
%---------------------------------------

  rng default

  abstol = 1.e-10;	% default precision
  dxmin  = 1.e-12;	% default resolution in x
  dymin  = 1.e-8;	% default resolution in y

  linesearch_itermax = 200;

%:::::::::::::::::::::::::::::::::::::::
% Murty problem
%:::::::::::::::::::::::::::::::::::::::

  scaling = 0;

% n  =   62;

% n  =   64;
% n  =  128;
% n  =  256;
% n  =  512;
% n  = 1024;
% n  = 2048;
% n  = 4096;
% n  = 8192;

% [M,q,xsol,x1] = murty(n,1);

% M = sparse(M);
% x1 = (rand(n,1)-2.0)*10;

%:::::::::::::::::::::::::::::::::::::::
% Fathi problem
%:::::::::::::::::::::::::::::::::::::::

% scaling = 0;

% n  =   44;

% n  =    4;
% n  =   64;
% n  =  128;
% n  =  256;
% n  =  512;
% n  = 1024;
% n  = 2048;
% n  = 4096;
% n  = 8192;

% np = 1;
% [M,q,xsol,x1] = fathi(n,np,1);

% x1 = [10/8;-7/8;5/8;-3/8;1/8]
% x1 = [1;-1;1;-1;0]
% M*x1+q
% x1 = [0;1;-1]
% x1 = [-1.0976153215320070266614038700936;0.37356087289743028501121102635807;0.12643912505645274468868421990919]
% x1 = [0;0;1/7]
% x1 = [-0.05;1/10;1/14]
% x1 = [0;0;8/9;0]
% rng(16,'twister');
% x1 = (rand(n,1)-0.5)*100;
% M = sparse(M);

%:::::::::::::::::::::::::::::::::::::::
% BenGharbia-Gilbert (2012) problems (cycling with unit stepsize, solved in 2 iterations with Armijo LS)
%:::::::::::::::::::::::::::::::::::::::

% n = 8192;
% [M,q,xsol,x1] = bg2012(n,0);
% M = full(M);

%:::::::::::::::::::::::::::::::::::::::
% Csizmadia problem
%:::::::::::::::::::::::::::::::::::::::

  % No scaling
% n =   32; scaling = 0;
% n =   64; scaling = 0;
% n =  128; scaling = 0;                            abstol = 1.e-15;
% n =  256; scaling = 0; linesearch_itermax =  300; abstol = 1.e-12;
% n =  512; scaling = 0; linesearch_itermax =  500; abstol = 1.e-14;
% n = 1024; scaling = 0; linesearch_itermax = 1100; abstol = 1.e-14; %fails
% n = 2048; scaling = 0;
% n = 4096; scaling = 0;
% n = 8192; scaling = 0;

  % Scaling 1
% n =   32; scaling = 1;	% OK
% n =   64; scaling = 1;	% fails

  % Scaling 2
% n =   32; scaling = 2;	% OK
% n =   64; scaling = 2;	% OK
% n =  128; scaling = 2;	% fails

  % Scaling 3
% n =    4; scaling = 3; scaling_prec = 1.e-3; scaling_norm = 2; abstol = 1.e-05;
% n =   64; scaling = 3; scaling_prec = 1.e-3; scaling_norm = 2; abstol = 1.e-26;
% n =  128; scaling = 3; scaling_prec = 1.e-3; scaling_norm = 2; abstol = 1.e-10;		% OK
% n =  256; scaling = 3; scaling_prec = 1.e-3; scaling_norm = 2; abstol = 1.e-10;		% OK
% n =  512; scaling = 3; scaling_prec = 1.e-2; scaling_norm = 2; abstol = 1.e-10;		% fails (positive slope)
%---
% n =  512; scaling = 3; scaling_prec = 1.e-3; scaling_norm = 2; abstol = 1.e-10;		% fails (cond(M) = 6, norm(q) = 3.e+119 and the solver cycles)
% n = 1024; scaling = 3; scaling_prec = 1.e-2; scaling_norm = 2; abstol = 1.e-10;		% theta=0 but not the right solution
% n = 2048; scaling = 3; scaling_prec = 1.e-2; scaling_norm = 2; abstol = 1.e-10;		% fails at the first iteration, cond of the scaling 8.e+117...

  % Starting at zero
% n =  128; scaling = 0; abstol = 1.e-15;	% OK
% n =  256; scaling = 0; abstol = 1.e-15;	% ko
% n =  512; scaling = 0;	% ko

% [M,q,xsol,x1] = csizmadia(n,1,1);

% x1 = zeros(n,1);

%:::::::::::::::::::::::::::::::::::::::
% Hungerländer-Rendl (2015) sparse problems
%:::::::::::::::::::::::::::::::::::::::

% n = 2000;
% small = 1%.e-14;	% should be set to 1, 1.e-5, 1.e-10, 1.e-14
% [M,q,xsol,x1] = hr2015(n,small);
% M = full(M);
% abstol = 1.e-9;

%:::::::::::::::::::::::::::::::::::::::
% Sparse random problems with positive definite matrix
% Curtis, Han, and Robinson (2015) generation style
%:::::::::::::::::::::::::::::::::::::::

% n =  4096; dens =   0.1; cn = 1.e02; abstol = 1.e-10;
% n =  4096; dens =   0.1; cn = 1.e06; abstol = 1.e-08;
% n =  4096; dens =  0.01; cn = 1.e02; abstol = 1.e-10;
% n =  4096; dens =  0.01; cn = 1.e06; abstol = 1.e-09;
% n =  4096; dens = 0.001; cn = 1.e02; abstol = 1.e-10;
% n =  4096; dens = 0.001; cn = 1.e06; abstol = 1.e-09;

% n =  8192; dens =   0.1; cn = 1.e02; abstol = 1.e-10;
% n =  8192; dens =   0.1; cn = 1.e06; abstol = 1.e-08;
% n =  8192; dens =  0.01; cn = 1.e02; abstol = 1.e-10;
% n =  8192; dens =  0.01; cn = 1.e06; abstol = 1.e-09;
% n =  8192; dens = 0.001; cn = 1.e02; abstol = 1.e-10;
% n =  8192; dens = 0.001; cn = 1.e06; abstol = 1.e-09;


% n =  5000; dens =   0.1; cn = 1.e02; abstol = 1.e-10;
% n =  5000; dens =   0.1; cn = 1.e06; abstol = 1.e-07;
% n =  5000; dens =   0.1; cn = 1.e10; abstol = 1.e+02;	% 9.e-9 pour S2
% n =  5000; dens =  0.01; cn = 1.e02; abstol = 1.e-12;
% n =  5000; dens =  0.01; cn = 1.e06; abstol = 1.e-08;
% n =  5000; dens =  0.01; cn = 1.e10; abstol = 1.e+01;	% 1.e-9 pour S2 ou alors test d'arret si pas unité et pas de changement d'ensemble actif
% n =  5000; dens = 0.001; cn = 1.e02; abstol = 1.e-10;
% n =  5000; dens = 0.001; cn = 1.e06; abstol = 1.e-04;	% 1.e-9 pour S2
% n =  5000; dens = 0.001; cn = 1.e10; abstol = 1.e+01;	% 1.e-6 pour S2
% n = 10000; dens = 0.001; cn = 1.e02; abstol = 1.e-10;
% n = 10000; dens = 0.001; cn = 1.e06; abstol = 1.e-07;
% n = 10000; dens = 0.001; cn = 1.e10; abstol = 1.e-05;	% 1.e+1;	% 1.e-6 pour S2

% [M,q,xsol,x1] = chr2015(n,dens,cn,1);
% if isempty(M); return; end
% M    = full(M);
% q    = full(q);
% xsol = full(xsol);
% x1   = full(x1);

%:::::::::::::::::::::::::::::::::::::::
% Curtis-Han-Robinson counter-example
%:::::::::::::::::::::::::::::::::::::::

% [M,q,xsol,x1] = curtis_han_robinson_2015;

%::::::::
% Lcprand
%::::::::

% [M,q,xsol,x1] = lcprand( 128, 40, 40,1); dymin = 1.e-12;
% [M,q,xsol,x1] = lcprand( 256, 70, 70,1); dymin = 1.e-11;
% [M,q,xsol,x1] = lcprand( 512,130,130,1); n = 512; dymin = 1.e-11;
% [M,q,xsol,x1] = lcprand(1024,250,250,1); dymin = 1.e-11;
% [M,q,xsol,x1] = lcprand(2048,400,400,1); abstol = 1.e-09; dymin = 1.e-10;
% [M,q,xsol,x1] = lcprand(4096,700,700,1); abstol = 1.e-09; dymin = 1.e-10;
% [M,q,xsol,x1] = lcprand(8192,1000,1000,1); abstol = 1.e-08; dymin = 1.e-09;

%:::::::::::::::::::::::::::::::::::::::

%M = M + 0.01*diag(ones(n-1,1),1);

% n = 16; na = 5; ni = 5;
% n = 2; na = 1; ni = 1;
% n = 4; na = 2; ni = 2;
% n = 5; na = 2; ni = 2;
% n = 256; na = 64; ni = 128;
% n = 512; na = 128; ni = 256;
% n = 1024; na = 256; ni = 512; abstol = 1.e-9;
% n = 2048; na = 512; ni = 1024; abstol = 1.e-8;
% n = 4096; na = 1024; ni = 2048; abstol = 1.e-8;
% n = 8192; na = 2048; ni = 4096; abstol = 1.e-7;
% n = 16384; na = 4096; ni = 8192; abstol = 1.e-7;
% [M,q,xsol,x1] = lcprand(n,na,ni);

% n = 2048; na =  512; ni = 1024; density = 1.e-0;
% n = 2048; na =  512; ni = 1024; density = 1.e-1;
% n = 2048; na =  512; ni = 1024; density = 1.e-2;
% n = 2048; na =  512; ni = 1024; density = 1.e-3;
% n = 2048; na =  512; ni = 1024; density = 1.e-4;
% n = 2048; na =  512; ni = 1024; density = 1.e-5;
% n = 4096; na = 1024; ni = 2048; density = 1.e-0;
% n = 4096; na = 1024; ni = 2048; density = 1.e-1;
% n = 4096; na = 1024; ni = 2048; density = 1.e-2;
% n = 4096; na = 1024; ni = 2048; density = 1.e-3;
% n = 4096; na = 1024; ni = 2048; density = 1.e-4;
% n = 4096; na = 1024; ni = 2048; density = 1.e-5;
% n = 8192; na = 2048; ni = 4096; density = 1.e-2;
% n = 8192; na = 2048; ni = 4096; density = 1.e-3;
% n = 8192; na = 2048; ni = 4096; density = 1.e-4;
% n = 8192; na = 2048; ni = 4096; density = 1.e-5;
% [M,q,xsol,x1] = lcprandsp(n,na,ni,density);

% [M,q,xsol,x1] = lcphard;		% Author Jean-Pierre Dussault
% x1 = [];

%:::::::::::::::::::::::::::::::::::::::
% Contact problem from https://github.com/erleben/num4lcp/blob/master/matlab/make_contact_matrix.m
%:::::::::::::::::::::::::::::::::::::::

% k =    1; s = 0.5; n = 6*k;	% n =    6
% k =    2; s = 0.5; n = 6*k;	% n =   12
% k =    3; s = 0.5; n = 6*k;	% n =   18
% k =    4; s = 0.5; n = 6*k;	% n =   24
  k =    5; s = 0.5; n = 6*k;	% n =   30
% k =   86; s = 0.5; n = 6*k;	% n =  516
% k =  171; s = 0.5; n = 6*k;	% n = 1026
% k =  342; s = 0.5; n = 6*k;	% n = 2052
% k =  683; s = 0.5; n = 6*k;	% n = 4098
% k = 1366; s = 0.5; n = 6*k;	% n = 8196

  % get M (dense matrix)

  M        = make_contact_matrix(k);		% full matrix
% M        = make_contact_matrix_sp(k);		% sparse matrix

  % force M to have a positive diagonal

  diagM    = diag(M);
  diagMmin = min(diagM(diagM>0));
  M        = M + 1.e-8*diagMmin*eye(n);

  [xsol,q] = make_lcp_corrected(M,s);		% xsol is the solution with sparcity s in (0,1)
  x1       = [];

% if ~ispmat(M,1)
%   fprintf('\n\n### the matrix M is not a P-matrix\n\n');
%   return
% else
%   fprintf('\n\n### the matrix M is a P-matrix\n\n');
% end

  % introductory message

  cl = clock;
  dstr = sprintf('%i-%i-%i, %i:%i:%i',cl(3),cl(2),cl(1),cl(4),cl(5),fix(cl(6)));
  fprintf('\nErleben contact problem');
  fprintf('\n. k = %i',k);
  fprintf('\n. n = %i\n',n);

%:::::::::::::::::::::::::::::::::::::::
% Fluid problem from https://github.com/erleben/num4lcp/blob/master/matlab/make_contact_matrix.m
%:::::::::::::::::::::::::::::::::::::::

% g = grid size, n = g^3

% g =    2; s = 0.5; n = g*g*g;		% n =      8
% g =    3; s = 0.5; n = g*g*g;		% n =     27
% g =    4; s = 0.5; n = g*g*g;		% n =     64
% g =    5; s = 0.5; n = g*g*g;		% n =    125
% g =    6; s = 0.5; n = g*g*g;		% n =    216
% g =    7; s = 0.5; n = g*g*g;		% n =    343
% g =    9; s = 0.5; n = g*g*g;		% n =    729

% g =    8; s = 0.5; n = g*g*g;		% n =    512
% g =   11; s = 0.5; n = g*g*g;		% n =   1331
% g =   13; s = 0.5; n = g*g*g;		% n =   2197
% g =   16; s = 0.5; n = g*g*g;		% n =   4096
% g =   21; s = 0.5; n = g*g*g;		% n =   9261
% g =   26; s = 0.5; n = g*g*g;		% n =  17576
% g =   32; s = 0.5; n = g*g*g;		% n =  32768
% g =   41; s = 0.5; n = g*g*g;		% n =  68921
% g =   51; s = 0.5; n = g*g*g;		% n = 132651
% g =   64; s = 0.5; n = g*g*g;		% n = 262144

  % get M

%>   M        = make_fluid_matrix(g);
%> 
%> % if ~ispmat(M,1)
%> %   fprintf('\n\n### the matrix M is not a P-matrix');
%> %   return
%> % end
%> 
%>   [xsol,q] = make_lcp_corrected(M,s);		% xsol is the solution with sparcity s in (0,1)
%>   x1       = [];
%> 
%>   % introductory message
%> 
%>   cl = clock;
%>   dstr = sprintf('%i-%i-%i, %i:%i:%i',cl(3),cl(2),cl(1),cl(4),cl(5),fix(cl(6)));
%>   fprintf('\nErleben fluid problem');
%>   fprintf('\n. n = %i\n',n);

%:::::::::::::::::::::::::::::::::::::::
% Diphasic (dense matrix)
%:::::::::::::::::::::::::::::::::::::::

% [M,q,xsol,x1] = diphasic('diphasic-101-a.tdt',0,0);
% [M,q,xsol,x1] = diphasic('diphasic-101-b.tdt',0,0);
% [M,q,xsol,x1] = diphasic('diphasic-101-c.tdt',0,0);
% [M,q,xsol,x1] = diphasic('diphasic-101-d.tdt',0,0);
% [M,q,xsol,x1] = diphasic('diphasic-201-a.tdt',0,0);
% [M,q,xsol,x1] = diphasic('diphasic-201-b.tdt',0,0);
% [M,q,xsol,x1] = diphasic('diphasic-501-a.tdt',0,0);
% [M,q,xsol,x1] = diphasic('diphasic-501-b.tdt',0,0);

% n = length(q)

% ispmat(M,1);
% return
% ispmat(M)
% keyboard
% M = M+1.e-3*eye(size(M));
% eigM = eig(M);
% mineigM = min(eigM)
% maxeigM = max(eigM)

% ispmat(M)

%:::::::::::::::::::::::::::::::::::::::

% tic
% if issparse(M)
%   S = zeros(n,1);
%   MT = M';		% since the indices of the columns are sorted
%   for i = 1:n
%     s = 1/norm(MT(:,i));
%     MT(:,i) = s*MT(:,i);
%     S(i) = s;
%   end
%   M = MT';
%   q = S.*q;
% else
%   S = zeros(n,1);
%   for i = 1:n
%     s = 1/norm(M(i,:));
%     M(i,:) = s*M(i,:);
%     S(i) = s;
%   end
%   q = S.*q;
% end
% fprintf('\nScaling in %g sec\n',toc);

% tic
% [M,DL,dr] = nmvar_biscale (M,2,5.e-1);
% fprintf('\nScalin in %g sec\n',toc);
% q = DL*q;

%> tic
%> lcpsolve(-M,-q,zeros(n,1),1.e10*ones(n,1))
%> fprintf('%g sec\n',toc);
%> 
%> tic
%> lcplemke(-M,-q,zeros(n,1),1.e10*ones(n,1))
%> fprintf('%g sec\n',toc);

% reset(RandStream.getDefaultStream);	% reset the default random stream to its initial state
% q = (rand(n,1)-0.5)*100;
% x1 = (rand(n,1)-0.5)*100;
% x1 = ones(n,1);

% % example with a minimum taken along the first segment
% M = [ 2  2
%      -3  4 ];
% q = [-1; 1];
% x1 = [-0.5; 0.75];

%---------------------------------------
% Scaling
%---------------------------------------

  tic_scaling = tic;	% Launch the timer for scaling

  if scaling == 0

    dymin = dymin*ones(n,1);
   
  elseif scaling == 1
   
    sigma = 1/norm(M,'fro');
    M = sigma*M;
    q = sigma*q;
%   y = M*x+q;
    dymin = sigma*dymin*ones(n,1);

  elseif scaling == 2

    options.fout = 1;
    options.verb = 1;
  
    [M,q] = hnm4lcp_lscale (M,q,options);
%   y = M*x+q;

    dymin = dl.*dymin;

  elseif scaling == 3

    options.fout = 1;
    options.verb = 1;
    [M,q,dl,dr] = hnm4lcp_biscale (M,q,scaling_norm,scaling_prec,options);

%   [M,q,dl,dr] = hnm4lcp_biscale_ruiz (M,q,1,scaling_prec,1,1);

%   [x,info] = hnm4lcp(M,q,x1,hnm4lcp_options);

% [M,q,dl,dr] = hnm4lcp_biscale (M,q,2,1.e+0,options);
% [M,q,dl,dr] = hnm4lcp_biscale (M,q,2,5.e-1,options);
% [M,q,dl,dr] = hnm4lcp_biscale (M,q,2,1.e-1,options);
% [M,q,dl,dr] = hnm4lcp_biscale (M,q,2,1.e-2,options);
% [M,q,dl,dr] = hnm4lcp_biscale (M,q,2,1.e-3,options);
% [M,q,dl,dr] = hnm4lcp_biscale (M,q,2,1.e-4,options);
% [M,q,dl,dr] = hnm4lcp_biscale (M,q,Inf,1.e-0,options);
% [M,q,dl,dr] = hnm4lcp_biscale (M,q,Inf,1.e-1,options);
% [M,q,dl,dr] = hnm4lcp_biscale (M,q,Inf,1.e-2,options);
% [M,q,dl,dr] = hnm4lcp_biscale (M,q,Inf,1.e-4,options);

    if ~isempty(x1); x1 = x1./dr; end
    if ~isempty(xsol); xsol = xsol./dr; end

    dymin = dl.*dymin;

%dldeb=dl(1:10)
%dlfin=dl(end-9:end)
%drdeb=dr(1:10)
%drfin=dr(end-9:end)

  elseif scaling == 4

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
    fprintf ('\n. Diagonal scaling of M and q by      %8.2e (average value)',sigma);
    if sigma2 ~= 0; fprintf ('\n. Additional scaling of q by          %8.2e',1/sigma2); end
    y = M*x+q;

  elseif scaling == 5	% obsolete

    options.fout = 1;
    options.verb = 1;

    [M,q,dl,dr] = hnm4lcp_biscale2 (M,q,2,scaling_prec,options);

  end

  scaling_time = 0;
  if scaling
    fprintf('\nScaling time = %g',toc(tic_scaling));
  end

%---------------------------------------
% Solver NMVAR (don't remember waht it does...)
%---------------------------------------

  if solver.nmvar

  % Set options of the solver

    nmvar_options.itermax             = Inf;
    nmvar_options.scaling             = 0;	% 0: no scaling, 1: scalar scaling; 2: diagonal scaling; 3: biscaling
    nmvar_options.direction           = 3;	% 0: B-Newton-min, 1: Newton-min (several options, including Harker-Pang), 2: Newton-min-var, 3: Newton-min-hybrid
  % nmvar_options.direction           = 3;	% 0: Newton-min, 1: Newton-min-var, 2: B-Newton-min, 3: Newton-min-hybrid
    nmvar_options.linesearch          = 2;	% 0: no LS, 1: LS, m>1: nonmonotone LS on the last m values of the merit function
    nmvar_options.linesearch_interpol = false;
    nmvar_options.linesearch_m        = 1;	% 1: monotone LS, >1: nonmonotone LS on the last linesearch_m values of the merit function
    nmvar_options.linesearch_special  = 0;	% 0: unit stepsize; 1: exact stepsize (not yet implemented); 2: Harker-Pang stepsize; 3: modified Harker-Pang stepsize
    nmvar_options.abstol              = abstol;
    nmvar_options.verb                = 3;	% 0: silent, 1:, 2: standard

  % nmvar_options.reltol              = 1.e-8;
    nmvar_options.dxmin               = dxmin;	% resolution in x

    % be more demanding without scaling, since otherwise the solution may not be found

    if (nmvar_options.scaling == 0) & isfield(nmvar_options,'reltol')
      nmvar_options.reltol = nmvar_options.reltol/norm(M,'fro');
    end
  % if (nmvar_options.scaling == 3) & isfield(nmvar_options,'abstol')
  %   nmvar_options.abstol = 1.e-12;
  % end

  % Call the solver nmvar

    tic
    [x,info] = nmvar(M,q,x1,nmvar_options);
    if info.flag == 0
      fprintf('\n(nmvar) output flag = %i\n(nmvar) number of iterations = %i\n(nmvar) number of QP solves = %i',info.flag,info.iter,info.qpsolve);
      if ~isempty(xsol) & ~isempty(x)
        xsoln = norm(xsol);
        if xsoln ~= 0
          fprintf('\n(nmvar) L2 relative error on the solution = %11.5e',norm(x-xsol)/xsoln);
        else
          fprintf('\n(nmvar) L2 absolute error on the solution = %11.5e',norm(x-xsol));
        end
      end
      fprintf('\n(nmvar) CPU time = %g sec\n',toc);
    end
  % x = lcpsolve(-M,-q,zeros(n,1),1.e20*ones(n,1),x1);
  % x = lemke(M,q,x1)
  % x = lcp(M,q,[],[],[],true)

  end

%---------------------------------------
% Solver HNM4LCP
%---------------------------------------

  if solver.hnm4lcp

    hnm4lcp_options.itermax            = Inf;		%Inf;
    hnm4lcp_options.direction          = 0;		% 0: hybrid, 1: secured PNM
    hnm4lcp_options.tau                = 1.e-7;	% >=0: parameter tau measuring the kink proximity; -1: automatic determination of tau
    hnm4lcp_options.qpmaxdim           = 10;		% desirable bound on the QP sizes
    hnm4lcp_options.qptol              = 1.e-7;	% desirable precision on the solution of the QP problems
    hnm4lcp_options.linesearch         = 1;		% 0: unit stepsize, 1: stesize halving; 2: quadratic interpolation
    hnm4lcp_options.linesearch_itermax = linesearch_itermax;	% max linesearch stepsize trial per iteration (may be high for ill-conditioned problems)
    hnm4lcp_options.linesearch_m       = 10;		% 1: monotone LS, >1: nonmonotone LS on the last linesearch_m values of the merit function
    hnm4lcp_options.abstol             = abstol;
    hnm4lcp_options.dxmin              = dxmin;	% resolution in x
    hnm4lcp_options.dymin              = dymin;	% resolution in y

    hnm4lcp_options.scaling = 0;	% 0: no scaling, 1: scalar scaling; 2: diagonal scaling; 3: biscaling
    hnm4lcp_options.verb    = 2;	% 0: silent, 1: one line per iteration, 2: detailed

    tic
    [x,info] = hnm4lcp(M,q,x1,hnm4lcp_options);
    if ~isempty(xsol) & ~isempty(x)
      xsoln = norm(xsol);
      if xsoln ~= 0
        fprintf('\n(hnm4lcp) L2 relative error on the solution = %11.5e',norm(x-xsol)/xsoln);
      else
        fprintf('\n(hnm4lcp) L2 absolute error on the solution = %11.5e',norm(x-xsol));
      end
    end
    fprintf('\n(hnm4lcp) CPU time = %g sec\n',toc);

  end

%---------------------------------------
% Solver PATHLCP
%---------------------------------------

  if solver.pathlcp

    scaling = 0;
    options.verb = 1;
    options.fout = 1;

%     tic_scaling = tic;	% Launch the timer for scaling
% 
%     if scaling == 0
% 
%       fprintf ('\n. No scaling');
% 
%     elseif scaling == 1
%      
%       sigma = 1/norm(M,'fro');
%       M = sigma*M;
%       q = sigma*q;
%       fprintf ('\n\nScaling M and q');
%       fprintf ('\n. Scalar scaling factor = %8.2e',sigma);
%       y = M*x+q;
% 
%     elseif scaling == 2
% 
%       fprintf ('\n\nScaling M and q');
%       [M,q] = hnm4lcp_lscale (M,q,options);
% %     y = M*x+q;
% 
%     elseif scaling == 3
% 
%       fprintf ('\n\nScaling M and q');
% 
% %   [M,q,dl,dr] = hnm4lcp_biscale (M,q,2,1.e+0,options);
% %   [M,q,dl,dr] = hnm4lcp_biscale (M,q,2,5.e-1,options);
%     [M,q,dl,dr] = hnm4lcp_biscale (M,q,2,1.e-1,options);
% %   [M,q,dl,dr] = hnm4lcp_biscale (M,q,2,1.e-2,options);
% %   [M,q,dl,dr] = hnm4lcp_biscale (M,q,2,1.e-3,options);
% %   [M,q,dl,dr] = hnm4lcp_biscale (M,q,2,1.e-4,options);
% %   [M,q,dl,dr] = hnm4lcp_biscale (M,q,Inf,1.e-0,options);
% %   [M,q,dl,dr] = hnm4lcp_biscale (M,q,Inf,1.e-1,options);
% %   [M,q,dl,dr] = hnm4lcp_biscale (M,q,Inf,1.e-2,options);
% %   [M,q,dl,dr] = hnm4lcp_biscale (M,q,Inf,1.e-4,options);
% 
%       if ~isempty(xsol); xsol = xsol./dr; end
% %     y = M*x+q;
% 
%     elseif scaling == 4
% 
%       S = zeros(n,1);
%       for i = 1:n
%         s = 1/norm(M(i,:));
%         M(i,:) = s*M(i,:);
%         S(i) = s;
%       end
%       q = S.*q;
%       sigma  = sum(S)/n;
%       sigma2 = norm(q);
%       if sigma2 ~= 0
%         q = q/sigma2;
%       end
%       fprintf ('\n. Diagonal scaling of M and q by      %8.2e (average value)',sigma);
%       if sigma2 ~= 0; fprintf ('\n. Additional scaling of q by          %8.2e',1/sigma2); end
%       y = M*x+q;
% 
%     end
% 
%     scaling_time = toc(tic_scaling);

    try
      tic
    % [x,mu] = pathlcp(M,q,[],[],100*rand(n,1));
      [x,mu] = pathlcp(M,q);
      if ~isempty(xsol) & ~isempty(x); fprintf('\n(pathlcp) L2 relative error on the solution = %11.5e',norm(x-xsol)/norm(xsol)); end
      fprintf('\n(pathlcp) CPU time = %g sec\n',toc);
    catch ME
      fprintf('\n### Pathlcp fails to solve the problem');
      ME_message = ME.message;
      if strcmp(ME_message(1:19),'Solution not finite')
        fprintf(' (the solution is not finite)\n');
      elseif strcmp(ME_message(1:19),'Path fails to solve problem')
        fprintf(' (unspecified reason)\n');
      else
        fprintf(' (not recognoized reason)\n');
      end
    end

  end

%---------------------------------------
% Solver LCPsolve
%---------------------------------------

  if solver.lcpsolve

    tic
    [y,x,retcode] = LCPsolve(M,q);
    if retcode(1) == 1
      fprintf('\nSolution found\nComputing time = %g\n',toc);
    elseif retcode(1) == 2
      fprintf('\nRay termination\nComputing time = %g\n',toc);
    elseif retcode(1) == 3
      fprintf('\nMax iteration reached\nComputing time = %g\n',toc);
    end

  end

return
