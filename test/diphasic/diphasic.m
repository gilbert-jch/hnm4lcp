function [M,q,A,a,B,b,C,c] = diphasic(file,flag,verb)

%%
% [M,q,A,a,B,b,C,c] = diphasic(file,flag,verb)
%
% flag is
% 0: export only the matrix M and the vector q
% 1: export also the matrices A,B,C and the vectors a,b,c
%
% verb
% == 0: work silently
% ~= 0: printings is allowed

%x  = [];
%x1 = [];

  A = [];
  B = [];
  C = [];
  a = [];
  b = [];
  c = [];

% Introductory message

  cl = clock;
  dstr = sprintf('%i-%i-%i, %i:%i:%i',cl(3),cl(2),cl(1),cl(4),cl(5),fix(cl(6)));
  if verb; fprintf('\nDiphasic problem ''%s'' (%s)',file,dstr); end

% Read data in the file

  fmat = fopen(file,'r');
% n = fread(fmat,1,'int32');
% fprintf('\n. n = %i\n',n);
% q = fread(fmat,n,'double');
% M = fread(fmat,[n,n],'double');
  n = fscanf(fmat,'%i',1);
  fgets(fmat);
  if verb; fprintf('\n. n = %i\n',n); end
  q = fscanf(fmat,'%e',n);
  fscanf(fmat,'%i',1);
  fgets(fmat);
  M = fscanf(fmat,'%e',[n,n]);
  if (flag==1)
      % vecteur a et matrice A => H(x) et dH(x)
      fgets(fmat);
      m = fscanf(fmat,'%i',1);
      fgets(fmat);
      %fprintf('\n. size(H(x))  = %i\n',m);
      a = fscanf(fmat,'%e',m);
      fscanf(fmat,'%i',1);
      fgets(fmat);
      %fprintf('\n. size(dH(x)) = %i  %i\n',m,n+m);
      A= fscanf(fmat,'%e',[m,n+m]);
      fgets(fmat);
      % vecteur b et matrice B => G(x) et dG(x)
      p = fscanf(fmat,'%i',1);
      fgets(fmat);
      %fprintf('\n. size((G(x)) = %i\n',p);
      b = fscanf(fmat,'%e',p);
      fscanf(fmat,'%i',1);
      fgets(fmat);
      %fprintf('\n. size(dG(x)) = %i  %i\n',p,n+m);
      B= fscanf(fmat,'%e',[p,n+m]);
      fgets(fmat);
      % vecteur c et matrice C => F(x) et dF(x)
      k = fscanf(fmat,'%i',1);
      fgets(fmat);
      %fprintf('\n. size((F(x)) = %i\n',k);
      c = fscanf(fmat,'%e',k);
      fscanf(fmat,'%i',1);
      fgets(fmat);
      %fprintf('\n. size(dF(x)) = %i  %i\n',k,n+m);
      C= fscanf(fmat,'%e',[k,n+m]);
      fgets(fmat);
  end

  fclose(fmat);

  return
