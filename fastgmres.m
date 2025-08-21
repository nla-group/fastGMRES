function [x,resvec,restime] = fastgmres(A,b,tol,maxit,M1,M2,x0,opts)
% FASTGMRES   Flexible and sketched GMRES with Arnoldi truncation
%
% [x,resvec,restime] = fastgmres(A,b,tol,maxit,M1,M2,x0,opts)
%
% Solves A*x = b to a relative residual tolerance tol.
% Maximal number of outer FGMRES iterations is maxit (default 100).
% Preconditioners M1 and M2 can be provided, effectively solving
% the problem  M2\(M1\(A*x)) = M2\(M1\b). In this case, tol will
% be the relative target residual norm for the precond. problem.
% An initial guess x0 can be provided as well. 
%
% Optional input structure 
%   opts.verbose (0 = turn off info, [1] = some info, 2 = all info)
%   opts.maxtime (maximal runtime of the outer FGMRES loop, [inf])
%   opts.t (Arnoldi truncation, default [0])
%   opts.Sfun (sketching function)
%   opts.m (maximal number of inner sGMRES iteration, default [500])
%   opts.s (embedding dimension, default [1000])
%
% The outputs are x (the approximate solution), and two vectors 
% with the absolute residual and runtime after each FGMRES iter.
%
% Reference: https://arxiv.org/abs/2506.18408

tstart = tic;
restime = []; % time since tstart when residual is measured
resvec = [];  % residual norms of flexible GMRES

if nargin==8 && isfield(opts,'verbose') 
    verbose = opts.verbose;
else
    verbose = 1;
end
if nargin==8 && isfield(opts,'maxtime') % timeout in seconds
    maxtime = opts.maxtime;
else
    maxtime = inf;
end
if nargin==8 && isfield(opts,'t') % Arnoldi truncation param for sGMRES
    t = opts.t;
else
    t = 0;
end
if nargin==8 && isfield(opts,'Sfun') % function returning sketching handle
    Sfun = opts.Sfun;                % must be callable as Sfun(N,s)
else
    Sfun = @(N,s) clarkson_woodruff_inner(N,s);
end
if nargin==8 && isfield(opts,'m') % max nr of inner sGMRES iterations
    m = opts.m;
else
    m = 500;
end
if nargin==8 && isfield(opts,'s') % embedding dimension
    s = opts.s;
else
    s = 2*m;
end
assert(m<s,'embedding dimension s must be > m inner iterations')
if nargin==8 && isfield(opts,'cndtol') % basis condition number tol
    cndtol = opts.cndtol;
else
    cndtol = 1e15;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isnumeric(A)
    A = @(x) A*x;
end

N = size(b,1);
if nargin >= 7 && ~isempty(x0)
    r0 = b - A(x0);
else
    x0 = zeros(N,1);
    r0 = b;
end

if nargin >=5 && ~isempty(M1) && ~isempty(M2)
    hA = @(x) M2\(M1\(A(x)));
    r0 = M2\(M1\r0);
    rhs = M2\(M1\b); 
    nrmrhs = norm(rhs); % norm of precond rhs.
else
    hA = A;
    rhs = b;
    nrmrhs = norm(rhs);
end

if nargin < 4
    maxit = 100;
end

%%
% r0 = b - A*x0
% solve A*xm ~ r0 with rm' = r0 - A*xm
% then rm = b - A*(x0+xm) = b - A*x0 - A*xm = r0 - Axm = rm'

W = zeros(N,maxit+1);
Z = zeros(N,maxit);
H = zeros(maxit+1,maxit);
nrmr0 = norm(r0);
W(:,1) = r0/nrmr0;
%solver = @(rhs) sgmres(hA,rhs,m,cndtol,t);
solver = @(rhs) sgmres_inner(hA, rhs, m, cndtol, t, s, Sfun, verbose);

for j = 1:maxit

    if verbose >=2, fprintf('  FGMRES iteration = %d\n',j); end

    wvec = W(:,j);
    z = solver(wvec);
    Z(:,j) = z;      % store
    w = hA(z);       % matvec

    if 0 % currently not using inner residual estimator
        pcres(j) = norm(wvec - w);
        if j > 1
            gam = H(1:j-1,1:j-1)\(nrmb*eye(j-1,1)); % FOM coeffs of previous iter
            bnd(j) = H(j,j-1)*abs(gam(end))*pcres(j);
        else
            bnd(1) = nrmb*pcres(1);  % nrmb = \|r_0\|
        end
    end

    for i = 1:j
        H(i,j) = W(:,i)'*w;
        w = w - W(:,i)*H(i,j);
    end
    % CGS
    %H(1:j,j) = W(:,1:j)'*w;
    %w = w - W(:,1:j)*H(1:j,j);

    H(j+1,j) = norm(w);
    W(:,j+1) = w/H(j+1,j);

    % solve LS problem for FGMRES
    y = H(1:j+1,1:j)\(nrmr0*eye(j+1,1));

    % get true solution x and true residual of unprecond problem
    %x = x0 + Z(:,1:j)*y;
    %resvec(j) = norm(b - A(x));

    % this is a cheaper way to compute the precond residual without x
    resvec(j) = norm(nrmr0*eye(j+1,1) - H(1:j+1,1:j)*y);

    if verbose >= 2, fprintf('  FGMRES residual  = %5.3e\n\n',resvec(j)); end
    restime(j) = toc(tstart);

    if resvec(j)/nrmrhs < 0.99*tol || restime(j) > maxtime
        break
    end
end

% form solution and validate true (precond.) residual
x = x0 + Z(:,1:j)*y;
res_true = norm(rhs - hA(x))/nrmrhs;
if res_true > tol
    warning('Target residual NOT reached after %d iterations: %5.3e > %5.3e \n         Consider increasing maxit or restart with x as initial guess.',j,res_true,tol);
end
if verbose && res_true <= tol
    fprintf('fastgmres reached target residual after %d iterations: %5.3e <= %5.3e \n',j,res_true,tol);
end

end

function [x,res,sres,cnd,scnd] = sgmres_inner(A, b, m, cndtol, t, s, Sfun, verbose)

N = size(b,1);
hS = Sfun(N,s);

res = norm(b); sres = []; cnd = []; scnd = [];
H = zeros(m+1,m); V = zeros(N,m+1);
v = b/norm(b); V(:,1) = v;
SAV = zeros(s,m+1);
SV = zeros(s,m+1);
SV(:,1) = hS(V(:,1));
Sb = hS(b);
Q = []; R = [];
if verbose >= 2, fprintf('  running sgmres '); end
for j = 1:m

    if verbose >= 2, fprintf('.'); end
    w = A(v);
    SAV(:,j) = hS(w);
    %[Q,R] = qr(SAV(:,1:j),0);
    [Q,R] = qrupdate_gs_inner(SAV(:,1:j),Q,R);

    if j > 1 && cond(R) > cndtol
        break
    end

    y = R\(Q'*Sb);   % sketched GMRES
    %sres(j+1) = norm(Sb - SAV(:,1:j)*y);

    for i = max(j-t+1,1):j
        H(i,j) = V(:,i)'*w;
        w = w - V(:,i)*H(i,j);
    end
    H(j+1,j) = norm(w);
    v = w/H(j+1,j);
    V(:,j+1) = v;

end
if verbose >= 2, fprintf('\n'); end
x = V(:,1:length(y))*y;
end % sgmres_inner

function S = clarkson_woodruff_inner(N, s)

i = randi(s,N,1);
j = 1:N;
el = 2*randi(2,N,1)-3;
S = sparse(i,j,el,s,N,N);
S = @(x) S*x;

end

function [Q1,R1] = qrupdate_gs_inner(A,Q,R)
% modified Gram-Schmidt update of QR factorization
% Example use:
% 
%  A = randn(10,6);
%  [Q,R] = qr(A(:,1:4),0); % partial QR
%  [Q1,R1] = qrupdate_gs(A,Q,R);
[m,n] = size(A);
k = size(Q,2);
assert(m>n & k<n, 'A needs to be tall and size(Q,2)<size(A,2)')
Q1 = zeros(m,n);
R1 = zeros(n,n);
Q1(1:m,1:k) = Q;
R1(1:k,1:k) = R;
for j = k+1:n
    w = A(:,j);
    for reo = 0:1
        for i = 1:j-1
            r = Q1(:,i)'*w;
            R1(i,j) = R1(i,j) + r;
            w = w - Q1(:,i)*r;
        end
    end
    R1(j,j) = norm(w);
    Q1(:,j) = w/R1(j,j);
end

end % qrupdate_gs_inner