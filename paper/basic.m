clc
load problems/vas_stokes_1M.mat
A = Problem.A;
[L,U] = ilu(A);
b = sum(A,2);
tol = 1e-6;

%%
disp('full gmres - 500 iterations max')
tic
x0 = gmres(A,b,[],tol,500,L,U);
toc                  % 127 seconds
norm(b-A*x0)/norm(b) % 9.9505e-04 (converged)

%%
disp('gmres(10) - 200 restarts max')
tic
x1 = gmres(A,b,10,tol,200,L,U);
toc                  % 124 seconds
norm(b-A*x1)/norm(b) % 0.0303 (not converged)

%%
disp('fastgmres')
rng('default')
tic
[x2,resvec,restime] = fastgmres(A,b,tol,100,L,U);
toc                  % 34 seconds
norm(b-A*x2)/norm(b) % 9.2501e-04 (converged)
