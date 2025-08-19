# fastGMRES
Flexible and sketched GMRES with truncated Arnoldi

These are the MATLAB files accompanying the paper https://arxiv.org/abs/2506.18408

Basic usage example (requires the matrix [vas_stokes_1M](https://sparse.tamu.edu/VLSI/vas_stokes_1M)):

```Matlab
load vas_stokes_1M.mat
A = Problem.A;
b = sum(A,2);
tol = 1e-3;

%%
disp('full gmres - 500 iterations max')
tic
x0 = gmres(A,b,[],tol,500);
norm(b-A*x0)/norm(b) % 9.9505e-04 (converged)
toc                  % 127 seconds

%%
disp('gmres(10) - 200 restarts max')
tic
x1 = gmres(A,b,10,tol,200);
norm(b-A*x1)/norm(b) % 0.0303 (not converged)
toc                  % 132 seconds

%%
disp('fastgmres')
tic
[x2,resvec,restime] = fastgmres(A,b,tol);
norm(b-A*x2)/norm(b) % 9.2365e-04 (converged)
toc                  % 34 seconds
```
