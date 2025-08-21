# fastGMRES [![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=nla-group/fastgmres)

Flexible and sketched GMRES with truncated Arnoldi 

These are the MATLAB files accompanying the paper https://arxiv.org/abs/2506.18408

Basic usage example (requires the matrix [vas_stokes_1M](https://sparse.tamu.edu/VLSI/vas_stokes_1M)):

```Matlab
clc
load vas_stokes_1M.mat
A = Problem.A;
b = sum(A,2);
tol = 1e-3;

%%
disp('full gmres - 500 iterations max')
tic
x0 = gmres(A,b,[],tol,500);
toc                  % 127 seconds
norm(b-A*x0)/norm(b) % 9.9505e-04 (converged)

%%
disp('gmres(10) - 200 restarts max')
tic
x1 = gmres(A,b,10,tol,200);
toc                  % 124 seconds
norm(b-A*x1)/norm(b) % 0.0303 (not converged)

%%
disp('fastgmres')
rng('default')
tic
[x2,resvec,restime] = fastgmres(A,b,tol);
toc                  % 34 seconds
norm(b-A*x2)/norm(b) % 9.2501e-04 (converged)
```

(Timings in MATLAB Online, August 2025, running on Intel(R) Xeon(R) Platinum 8375C CPU @ 2.90GHz.)
