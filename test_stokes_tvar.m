% Exploring different Arnoldi truncation parameters for the precond. Stokes
% problem.

mydefaults
clear all
close all
clc
rng('default')
load 'problems/vas_stokes_1M.mat'
A = Problem.A;
N = size(A,1);
b = randn(N,1);
[L,U] = ilu(A);

tol = 1e-13; 
maxit = 90;
x0 = zeros(N,1);

%%
for t = 0:3
    fprintf('\n--------------------------\nArnoldi truncation t = %d\n',t)
    opts.t = t;
    [x,resvec,restime] = fastgmres(A,b,tol,maxit,L,U,x0,opts);
    semilogy(restime,resvec,'-o','LineWidth',2)
    hold on
    shg
end

xlabel('time in seconds','interpreter','latex')
ylabel('residual norm','interpreter','latex')
xlim([0,500])
ylim([1e-12,1e4]), shg
legend('fastGMRES($t = 0$)','fastGMRES($t = 1$)','fastGMRES($t = 2$)','fastGMRES($t = 3$)','Location','southwest','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
mypdf('stokes_tvar',0.6,1.0)

%% GMRES-SDR
% Code from https://github.com/nla-group/GMRES-SDR
addpath('other_solvers')
param.k = 20;
param.max_it = 100;
param.max_restarts = 20;
param.t = 2;       % Arnoldi truncation parameter
param.s = 1000;
param.hS = clarkson_woodruff(N, 1000);
param.pert = 0;    % matrix A stays constant
param.verbose = 2; % no debug info computed/printed
rng('default')     % Re-initialize for randomized sketching
param.U = []; param.SU = []; param.SAU = [];

Ap = @(x) U\(L\(A*x));
bp = U\(L\b);
param.tol = tol*norm(bp);
[x,out] = gmres_sdr(Ap,bp,param);
semilogy(out.timevec(2:end), out.residuals(2:end),'-s','LineWidth',2)
shg
legend('fastGMRES($t = 0$)','fastGMRES($t = 1$)','fastGMRES($t = 2$)','fastGMRES($t = 3$)','GMRES-SDR(100,20)','Location','southwest','interpreter','latex','NumColumns',1,'Box','off')
mypdf('stokes_tvar',0.6,1.0)
