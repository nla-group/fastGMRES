mydefaults
clear all
close all
clc
rng('default')
addpath('other_solvers')

%%
description = 't2em without precond, b = randn';
load 'problems/t2em.mat'
A = Problem.A;
N = size(A,1);
b = randn(N,1); b = b/norm(b);
hA = @(x) A*x;
tol = 1e-6; 
x0 = zeros(N,1);
maxtime = 300; % maximal runtime
maxit = 1e6;   % all methods terminated by timeout or convergence

%%
tstart_total = tic;
mylog(sprintf('\n---------------------------------------------------------------------------\n'))
mylog(mfilename)
mylog(sprintf('\n'))
mylog(string(datetime))
mylog(' -- ')
mylog(string(feature('GetCPU')))
mylog(sprintf('\n'))
mylog(description)
mylog(sprintf('\n'))
mylog(sprintf('tol: %3.1e, maxtime: %3.1f\n\n',tol,maxtime))

%%
% fastGMRES
mylog(sprintf('fastgmres\n'))
opts.maxtime = maxtime;
opts.verbose = 0;
tstart = tic;
x = x0; % reset to detect failure
[x,resvec,restime] = fastgmres(hA,b,tol,500,[],[],x0,opts);
truetime = toc(tstart);
trueres = norm(A*x - b)/norm(b);
semilogy(restime,resvec,'-o')
mylog(sprintf('  residual: %3.1e, time: %3.1f\n\n',trueres,truetime))
hold on
shg

%%
%maxtime = 2*truetime; % new maxtime = twice that of fastgmres
%mylog(sprintf('resetting maxtime to %3.1f\n\n',maxtime))

%% bicgstab
mylog(sprintf('bicgstab\n'))
opts = struct();
opts.reltol = tol; 
opts.maxit = maxit;
opts.maxtime = maxtime;
tstart = tic;
x = x0;
[x,out] = mybicgstab(hA,b,opts);
truetime1 = toc(tstart);
trueres1 = norm(A*x - b)/norm(b);
semilogy(out.restime, out.resvalue(1:length(out.restime)),'--')
mylog(sprintf('  residual: %3.1e, time: %3.1f\n\n',trueres1,truetime1))
shg

%% gmres(50)
mylog(sprintf('gmres(50)\n'))
restart = 50; 
tstart = tic;
x = x0;
[x,flag,relres,iter,resvec,restime] = mygmres(hA,b,restart,tol,maxit,[],[],maxtime);
truetime3 = toc(tstart);
trueres3 = norm(A*x - b)/norm(b);
semilogy(restime(1:restart:end),resvec(1:restart:end),'-s')
mylog(sprintf('  residual: %3.1e, time: %3.1f\n\n',trueres3,truetime3))
shg

% gmres(100)
mylog(sprintf('gmres(100)\n'))
restart = 100; 
tstart = tic;
x = x0;
[x,flag,relres,iter,resvec,restime] = mygmres(hA,b,restart,tol,maxit,[],[],maxtime);
truetime4 = toc(tstart);
trueres4 = norm(A*x - b)/norm(b);
semilogy(restime(1:restart:end),resvec(1:restart:end),'-s')
mylog(sprintf('  residual: %3.1e, time: %3.1f\n\n',trueres4,truetime4))
shg

%% gmres_dr(50,10)
mylog(sprintf('gmres_dr(50,10)\n'))
opts = struct();
opts.tol = tol/norm(b); 
opts.k = 10;
opts.max_it = 50; % careful, max_it is restart length here
opts.max_restarts = maxit;
opts.maxtime = maxtime;
tstart = tic;
x = x0;
[x,out] = gmres_dr(hA,b,opts);
truetime5 = toc(tstart);
trueres5 = norm(A*x - b)/norm(b);
semilogy(out.restime(1:opts.max_it:end), out.resvec(1:opts.max_it:end),'d-')
mylog(sprintf('  residual: %3.1e, time: %3.1f\n\n',trueres5,truetime5))
shg

% gmres_dr(100,20)
mylog(sprintf('gmres_dr(100,20)\n'))
opts = struct();
opts.tol = tol/norm(b); 
opts.k = 20;
opts.max_it = 100; % careful, max_it is restart length here
opts.max_restarts = maxit;
opts.maxtime = maxtime;
tstart = tic;
x = x0;
[x,out] = gmres_dr(hA,b,opts);
truetime6 = toc(tstart);
trueres6 = norm(A*x - b)/norm(b);
semilogy(out.restime(1:opts.max_it:end), out.resvec(1:opts.max_it:end),'d-')
mylog(sprintf('  residual: %3.1e, time: %3.1f\n\n',trueres6,truetime6))
shg

%% gmres_sdr(50,10)
mylog(sprintf('gmres_sdr(50,10)\n'))
param = struct();
param.k = 10;
param.max_it = 50; % careful, max_it is restart length here
param.maxtime = maxtime;
param.max_restarts = maxit;
param.tol = tol/norm(b);
param.t = 2;       % Arnoldi truncation parameter
param.s = 1000;
param.hS = clarkson_woodruff(N, 1000);
param.pert = 0;    % matrix A stays constant
param.verbose = 0; % no debug info computed/printed
rng('default')     % Re-initialize for randomized sketching
param.U = []; param.SU = []; param.SAU = [];
param.maxtime = maxtime;
tstart = tic;
x = x0;
[x,out] = gmres_sdr(hA,b,param);
truetime7 = toc(tstart);
trueres7 = norm(A*x - b)/norm(b);
semilogy(out.timevec, out.residuals,'*-')
mylog(sprintf('  residual: %3.1e, time: %3.1f\n\n',trueres7,truetime7))
shg

%% gmres_sdr(100,20)
mylog(sprintf('gmres_sdr(100,20)\n'))
param = struct();
param.k = 20;
param.max_it = 100; % careful, max_it is restart length here
param.maxtime = maxtime;
param.max_restarts = maxit;
param.tol = tol/norm(b);
param.t = 2;       % Arnoldi truncation parameter
param.s = 1000;
param.hS = clarkson_woodruff(N, 1000);
param.pert = 0;    % matrix A stays constant
param.verbose = 0; % no debug info computed/printed
rng('default')     % Re-initialize for randomized sketching
param.U = []; param.SU = []; param.SAU = [];
param.maxtime = maxtime;
tstart = tic;
x = x0;
[x,out] = gmres_sdr(hA,b,param);
truetime8 = toc(tstart);
trueres8 = norm(A*x - b)/norm(b);
semilogy(out.timevec, out.residuals,'*-')
mylog(sprintf('  residual: %3.1e, time: %3.1f\n\n',trueres8,truetime8))
shg

%%
mylog(sprintf('Done. Total runtime: %5.1f seconds\n\n',toc(tstart_total)))
mylog(sprintf('%s & %3.1f (%3.1e) & %3.1f (%3.1e) & %3.1f (%3.1e) & %3.1f (%3.1e) & %3.1f (%3.1e) & %3.1f (%3.1e) & %3.1f (%3.1e) & %3.1f (%3.1e)  \\\\ \n\n',mfilename,truetime,trueres,truetime1,trueres1,truetime3,trueres3,truetime4,trueres4,truetime5,trueres5,truetime6,trueres6,truetime7,trueres7,truetime8,trueres8))

%%
legend('fastgmres','bicgstab','gmres(50)','gmres(100)','gmres-dr(50,10)','gmres-dr(100,20)','gmres-sdr(50,10)','gmres-sdr(100,20)')
title(description)
xlabel('runtime in seconds')
ylabel('residual norm')
axis tight
shg
mypdf(mfilename,0.6,1.0)