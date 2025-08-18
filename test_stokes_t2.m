% Stokes problem with ILU preconditioning, using maximal Krylov basis
% size of 500 for sGMRES (never reached), and halting sGMRES when
% cond(SAV) > cndtol = 1e15.
%
% This illustrates that FGMRES is better than restarting sGMRES(cndtol). 
%
% Uses Arnoldi truncation parameter t=2.

mydefaults
clear all
close all
clc
addpath('other_solvers')
rng('default')
load 'problems/vas_stokes_1M.mat'
A = Problem.A;
N = size(A,1);

[L,U] = ilu(A);
A = @(x) U\(L\(A*x));

b = randn(N,1); 
b = U\(L\b);
%b = b/norm(b);

%%
maxit = 30;
W = zeros(N,maxit+1);
Z = zeros(N,maxit);
H = zeros(maxit+1,maxit);
nrmb = norm(b);
W(:,1) = b/nrmb;
FGMRES = []; % residuals of flexible GMRES
FFOM = [];   % residuals of flexible FOM

solver = @(rhs) sgmres(A,rhs,500,1e15,2);

tstart = tic;
TIME_FGMRES = [];

for j = 1:maxit

    fprintf('FGMRES iteration   = %d\n',j)

    rhs = W(:,j);
    z = solver(rhs);    
    Z(:,j) = z;      % store
    w = A(z);        % matvec
    pcres(j) = norm(rhs - w); 

    if j > 1
        gam = H(1:j-1,1:j-1)\(nrmb*eye(j-1,1)); % FOM coeffs of previous iter
        bnd(j) = H(j,j-1)*abs(gam(end))*pcres(j);
    else
        bnd(1) = nrmb*pcres(1);  % nrmb = \|r_0\|
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
    y = H(1:j+1,1:j)\(nrmb*eye(j+1,1));
    x = Z(:,1:j)*y;

    FGMRES(j) = norm(b - A(x));
    %FGMRES2(j) = norm(nrmb*eye(j+1,1) - H(1:j+1,1:j)*y);
    fprintf('  FGMRES residual = %5.3e\n',FGMRES(j))

    % NOTE: We're computing the true residuals for the plots;
    % The other expressions are easier to evaluate but they can be
    % misleading close to stagnation. This is also for fairness 
    % to the other methods.

    %norm(A*Z(:,1:j) - W(:,1:j+1)*H(1:j+1,1:j))

    TIME_FGMRES(j) = toc(tstart);
end

figure
semilogy(FGMRES,'-o','LineWidth',2), hold on
shg

%%
tstart = tic;
TIME_FFOM = [];

for j = 1:maxit

    fprintf('FFOM iteration   = %d\n',j)

    rhs = W(:,j);
    z = solver(rhs);
    Z(:,j) = z;      % store
    w = A(z);        % matvec
    pcres(j) = norm(rhs - w); 

    for i = 1:j
        H(i,j) = W(:,i)'*w;
        w = w - W(:,i)*H(i,j);
    end
    % CGS
    %H(1:j,j) = W(:,1:j)'*w;
    %w = w - W(:,1:j)*H(1:j,j);

    H(j+1,j) = norm(w);
    W(:,j+1) = w/H(j+1,j);

    % solve linear system for FFOM
    y_fom = H(1:j,1:j)\(nrmb*eye(j,1));
    x_fom = Z(:,1:j)*y_fom;
    FFOM(j) = norm(b - A(x_fom));
    %FFOM2(j) = abs(H(j+1,j)*y_fom(end));

    fprintf('  FFOM residual   = %5.3e\n',FFOM(j))

    % NOTE: We're computing the true residuals for the plots;
    % The other expressions are easier to evaluate but they can be
    % misleading close to stagnation. This is also for fairness 
    % to the other methods.

    %norm(A*Z(:,1:j) - W(:,1:j+1)*H(1:j+1,1:j))

    TIME_FFOM(j) = toc(tstart);
end

semilogy(FFOM,'-*','LineWidth',2)
shg

%%
x = zeros(N,1); r = b - A(x);

tstart = tic;
TIME_REST = [];

for j = 1:maxit
    fprintf('restart cycle   = %d\n',j)

    z = solver(r);
    x = x + z;
    r = b - A(x);
    REST(j) = norm(r);
    fprintf('  RESTART residual = %5.3e\n',REST(j))

    TIME_REST(j) = toc(tstart);
end

semilogy(REST,'-x','LineWidth',2)
semilogy(bnd,'k--','LineWidth',2)
xlabel('outer iteration $j$','interpreter','latex')
ylabel('residual norm','interpreter','latex')
ylim([1e-13,1e4]), shg
legend('FGMRES with sGMRES($10^{15}$)','FFOM with sGMRES($10^{15}$)','restarted sGMRES($10^{15}$)','FGMRES bound','Location','southwest','interpreter','latex','NumColumns',1,'Box','off')
set(gca,'TickLabelInterpreter','latex')
mypdf('stokes_t2_res',0.6,1.0)

%%
figure
semilogy(TIME_FGMRES,FGMRES,'-o','LineWidth',2), hold on
semilogy(TIME_FGMRES,FFOM,'-*','LineWidth',2)
semilogy(TIME_REST,REST,'-x','LineWidth',2)
xlabel('time in seconds','interpreter','latex')
ylabel('residual norm','interpreter','latex')
xlim([0,700])
ylim([1e-12,1e4]), shg
legend('FGMRES with sGMRES($10^{15}$)','FFOM with sGMRES($10^{15}$)','restarted sGMRES($10^{15}$)','Location','southwest','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
mypdf('stokes_t2_time',0.6,1.0)


%% GMRES-SDR
% Code from https://github.com/nla-group/GMRES-SDR

param.k = 20;
param.max_it = 100;
param.max_restarts = 20;
param.tol = 0;
param.t = 2;       % Arnoldi truncation parameter
param.s = 1000;
param.hS = clarkson_woodruff(N, 1000);
param.pert = 0;    % matrix A stays constant
param.verbose = 2; % no debug info computed/printed
rng('default')     % Re-initialize for randomized sketching
param.U = []; param.SU = []; param.SAU = [];
[x,out] = gmres_sdr(A,b,param);

semilogy(out.timevec(2:end), out.residuals(2:end),'-s','LineWidth',2)
legend('FGMRES with sGMRES($10^{15}$)','FFOM with sGMRES($10^{15}$)','restarted sGMRES($10^{15}$)','GMRES-SDR(100,20)','Location','southwest','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
shg
mypdf('stokes_t2_time',0.6,1.0)

%% RGS_GMRES_mdr
addpath('other_solvers/RandomizedGMRES/example2/')

% This is the code from https://github.com/Yongseok7717/RandomizedGMRES
% using the functions in the 'example2' subfolder

%[x,residual,iter]=RGS_GMRES_mdr(A,b,tol,m,k,t,M1,M2)
%   t sketching dimension
%   k deflated harmonic Ritz vectors between FGMRES restarts
%   m FGMRES basis to be computed
%
% WE'VE MADE THESE CHANGES: 
% - there were at most 500 total iterations (hard coded), 
%   changed to 20 outer GMRES restarts
% - we have changed every occurence of A*x to A(x), 
%   so that the code can work with function handles, and provide the
%   preconditioner in form of a modified matrix to be consistent with 
%   our other codes
% - the sketching operator theta was implemented as a dense matrix 
%   in RGS_GMRES_mdr; changed to using srft for better performance
% - residual norm output to absolute residuals
% 
% (???It looks like RGS_GMRES_mdr uses m-k vectors from the second cycle;
% that may be different from GMRES-SDR.???)
% 

L = []; U = []; % preconditioner already provided with A
tol = 1e-16;
[x,Rres_mdr,Riter_mdr,timevec] = RGS_GMRES_mdr_maxit20(A,b,tol,100,20,1000,L,U);
semilogy(timevec(2:end), Rres_mdr(2:end),'-d','LineWidth',2)
legend('FGMRES w/ sGMRES($10^{15}$)$~~$','FFOM w/ sGMRES($10^{15}$)','restarted sGMRES($10^{15}$)','GMRES-SDR(100,20)','RGS-GMRES-MDR(100,20)','Location','northeast','interpreter','latex','NumColumns',1,'Box','off')
shg
xlim([0,600])
set(gca,'TickLabelInterpreter','latex')
mypdf('stokes_t2_time',0.6,1.0)

%save stokes_t2.mat