% Naca12 problem with ILU preconditioning, using maximal Krylov size 500 
% Using sGMRES with the cond(SAV)<1e15 condition.
%

mydefaults
clear all
close all
clc
rng('default')

load('problems/naca12A.mat');
N = length(A);
%b = load('problems/RHS.dat');
[L,U] = ilu(A);
A = @(x) U\(L\(A*x));
b = U\(L\b);

%%
for t = 0:3
    fprintf('\n\ntruncation t = %d\n\n',t)
    maxit = 100;
    W = zeros(N,maxit+1);
    Z = zeros(N,maxit);
    H = zeros(maxit+1,maxit);
    nrmb = norm(b);
    W(:,1) = b/nrmb;
    FGMRES = []; % residuals of flexible GMRES
    
    solver = @(rhs) sgmres(A,rhs,500,1e15,t);
    
    tstart = tic;
    TIME_FGMRES = [];
    
    for j = 1:maxit
    
        fprintf('  FGMRES iteration   = %d\n',j)
    
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

        if FGMRES(j)/norm(b) < 1e-13
            break
        end
    end
    
    figure(1)
    semilogy(TIME_FGMRES,FGMRES,'-o','LineWidth',2), hold on
    shg
end

xlabel('time in seconds','interpreter','latex')
ylabel('residual norm','interpreter','latex')
%xlim([0,600])
%ylim([1e-12,1e4]), shg
legend('fastGMRES($t = 0$)','fastGMRES($t = 1$)','fastGMRES($t = 2$)','fastGMRES($t = 3$)','Location','northeast','interpreter','latex')
title('NACA12 with ILU preconditioner','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
mypdf('naca_ilu_time',0.9,0.6)

%% RGS_GMRES_mdr
addpath('other_solvers/RandomizedGMRES/example2')

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
%   in RGS_GMRES_mdr; changed to using srft.
% - residual norm output to absolute residuals
% 

tol = 1e-13;
[x,Rres_mdr,Riter_mdr,timevec] = RGS_GMRES_mdr(A,b,tol,40,20,1000,[],[]);
semilogy(timevec(2:end), Rres_mdr(2:end),'k-','LineWidth',2)
legend('fastGMRES($t = 0$)','fastGMRES($t = 1$)','fastGMRES($t = 2$)','fastGMRES($t = 3$)','RGS-GMRES-MDR(40,20)','Location','northeast','interpreter','latex','NumColumns',1,'Box','off')
xlim([0,12])
set(gca,'TickLabelInterpreter','latex','FontSize',18) % JP: I changed font-size in legend and axes
mypdf('naca_ilu_time',1,0.52)
