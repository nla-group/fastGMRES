% ML_Geer problem
% Using sGMRES with the cond(SAV)<1e15 condition.
%
% ILU preconditioning

mydefaults
clear all
close all
clc
rng('default')

load('problems/ML_Geer.mat')
A = Problem.A;
N = length(A);
[L,U] = ilu(A);
A = @(x) U\(L\(A*x));
b = randn(N,1);
b = U\(L\b);

%%
for t = 0:3
    fprintf('\n\ntruncation t = %d\n\n',t)
    maxit = 200;
    W = zeros(N,maxit+1);
    Z = zeros(N,maxit);
    H = zeros(maxit+1,maxit);
    nrmb = norm(b);
    W(:,1) = b/nrmb;
    FGMRES = []; % residuals of flexible GMRES

    solver = @(rhs) sgmres(A,rhs,100,1e15,t);

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

        if FGMRES(j)/norm(b) < 1e-6
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
legend('fastGMRES($t = 0$)','fastGMRES($t = 1$)','fastGMRES($t = 2$)','fastGMRES($t = 3$)','Location','northeast','interpreter','latex','NumColumns',1,'Box','off')
title('ML\_Geer with ILU preconditioner','interpreter','latex')
xlim([0,1200])
set(gca,'TickLabelInterpreter','latex','FontSize',18) % JP: I changed font-size in legend and axes
mypdf('ml_geer_ilu_time',1,0.5)

%% Two-level
hold on

maxit = 50;
W = zeros(N,maxit+1);
Z = zeros(N,maxit);
H = zeros(maxit+1,maxit);
nrmb = norm(b);
W(:,1) = b/nrmb;
FGMRES = []; % residuals of flexible GMRES

%solver = @(rhs) sgmres(A,rhs,100,1e15,t);
zer = zeros(N,1);
opts.m = 200;
solver = @(rhs) fastgmres(A,rhs,1e-12,12,[],[],zer,opts);

tstart = tic;
TIME_FGMRES = [];

for j = 1:maxit

    fprintf('OUTER FGMRES iteration   = %d\n',j)

    rhs = W(:,j);
    z = solver(rhs);
    Z(:,j) = z;      % store
    w = A(z);        % matvec

    %
    pcres(j) = norm(rhs - w);
    if j > 1
        gam = H(1:j-1,1:j-1)\(nrmb*eye(j-1,1)); % FOM coeffs of previous iter
        bnd(j) = H(j,j-1)*abs(gam(end))*pcres(j);
    else
        bnd(1) = nrmb*pcres(1);  % nrmb = \|r_0\|
    end
    %

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
    fprintf('OUTER FGMRES residual = %5.3e\n',FGMRES(j))

    % NOTE: We're computing the true residuals for the plots;
    % The other expressions are easier to evaluate but they can be
    % misleading close to stagnation. This is also for fairness
    % to the other methods.

    %norm(A*Z(:,1:j) - W(:,1:j+1)*H(1:j+1,1:j))

    TIME_FGMRES(j) = toc(tstart);

    if FGMRES(j)/norm(b) < 1e-6
        break
    end
end

figure(1)
semilogy(TIME_FGMRES,FGMRES,'k-s','LineWidth',2), hold on
shg
