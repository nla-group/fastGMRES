% Demonstration of the FGMRES residual bound.
% Using GMRES(100) as the inner solver.

mydefaults
clear all
close all

rng('default')
N = 1000;
A = randn(N) + 30*eye(N); 
b = randn(N,1); b = b/norm(b);

maxit = 15;
W = zeros(N,maxit+1);
Z = zeros(N,maxit);
H = zeros(maxit+1,maxit);
nrmb = norm(b);
W(:,1) = b/nrmb;
FGMRES = []; % residuals of flexible GMRES
FFOM = [];   % residuals of flexible FOM

solver = @(rhs) gmres(A,rhs,[],1e-14,100);
x = 0*b;    % initial guess

for j = 1:maxit

    rhs = W(:,j);
    z = solver(rhs);

    pcres(j) = norm(rhs - A*z); 
    Z(:,j) = z;      % store
    w = A*z;         % matvec

    if j > 1
        gam = H(1:j-1,1:j-1)\(nrmb*eye(j-1,1));
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
    FGMRES(j) = norm(b - A*x);
    fprintf('FGMRES residual = %5.3e\n',FGMRES(j))

    % solve linear system for FFOM
    y_fom = H(1:j,1:j)\(nrmb*eye(j,1));
    x_fom = Z(:,1:j)*y_fom;
    FFOM(j) = norm(b - A*x_fom);
    fprintf('FFOM residual   = %5.3e\n',FFOM(j))

    %norm(A*Z(:,1:j) - W(:,1:j+1)*H(1:j+1,1:j))

end

%% restarted GMRES
x = zeros(N,1); r = b - A*x;
for j = 1:maxit
    z = solver(r);
    x = x + z;
    r = b - A*x;
    REST(j) = norm(b - A*x);
    fprintf('RESTART residual = %5.3e\n',REST(j))
end

%%
semilogy(FGMRES,'-o','LineWidth',2), hold on
semilogy(FFOM,'-*','LineWidth',2)
semilogy(REST,'-x','LineWidth',2)
semilogy(bnd,'k--','LineWidth',2)
xlabel('outer iteration $j$','interpreter','latex')
ylabel('residual norm','interpreter','latex')
ylim([1e-7,1]), shg
legend('FGMRES','FFOM','GMRES(100)','FGMRES bound','Location','southwest','interpreter','latex','NumColumns',1,'Box','off')
set(gca,'TickLabelInterpreter','latex')
mypdf('fgmres_bnd',0.6,1.0)

