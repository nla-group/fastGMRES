% Demonstration of FGMRES deflation effect.

mydefaults
clear all
close all

rng('default')
N = 1000;
%A = randn(N) + 30*eye(N); 
A = diag(linspace(1,N,N));
b = randn(N,1); b = b/norm(b);

maxit = 70;
W = zeros(N,maxit+1);
Z = zeros(N,maxit);
H = zeros(maxit+1,maxit);
nrmb = norm(b);
W(:,1) = b/nrmb;
FGMRES = []; % residuals of flexible GMRES
FFOM = [];   % residuals of flexible FOM

solver = @(rhs) gmres(A,rhs,[],1e-14,5);
x = 0*b;    % initial guess

rgam = [0];  % gamma from recursive bound
rbnd = [nrmb]; % recursive bnd

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

    rgam(j+1) = pcres(j)/sqrt(1-rgam(j)^2);
    rbnd(j+1) = rgam(j+1)*rbnd(j);

    for reo = 0:1
        for i = 1:j
            h = W(:,i)'*w;
            H(i,j) = H(i,j) + h;
            w = w - W(:,i)*h;
        end
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
figure
semilogy(FGMRES); 
hold on
semilogy(REST)
xlabel('outer iteration $j$','interpreter','latex','fontsize',18)
ylabel('residual norm','interpreter','latex','fontsize',18)
legend('FGMRES-GMRES(5)','GMRES(5)','Location','southwest','interpreter','latex','NumColumns',1,'Box','off')
lgd = legend('show'); lgd.FontSize = 18;
axis([0,maxit,1e-16,1])
set(gca,'TickLabelInterpreter','latex','fontsize',18)
mypdf('fgmres_deflation1',1,0.5)

%%
figure
%plot(eig(full(A))+eps*1i,'ko')
for j = 1:maxit
    Aj = pinv(Z(:,1:j))*(A*Z(:,1:j));
    %Am = W'*(A*W);
    ritz = eig(Aj);
    dist = abs(ritz - round(ritz));
    ind1 = find(dist < 1e-3);
    ind2 = find(1e-3 <= dist & dist < 1e-2);
    ind3 = find(1e-2 <= dist & dist < 1e-1);
    ind4 = find(1e-1 <= dist);
    if ~isempty(ind1), plot(j,ritz(ind1),'r.'); end
    if ~isempty(ind2), plot(j,ritz(ind2),'y.'); end
    if ~isempty(ind3), plot(j,ritz(ind3),'g.'); end
    if ~isempty(ind4), plot(j,ritz(ind4),'b.'); end
    hold on
end
axis tight
ylim([0 200])
xlabel('outer iteration $j$','interpreter','latex','fontsize',18)
ylabel('order $j$ Ritz values','interpreter','latex','fontsize',18)
set(gca,'TickLabelInterpreter','latex','fontsize',18)
mypdf('fgmres_deflation2',1,0.5)
