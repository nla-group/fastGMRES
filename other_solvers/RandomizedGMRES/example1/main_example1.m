clear
clc
f=@(x,y) sin(10*(x+y))./(cos(100*(y-x))+1.1);
n=1e5;
m=250;
x=linspace(0,1,n)';
y=linspace(0,1,m)';
[X,Y]=meshgrid(x,y);

W=f(X,Y)';

% Perform QR decomposition on W w.r.t GS process
[Q_rgs_t1,R_rgs_t1,iter_qr,n_cond_rgs_t1,n_condS_rgs_t1]=RGS_cond(W,200);
[Q_rgs_t2,R_rgs_t2,~,n_cond_rgs_t2,n_condS_rgs_t2]=RGS_cond(W,1000);
[Q_rgs_t3,R_rgs_t3,~,n_cond_rgs_t3,n_condS_rgs_t3]=RGS_cond(W,5000);
[Q_mgs,R_mgs,~,n_cond_mgs]=MGS_cond(W);
[Q_cgs,R_cgs,~,n_cond_cgs]=CGS_cond(W);

%%Graphics
% Stability figure
figure(1)
semilogy(iter_qr,n_cond_cgs,'-*',iter_qr,n_cond_mgs,'-o',...
    iter_qr,n_cond_rgs_t1,'-s',iter_qr,n_cond_rgs_t2,'-d',...
    iter_qr,n_cond_rgs_t3,'-<',iter_qr,n_condS_rgs_t1,'-+',...
    iter_qr,n_condS_rgs_t2,'-.',iter_qr,n_condS_rgs_t3,'-^');
legend('CGS Q','MGS Q','RGS Q, t=200','RGS Q, t=1000',...
    'RGS Q, t=5000','RGS S, t=200','RGS S, t=1000','RGS S, t=5000',...
    'Location','best')
xlim([0,250])
ylim([0,1e5])
title('Condition numbers of $$Q_i$$ and $$S_i$$','interpreter','latex' )
xlabel('iterations $$i$$','interpreter','latex')

% Approximation error
error_rgs_t1=zeros(1,25);
error_rgs_t2=zeros(1,25);
error_rgs_t3=zeros(1,25);
error_mgs=zeros(1,25);
error_cgs=zeros(1,25);

for i=1:25
    normW=norm(W(:,1:10*i));
    idx=1:10*i;
    error_rgs_t1(i)=norm(W(:,idx)-Q_rgs_t1(:,idx)*R_rgs_t1(idx,idx))/normW;
    error_rgs_t2(i)=norm(W(:,idx)-Q_rgs_t2(:,idx)*R_rgs_t2(idx,idx))/normW;
    error_rgs_t3(i)=norm(W(:,idx)-Q_rgs_t3(:,idx)*R_rgs_t3(idx,idx))/normW;
    error_mgs(i)=norm(W(:,idx)-Q_mgs(:,idx)*R_mgs(idx,idx))/normW;
    error_cgs(i)=norm(W(:,idx)-Q_cgs(:,idx)*R_cgs(idx,idx))/normW;
end

figure(2)
semilogy(iter_qr,error_cgs,'-*',iter_qr,error_mgs,'-o',...
    iter_qr,error_rgs_t1,'-s',iter_qr,error_rgs_t2,'-d',...
    iter_qr,error_rgs_t3,'-<');
xlim([0,250])
legend('CGS','MGS','RGS, t=200','RGS, t=100','RGS, t=5000')
title('Error $$|| W_i-Q_iR_i||/|| W_i||$$','interpreter','latex')
xlabel('iterations $$i$$','interpreter','latex')

