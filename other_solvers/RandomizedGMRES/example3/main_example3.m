load('LS89.mat');
n=length(A);
b=load('RHS.dat');
tol=1e-10;
t=1000;
theta=Rademacher_r(t,n);
m=100;
minner=15;

[V,H,Z]=inner_outer_Arnoldi_RGS_ras(A,b,m,minner,theta,6,6);
[Vc,Hc,Zc]=CGS_ras(A,b,m,minner,6,6);
[Vm,Hm,Zm]=MGS_ras(A,b,m,minner,6,6);

cond_r=zeros(m+1,1);
cond_m=zeros(m+1,1);
cond_c=zeros(m+1,1);

error_r=zeros(m,1);
error_m=zeros(m,1);
error_c=zeros(m,1);

cond_r(1)=cond(V(:,1));
cond_m(1)=cond(Vm(:,1));
cond_c(1)=cond(Vc(:,1));
    
for i=1:m
    cond_r(i+1)=cond(V(:,1:i+1));
    error_r(i)=norm(A*Z(:,1:i)-V(:,1:i+1)*H(1:i+1,1:i))/norm(A*Z(:,1:i));


    cond_m(i+1)=cond(Vm(:,1:i+1));
    error_m(i)=norm(A*Zm(:,1:i)-Vm(:,1:i+1)*Hm(1:i+1,1:i))/norm(A*Zm(:,1:i));


    cond_c(i+1)=cond(Vc(:,1:i+1));
    error_c(i)=norm(A*Zc(:,1:i)-Vc(:,1:i+1)*Hc(1:i+1,1:i))/norm(A*Zc(:,1:i));
end

% Stability figure
figure(1)
iter=1:m+1;
semilogy(iter,cond_c,'-*',iter,cond_m,'-o',...
    iter,cond_r,'-s');
xlim([0,m]);
legend('CGS','MGS','RGS, t=1000', 'Location','best')
title('Condition numbers of $$V_i$$','interpreter','latex' )
xlabel('iterations $$i$$','interpreter','latex')

% Approximation errors figure
figure(2)
iter=1:m;
semilogy(iter,error_c,'-*',iter,error_m,'-o',iter,error_r,'-s');
xlim([0,m])
legend('CGS','MGS','RGS, t=1000', 'Location','NorthWest')
title('Error $$|| AZ_i-V_{i+1}H_i||/|| AZ_i||$$','interpreter','latex')
xlabel('iterations $$i$$','interpreter','latex')

%% Compare CGS- and RGS-FGMRES-(M)DR

% CGS
[x_rl,res_rl,iter_rl]=FGMRES_dr_ras(A,b,tol,40,15,0,6,12);
[x_drl,res_drl,iter_drl]=FGMRES_dr_ras(A,b,tol,40,15,10,6,12);
[x_mdrl,res_mdrl,iter_mdrl]=FGMRES_mdr_ras(A,b,tol,40,15,10,6,12);

% RGS
[x_Rrl,res_Rrl,iter_Rrl]=RGS_FGMRES_dr_ras(A,b,tol,40,15,0,1000,6,12);
[x_Rdrl,res_Rdrl,iter_Rdrl]=RGS_FGMRES_dr_ras(A,b,tol,40,15,10,1000,6,12);
[x_Rmdrl,res_Rmdrl,iter_Rmdrl]=RGS_FGMRES_mdr_ras(A,b,tol,40,15,10,1000,6,12);

figure(3)
semilogy(iter_rl,res_rl,'-',iter_drl,res_drl,':',...
    iter_mdrl,res_mdrl,'-.',iter_Rrl,res_Rrl,'-x',...
    iter_Rdrl,res_Rdrl,'-s',iter_Rmdrl,res_Rmdrl,'-o')
% 
legend('FGMRES(40,15)','FGMRES DR(40,15,10)','FGMRES MDR(40,15,10)',...
'RGS FGMRES(40,15)','RGS FGMRES DR(40,15,10)','RGS FGMRES MDR(40,15,10)',...
'Location','southwest')
xlim([0,200])
xlabel('Iterations')
ylabel('Relative residual norms')