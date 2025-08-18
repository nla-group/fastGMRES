clear 
clc 
load('naca12A.mat');
n=length(A);
b=load('RHS.dat');
[L,U]=ilu(A);
tol=1e-8;


%% CGS 
% GMRES(40)
[~,res_40,iter_40]=GMRES_dr(A,b,tol,40,0,L,U);

% GMRES(40,20)
[~,res_dr40,iter_dr40]=GMRES_dr(A,b,tol,40,20,L,U);

% GMRES_SVD(40,20)
[~,res_mdr40,iter_mdr40]=GMRES_mdr(A,b,tol,40,20,L,U);


%% RGS
% GMRES(40)
[~,Rres_40,Riter_40]=RGS_GMRES_dr(A,b,tol,40,0,1000,L,U);

% GMRES(40,20)
[~,Rres_dr40,Riter_dr40]=RGS_GMRES_dr(A,b,tol,40,20,1000,L,U);

% GMRES_SVD(40,20)
[~,Rres_mdr40,Riter_mdr40]=RGS_GMRES_mdr(A,b,tol,40,20,1000,L,U);

%% Figure
figure(1)
semilogy(iter_40,res_40,'-',iter_dr40,res_dr40,':',...
    iter_mdr40,res_mdr40,'--',Riter_40,Rres_40,'-*',...
    Riter_dr40,Rres_dr40,'-s',Riter_mdr40,Rres_mdr40,'-o')
legend('GMRES(40)','GMRES DR(40,20)','GMRES MDR(40,20)',...
'RGS GMRES(40)','RGS GMRES DR(40,20)','RGS GMRES MDR(40,20)')
xlabel('Iterations')
ylabel('Relative residual norms')
xlim([0,400])

%% Various deflation numbers

% RGS FGMRES(60,20)
[~,Rres_dr20,Riter_dr20]=RGS_GMRES_dr(A,b,tol,60,20,1000,L,U);
% RGS FGMRES(60,30)
[~,Rres_dr30,Riter_dr30]=RGS_GMRES_dr(A,b,tol,60,30,1000,L,U);
% RGS FGMRES(60,30)
[~,Rres_dr40,Riter_dr40]=RGS_GMRES_dr(A,b,tol,60,40,1000,L,U);



% RGS FGMRES MDR(60,20)
[~,Rres_mdr20,Riter_mdr20]=RGS_GMRES_mdr(A,b,tol,60,20,1000,L,U);
% RGS FGMRES MDR(60,30)
[~,Rres_mdr30,Riter_mdr30]=RGS_GMRES_mdr(A,b,tol,60,30,1000,L,U);
% RGS FGMRES MDR(60,40)
[~,Rres_mdr40,Riter_mdr40]=RGS_GMRES_mdr(A,b,tol,60,40,1000,L,U);

%%
figure(2)
semilogy(Riter_dr20,Rres_dr20,':',...
    Riter_dr30,Rres_dr30,'--',Riter_dr40,Rres_dr40,'-*',...
    Riter_mdr20,Rres_mdr20,'-x',Riter_mdr30,Rres_mdr30,'-d',...
    Riter_mdr40,Rres_mdr40,'-s')

legend('DR(60,20)','DR(60,30)','DR(60,40)',...
    'MDR(60,20)','MDR(60,30)','MDR(60,40)')
xlabel('Iterations')
ylabel('Relative residual norms')
xlim([0,220])

