load('LS89.mat');
n=length(A);
b=load('RHS.dat');
tol=1e-10;
t=1000;
theta=Rademacher_r(t,n);
m=100;
minner=15;

% m is outer FGMRES dimension
% m_inner is 
%RGS_FGMRES_mdr_ras(A,b,tol,m,m_inner,k,t,n_lap,n_subd)
[x_Rmdrl,res_Rmdrl,iter_Rmdrl]=RGS_FGMRES_mdr_ras(A,b,tol,40,15,10,1000,6,12);