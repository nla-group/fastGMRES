clear 
clc 
%load('naca12A.mat');
load vas_stokes_1M.mat
A = Problem.A;
n=length(A);
%b=load('RHS.dat');
b = randn(n,1); b = b/norm(b);
[L,U]=ilu(A);
tol=1e-16;

%[x,residual,iter]=RGS_GMRES_mdr(A,b,tol,m,k,t,M1,M2)
% t sketching dimension
% k deflated harmonic Ritz vectors between FGMRES restarts
% m FGMRES basis to be computed
% there were at most 500 total iterations (hard coded), 
% changed to 20 outer FGMRES restarts
% we have also changed every occurence of A*x to A(x), 
% so that the code can work with function handles, and provide the
% preconditioner in form of a modified matrix to be consistent with 
% our other codes. 

A = @(x) U\(L\(A*x));
L = [];
U = [];
[x,Rres_mdr40,Riter_mdr40,timevec]=RGS_GMRES_mdr(A,b,tol,120,20,1000,L,U);

