function [x,residual,iter,timevec]=RGS_GMRES_mdr(A,b,tol,m,k,t,M1,M2)
n=length(b);
r0=b;
nb=norm(b);
x=zeros(n,1);

%theta=Rademacher_r(t,n);
%theta=@(x) theta*x;

%theta = srft(n, t); % SG
theta = clarkson_woodruff(n, t);

iterMax=500; % SG: this is the total number of matvecs
residual=zeros(iterMax,1);
iter=zeros(iterMax,1);
residual(1)=1;
ek=zeros(m+1,1);
ek(k+1)=1;

timevec = [0];
tstart = tic;

% first m Arnoldi iterations
fprintf('RGS_GMRES_mdr iter = 1\n')
[V,H,Z,S]=RGS_Arnoldi(A,r0,theta,m,M1,M2);
c=[norm(theta(r0)),zeros(1,m)]';
y=H\c;
x=x+Z*y;
r0=b-A(x);
beta=norm(r0);
W=S(:,1:m);
residual(2)=beta;
fprintf('  residual norm = %e\n', residual(2))
timevec(2) = toc(tstart);
iter(2)=m;
it=1;


if beta<tol*nb
    residual=residual(1:it+1);
    iter=iter(1:it+1);
    return
end
    
%while beta>tol*nb&&iter(it+1)<iterMax
while beta>tol*nb&&it<800
    it=it+1;
    fprintf('RGS_GMRES_mdr iter = %d\n', it)

    % compute singular vectors
    T1=H'*H;
    if it==2
        T2=eye(m,m);
    else
        T2=[W_n'*W_n,zeros(k,m-k);zeros(m-k,k),eye(m-k,m-k)];
    end
    [E,EV]=eig(T1,T2);
    e=diag(EV);
    [~,idx]=sort(e,'ascend');
    id=idx(1:k);
    G=E(:,id);
    [Q,R]=qr(H*G,0);
    invR=pinv(R);
    V_n=V*Q;
    S_n=S*Q;
    Z_n=Z*G*invR;
    W_n=W*G*invR;

    % (m-k) Arnoldi iterations
    [V2,H2,Z2,S2]=RGS_Arnoldi_mdr(A,r0,theta,m-k,M1,M2,V_n,S_n);
    V=[V_n,V2];
    S=[S_n,S2];
    Z=[Z_n,Z2];
    W=[W_n,S2(:,1:end-1)];
    H=[eye(k,k),(S_n')*theta(A(Z2));zeros(m-k+1,k),H2];

    c=norm(theta(r0))*ek;
    y=H\c;    
    x=x+Z*y;   
    r0=b-A(x);    
    beta=norm(r0);
    residual(it+1)=beta;
    fprintf('  residual norm = %e\n', residual(it+1))
    timevec(it+1) = toc(tstart);
    if it==1
        iter(it+1)=iter(it)+(m-k)+m;
    else
        iter(it+1)=iter(it)+(m-k);
    end
    
end
residual=residual(1:it+1);
iter=iter(1:it+1);