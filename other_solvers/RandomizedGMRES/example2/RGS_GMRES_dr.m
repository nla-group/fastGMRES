function [x,residual,iter]=RGS_GMRES_dr(A,b,tol,m,k,t,M1,M2)
n=length(b);
r0=b;
nb=norm(b);
x=zeros(n,1);
theta=Rademacher_r(t,n);

iterMax=500;
residual=zeros(iterMax,1);
iter=zeros(iterMax,1);
residual(1)=1;
    
% first m Arnoldi iterations
[V,H,Z,S]=RGS_Arnoldi(A,r0,theta,m,M1,M2);
c=[norm(theta*r0),zeros(1,m)]';
y=H\c;
x=x+Z*y;
r0=b-A*x;
beta=norm(r0);
iter(2)=m;
it=1;
residual(2)=beta/nb;

if beta<tol*nb
    residual=residual(1:it+1);
    iter=iter(1:it+1);
    return
end
em=zeros(m,1);
em(m)=1;
    
while beta>tol*nb&&iter(it+1)<iterMax    
    it=it+1;
    % new computation with harmonic Ritz forms
    if k~=0
        rho=c-H*y;
        H2=H(1:m,1:m);
        h=H(m+1,m);
        f=pinv(H2')*em;
        T=H2+(h*h*f*(em'));
        [E,D]=eig(T);
	    e=diag(abs(D));
	    [~,idx]=sort(e,'ascend');
%         fprintf('Eigenvalues\n')
%         disp(e(idx(1:k)))
        G=E(:,idx(1:k));
        P2=[G;zeros(1,k)];
        AP=[P2,rho];
        [QQ,~]=qr(AP,0);
    
        V_n=V*QQ;
        Z_n=Z*QQ(1:m,1:k);
        H_n=QQ'*H*QQ(1:m,1:k);
        S_n=S*QQ;
        V(:,1:k+1)=V_n;
        Z(:,1:k)=Z_n;
        H(1:k+1,1:k)=H_n;
        S(:,1:k+1)=S_n;
    end
    

    % (m-k) Arnoldi iterations
    if k==0
        [V,H,Z,~]=RGS_Arnoldi(A,r0,theta,m,M1,M2);
        c=[norm(theta*r0),zeros(1,m)]';
        
    else
        [V,H,Z,S]=RGS_Arnoldi_dr(A,theta,m,k,M1,M2,V,H,Z,S);        
        c=[QQ'*rho;zeros(m-k,1)];
    end
    y=H\c;
    x=x+Z*y;
    r0=b-A*x;
    beta=norm(r0);
    residual(it+1)=beta/nb;
    
    if it==1
        iter(it+1)=iter(it)+(m-k)+m;
    else
        iter(it+1)=iter(it)+(m-k);
    end
    
end

residual=residual(1:it+1);
iter=iter(1:it+1);
