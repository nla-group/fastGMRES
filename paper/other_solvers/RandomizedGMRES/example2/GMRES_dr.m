function [x,residual,iter]=GMRES_dr(A,b,tol,m,k,M,M2)
n=length(b);
r0=b;
beta=norm(r0);
nb=beta;
c=[beta,zeros(1,m)]';

x=zeros(n,1);
iterMax=500;
iter=zeros(iterMax,1);
residual=zeros(iterMax,1);
V=zeros(n,m+1);
H=zeros(m+1,m);
Z=zeros(n,m);
V(:,1)=r0/beta;
residual(1)=1;
iter(1)=0;

% first m-Arnoldi iterations
for j=1:m
    if isempty(M)==1
        Z(:,j)=V(:,j);
    else
    Z(:,j)=M\V(:,j);
        if isempty(M2)==0
            Z(:,j)=M2\Z(:,j);
        end
    end
    w=A*Z(:,j);
    
    for i=1:j
        H(i,j)=V(:,i)'*w;
    end
    for i=1:j
        w=w-H(i,j)*V(:,i);
    end
    
    H(j+1,j)=norm(w);
    
    if H(j+1,j)==0
        V(:,j+1)=0;
    else
        V(:,j+1)=w/H(j+1,j);
    end
end

y=H\c;
x=x+Z*y;
r0=b-A*x;

beta=norm(r0);
residual(2)=beta/nb;
iter(2)=m;
it=1;
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
        e=diag(D);
        e=abs(e);
        [~,idx]=sort(e,'ascend');
        id=idx(1:k);
        G=E(:,id);
        [P,~]=qr(G,0);
        P2=[P;zeros(1,k)];
        AP=[P2,rho];
        
        [QQ,~]=qr(AP,0);
    
        V_n=V*QQ;
        Z_n=Z*QQ(1:m,1:k);
        H_n=QQ'*H*QQ(1:m,1:k);
        V(:,1:k+1)=V_n;
        Z(:,1:k)=Z_n;
        H(1:k+1,1:k)=H_n;
        
    end
    
    % (m-k) Arnoldi iterations
    if k==0
        V(:,1)=r0/beta;
    end
    for j=k+1:m
        if isempty(M)==1
            Z(:,j)=V(:,j);
        else
            Z(:,j)=M\V(:,j);
            if isempty(M2)==0
                Z(:,j)=M2\Z(:,j);
            end
        end
        w=A*Z(:,j);
        for i=1:j
            H(i,j)=V(:,i)'*w;
        end
        for i=1:j
            w=w-H(i,j)*V(:,i);
        end
        H(j+1,j)=norm(w);
        if H(j+1,j)==0
            V(:,j+1)=0;
            fprintf('A zero basis vector appears at %d-th\n',j+1)
    
        else
            V(:,j+1)=w/H(j+1,j);
        end
    end
    
    if k==0
       c=[beta,zeros(1,m)]';
    else
        c=[QQ'*rho;zeros(m-k,1)];
    end
    y=H\c;
    x=x+Z*y;
    r0=b-A*x;    
    beta=norm(r0);   
    residual(it+1)=norm(r0)/nb;
    if it==1
        iter(it+1)=iter(it)+(m-k)+m;
    else
        iter(it+1)=iter(it)+(m-k);
    end
end
residual=residual(1:it+1);
iter=iter(1:it+1);
end