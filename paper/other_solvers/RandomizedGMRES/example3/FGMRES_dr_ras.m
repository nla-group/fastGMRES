function [x,residual,iter]=FGMRES_dr_ras(A,b,tol,m,m_inner,k,n_lap,n_sub)
n=length(b);
r0=b;
beta=norm(r0);
nb=beta;
c=[beta,zeros(1,m)]';

x=zeros(n,1);
iterMax=200;
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
    Z(:,j)=inner_GMRES_ras(A,V(:,j),m_inner,n_lap,n_sub);
    w=A*Z(:,j);

    for i=1:j
        H(i,j)=w'*V(:,i);       
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
loss=norm(eye(m+1)-V'*V);

fprintf('======================================\n')
fprintf('Loss of orthogonality in the first cycle = %.3e\n',loss);

y=H\c;
x=x+Z*y;
r0=b-A*x;
beta=norm(r0);
residual(2)=beta/nb;
iter(2)=m;
it=1;
if beta<nb*tol
    residual=residual(1:it+1);
    iter=iter(1:it+1);
    return
end

em=zeros(m,1);
em(m)=1;
  
    
while beta>tol*nb&&iter(it+1)<iterMax
    it=it+1;
    
    fprintf('%d -th cycle\n',it);
    % new computation with harmonic Ritz forms
    if k~=0
        rho=c-H*y;
        H2=H(1:m,1:m);
        h=H(m+1,m);
        f=pinv(H2')*em;
        T=H2+(h*h*f*(em'));
        [E,D]=eig(T);
        e=diag(D);
        [~,idx]=sort(abs(e),"ascend");
        id=idx(1:k);
        G=E(:,id);
        P2=[G;zeros(1,k)];
        AP=[P2,rho];
        
        [QQ,~]=qr(AP,0);
    
        V_n=V*QQ;
        Z_n=Z*QQ(1:m,1:k);
        H_n=QQ'*H*QQ(1:m,1:k);
        V(:,1:k+1)=V_n;
        Z(:,1:k)=Z_n;
        H(1:k+1,1:k)=H_n;
        loss_k=norm(eye(k+1)-V_n'*V_n);
        fprintf('Loss of orthogonality of Vk+1 = %.3e\n',loss_k);
    end
    
    % (m-k) Arnoldi iterations
    if k==0
        V(:,1)=r0/beta;
    end
    for j=k+1:m
        Z(:,j)=inner_GMRES_ras(A,V(:,j),m_inner,n_lap,n_sub);
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
    loss=norm(eye(m+1)-V'*V);
    fprintf('Loss of orthogonality of Vm+1 = %.3e\n',loss);
    if k==0
       c=[beta,zeros(1,m)]';
    else
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
end