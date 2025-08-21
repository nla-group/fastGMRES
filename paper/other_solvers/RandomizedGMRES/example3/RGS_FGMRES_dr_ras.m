function [x,residual,iter]=RGS_FGMRES_dr_ras(A,b,tol,m,m_inner,k,t,...
                                             n_lap,n_subd)
n=length(b);
r0=b;
beta=norm(r0);
nb=beta;
rtol=tol*nb;
x=zeros(n,1);
theta=Rademacher_r(t,n);
iterMax=200;
iter=zeros(iterMax,1);
residual=zeros(iterMax,1);
residual(1)=1;
iter(1)=0;

c=[norm(theta*r0),zeros(1,m)]';

% first m-Arnoldi iterations
[V,H,Z,S]=inner_outer_Arnoldi_RGS_ras(A,r0,m,m_inner,theta,n_lap,n_subd);
loss=norm(eye(m+1)-S'*S);

fprintf('======================================\n')
fprintf('First cycle\n');
fprintf('Loss of orthogonality of Sm+1 = %.3e\n',loss);
y=H\c;
x=x+Z*y;
r0=b-A*x;
beta=norm(r0);
residual(2)=beta/nb;
iter(2)=m;
it=1;
if beta<rtol
    residual=residual(1:it+1);
    iter=iter(1:it+1);
    return
end
em=zeros(m,1);
em(m)=1;
while beta>rtol&&iter(it+1)<iterMax
    it=it+1;
    fprintf('%d-cycle cycle\n',it);
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
        S_n=S*QQ;
        V(:,1:k+1)=V_n;
        H(1:k+1,1:k)=H_n;
        Z(:,1:k)=Z_n;
        S(:,1:k+1)=S_n;

        loss_k=norm(eye(k+1)-S_n'*S_n);
        fprintf('Loss of orthogonality of Sk+1 = %.3e\n',loss_k);

        % (m-k) Arnoldi iterations
        [V,H,Z,S]=inner_outer_Arnoldi_RGS_ras_dr(A,r0,m,m_inner,k,...
            theta,n_lap,n_subd,V,H,Z,S);
        
        loss=norm(eye(m+1)-S'*S);
        fprintf('Loss of orthogonality of Sm+1 = %.3e\n',loss);
    else
        [V,H,Z,S]=inner_outer_Arnoldi_RGS_ras(A,r0,m,m_inner,theta,...
            n_lap,n_subd);

        loss=norm(eye(m+1)-S'*S);
        fprintf('Loss of orthogonality of Sm+1 = %.3e\n',loss);
    end
    
    if k==0
        c=[norm(theta*r0),zeros(1,m)]';
    else
        c=[QQ'*rho;zeros(m-k,1)];
    end

    y=H\c;
    x=x+Z*y;
    r0=b-A*x;    
    beta=norm(r0);   
    residual(it+1)=beta/nb; 
    iter(it+1)=iter(it)+(m-k);
end
residual=residual(1:it+1);
iter=iter(1:it+1);
end