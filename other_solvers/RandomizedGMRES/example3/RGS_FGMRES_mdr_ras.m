function [x,residual,iter]=RGS_FGMRES_mdr_ras(A,b,tol,m,m_inner,k,t,n_lap,n_subd)
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
ek=zeros(m+1,1);
ek(k+1)=1;

% first m-Arnoldi iterations
c=[norm(theta*r0),zeros(1,m)]';
[V,H,Z,S]=inner_outer_Arnoldi_RGS_ras(A,r0,m,m_inner,theta,n_lap,n_subd);

loss=norm(eye(m+1)-S'*S);
fprintf('======================================\n')
fprintf('First cycle\n');
fprintf('Loss of orthogonality of Sm+1 = %.3e\n',loss);

y=H\c;
x=x+Z*y;
r0=b-A*x;
W=S(:,1:m);
beta=norm(r0);
residual(2)=beta/nb;
iter(2)=m;
it=1;
if beta<rtol
    residual=residual(1:it+1);
    iter=iter(1:it+1);
    return
end

while beta>rtol&&iter(it+1)<iterMax
    it=it+1;
    fprintf('%d-cycle cycle\n',it);

     % new computation with singular vector forms
    if k~=0       
        T1=H'*H;
    if it==2
        T2=eye(m,m);
    else
        T2(1:k,1:k)=W_n'*W_n;
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
    loss_k=norm(eye(k)-S_n'*S_n);
    fprintf('Loss of orthogonality of Sk = %.3e\n',loss_k);    
    end
    
    % (m-k) Arnoldi iterations
    r0=r0-V_n*(S_n'*theta*r0);
    [V2,H2,Z2,S2]=inner_outer_Arnoldi_RGS_ras_mdr(A,r0,m-k,m_inner,theta,...
        n_lap,n_subd,V_n,S_n);
    
    H=[eye(k,k),(S_n')*theta*A*Z2;zeros(m-k+1,k),H2];
    V=[V_n,V2];
    Z=[Z_n,Z2];
    S=[S_n,S2];    
    W=[W_n,S2(:,1:end-1)];

    loss=norm(eye(m+1)-S'*S);
    fprintf('Loss of orthogonality of Sm+1 = %.3e\n',loss);
    
    if k==0
        c=[norm(theta*r0),zeros(1,m)]';
    else
        c=norm(theta*r0)*ek;
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