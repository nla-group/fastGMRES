function [x,residual,iter]=GMRES_mdr(A,b,tol,m,k,M,M2)
n=length(b);
r0=b;
beta=norm(r0);
nb=beta;
c=[beta,zeros(1,m)]';

x=zeros(n,1);
iterMax=500;
iter=zeros(iterMax,1);
residual=zeros(iterMax,1);
V2=zeros(n,m-k+1);
V=zeros(n,m+1);
H=zeros(m+1,m);
H2=zeros(m-k+1,m-k);
Z=zeros(n,m);
Z2=zeros(n,m-k);
V(:,1)=r0/beta;
residual(1)=1;
iter(1)=0;
ek=zeros(m+1,1);
ek(k+1)=1;
V_n=0;
W_n=zeros(k,k);

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
W=V(:,1:m);

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


while beta>tol*nb&&iter(it+1)<iterMax
    it=it+1;
    
    % new computation for singular vectors
    if k~=0
        T1=H'*H;
        if it==1
        T2=eye(m,m);
        else
        T2=[W_n'*W_n,zeros(k,m-k);zeros(m-k,k),eye(m-k,m-k)];
        end
        [E,D]=eig(T1,T2);
        e=diag(D);
        [~,idx]=sort(e);
        id=idx(1:k);
        P=E(:,id);
        [Q,R]=qr(H*P,0);
        invR=pinv(R);
        V_n=V*Q;
        Z_n=Z*P*invR;
        W_n=W*P*invR;
    end

    % (m-k) Arnoldi iterations
    V2(:,1)=r0/norm(r0);
    
    for j=1:m-k
        Z2(:,j)=V2(:,j);
        if isempty(M)==0
            Z2(:,j)=M\Z2(:,j);
            if isempty(M2)==0
                Z2(:,j)=M2\Z2(:,j);
            end
        end
        w=A*Z2(:,j);
        w=w-V_n*(V_n'*w);
        for i=1:j
            H2(i,j)=V2(:,i)'*w;
        end
        for i=1:j
            w=w-H2(i,j)*V2(:,i);
        end

        H2(j+1,j)=norm(w);
        if H2(j+1,j)==0
            V2(:,j+1)=0;
            fprintf('A zero basis vector appears at %d-th\n',j+1)    
        else
            V2(:,j+1)=w/H2(j+1,j);
        end
    end
    if k>1
    H=[eye(k,k),V_n'*A*Z2;zeros(m-k+1,k),H2];
    V=[V_n,V2];
    Z=[Z_n,Z2];
    W=[W_n,V2(:,1:end-1)];
    else
        H=H2;
        V=V2;
        Z=Z2;
    end
    if k==0
       c=[beta,zeros(1,m)]';
    else
       c=beta*ek;
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