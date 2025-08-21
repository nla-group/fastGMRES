function [Q,R,iter,n_cond,n_condS]=RGS_cond(W,k)
[n,m]=size(W);

m2=floor(m/10);
iter=zeros(m2,1);
n_cond=zeros(m2,1);
n_condS=zeros(m2,1);

% define the Rademacher matrix for sketch
theta=Rademacher_r(k,n);

Q=zeros(n,m);
R=zeros(m,m);
S=zeros(k,m);

% orthogonalization with RGS
p=theta*W(:,1);
R(1,1)=norm(p,2);
S(:,1)=p/R(1,1);
Q(:,1)=W(:,1)/R(1,1);

for i=2:m
    p=theta*W(:,i);
    R(1:i-1,i)=S(:,1:i-1)\p;
    q=W(:,i)-Q(:,1:i-1)*R(1:i-1,i);
    s=p-S(:,1:i-1)*R(1:i-1,i);
    R(i,i)=norm(s,2);
    S(:,i)=s/R(i,i);
    Q(:,i)=q/R(i,i);
end

% compute condition numbers
for it=1:m2
    iter(it)=10*it;
    n_cond(it)=cond(Q(:,1:10*it));
    n_condS(it)=cond(S(:,1:10*it));
end

end