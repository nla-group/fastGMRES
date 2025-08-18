function [Q,R,Z,S]=inner_outer_Arnoldi_RGS_ras(A,b,m,m_inner,theta,n_lap,n_subd)
n=length(b);
Q=zeros(n,m+1);
R=zeros(m+1,m);
Z=zeros(n,m);

w=b;

p=theta*w;
newq=w;
news=p;
h=norm(news,2);
S(:,1)=news/h;
Q(:,1)=newq/h;

for i=1:m
    Z(:,i)=RGS_inner_GMRES_ras(A,Q(:,i),m_inner,theta,n_lap,n_subd); 
    w=A*Z(:,i);
    p=theta*w;
    R(1:i,i)=S(:,1:i)\p;
    newq=w-Q(:,1:i)*R(1:i,i);
    news=theta*newq;
    R(i+1,i)=norm(news,2);
    S(:,i+1)=news/R(i+1,i);
    Q(:,i+1)=newq/R(i+1,i);
end

end