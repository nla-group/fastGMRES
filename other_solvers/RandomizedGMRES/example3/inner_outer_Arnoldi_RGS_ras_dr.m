function [Q,R,Z,S]=inner_outer_Arnoldi_RGS_ras_dr(A,b,m,m_inner,k,theta,...
                                                 n_lap,n_subd,Q_in,R_in,...
                                                 Z_in,S_in)
n=length(b);
Q=zeros(n,m+1);
R=zeros(m+1,m);
R(1:k+1,1:k)=R_in(1:k+1,1:k);
Z=zeros(n,m);

Z(:,1:k)=Z_in(:,1:k);
Q(:,1:k+1)=Q_in(:,1:k+1);
S(:,1:k+1)=S_in(:,1:k+1);

for i=k+1:m
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