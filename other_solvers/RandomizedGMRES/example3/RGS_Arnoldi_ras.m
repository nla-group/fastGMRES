function [Q,R,Z]=RGS_Arnoldi_ras(A,b,theta,m,n_lap,n_subd)
n=length(b);
Q=zeros(n,m+1);
Z=zeros(n,m);
R=zeros(m+1,m);


w=b;
p=theta*w;
newq=w;
news=p;
h=norm(news,2);
S(:,1)=news/h;
Q(:,1)=newq/h;

for i=1:m    
    Z(:,i)=prec_ras(A,Q(:,i),n_lap,n_subd);
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