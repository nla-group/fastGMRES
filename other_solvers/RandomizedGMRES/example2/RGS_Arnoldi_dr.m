function [Q,R,Z,S]=RGS_Arnoldi_dr(A,theta,m,k,M1,M2,Q_in,R_in,Z_in,S_in)
n=length(A);
Q=zeros(n,m+1);
R=R_in;
Z=zeros(n,m);

Z(:,1:k)=Z_in(:,1:k);
Q(:,1:k+1)=Q_in(:,1:k+1);
S(:,1:k+1)=S_in(:,1:k+1);


for i=k+1:m
    if isempty(M1)==0
        Z(:,i)=M1\Q(:,i);
        if isempty(M2)==0
            Z(:,i)=M2\Z(:,i);
        end
    else
       Z(:,i)=Q(:,i); 
    end
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