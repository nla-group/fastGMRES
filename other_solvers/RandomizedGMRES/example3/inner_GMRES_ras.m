function x=inner_GMRES_ras(A,b,m,n_lap,n_sub)
n=length(b);
r0=b;
beta=norm(r0);
c=[beta,zeros(1,m)]';

x=zeros(n,1);
V=zeros(n,m+1);
H=zeros(m+1,m);
Z=zeros(n,m);
V(:,1)=r0/beta;
if m==0
    x=prec_ras(A,V,n_lap,n_sub);
else
    for j=1:m
        Z(:,j)=prec_ras(A,V(:,j),n_lap,n_sub);
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
    y=H\c;
    x=x+Z*y;
end