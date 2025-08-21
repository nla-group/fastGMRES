function [V,H,Z]=MGS_ras(A,r,m,m_inner,n_lap,n_sub)
n=length(r);
V=zeros(n,m+1);
H=zeros(m+1,m);
Z=zeros(n,m);
V(:,1)=r/norm(r);

for j=1:m
    Z(:,j)=inner_GMRES_ras(A,V(:,j),m_inner,n_lap,n_sub);
    w=A*Z(:,j);

    for i=1:j
        H(i,j)=w'*V(:,i);
        w=w-H(i,j)*V(:,i);       
    end

    H(j+1,j)=norm(w);

    if H(j+1,j)==0
        V(:,j+1)=0;
    else
        V(:,j+1)=w/H(j+1,j);
    end
end