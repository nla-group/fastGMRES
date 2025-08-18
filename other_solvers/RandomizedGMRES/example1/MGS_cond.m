function [Q,R,iter,n_cond]=MGS_cond(W)
[n,m]=size(W);

m2=floor(m/10);
iter=zeros(m2,1);
n_cond=zeros(m2,1);

Q=zeros(n,m);
R=zeros(m,m);

% orthogonalization with MGS
for j=1:m
    v=W(:,j); 
    for i=1:j-1
        R(i,j)=Q(:,i)'*v;
        v=v-R(i,j)*Q(:,i);
    end
    R(j,j)=norm(v);
    Q(:,j)=v/R(j,j);
end


% compute condition numbers
for it=1:m2
    iter(it)=10*it;
    n_cond(it)=cond(Q(:,1:10*it));
end

end