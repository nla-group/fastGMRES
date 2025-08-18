function S = clarkson_woodruff(n, s)

i = randi(s,n,1);
j = 1:n;
el = 2*randi(2,n,1)-3;
S = sparse(i,j,el,s,n,n);
S = @(x) S*x;