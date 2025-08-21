function [x,res,sres,cnd,scnd] = sgmres(A, b, m, cndtol, t)
% plain sgmres implementation that breaks when cond(SAV)>1e15 (by default)
% optional condition number tolerance cndtol and truncation parameter t

N = size(b,1);
s = 2*m; 
%hS = srft(N, s); 
hS = clarkson_woodruff(N,s);

if nargin < 4
    cndtol = 1e15; 
end

if nargin < 5
    t = 2; % t = 0: monomial; t = 2: Lanczos-like
end

res = norm(b); sres = []; cnd = []; scnd = [];
H = zeros(m+1,m); V = zeros(N,m+1); 
v = b/norm(b); V(:,1) = v; 
SAV = zeros(s,m+1); 
SV = zeros(s,m+1); 
SV(:,1) = hS(V(:,1));
Sb = hS(b);
fprintf('  running sgmres ')
for j = 1:m

    fprintf('.')
    %cnd(j) = cond(V(:,1:j));
    %scnd(j) = cond(SV(:,1:j));
    
    w = A(v);
    SAV(:,j) = hS(w);

    [Q,R] = qr(SAV(:,1:j),0);

    if j > 1 && cond(R) > cndtol
        break
    end

    y = R\(Q'*Sb);   % sketched GMRES
    %sres(j+1) = norm(hS(b - A*x));
    %sres(j+1) = norm(Sb - SAV(:,1:j)*y);

    for i = max(j-t+1,1):j
        H(i,j) = V(:,i)'*w;
        w = w - V(:,i)*H(i,j);
    end
    H(j+1,j) = norm(w);
    v = w/H(j+1,j);
    V(:,j+1) = v;

end
fprintf('\n')

x = V(:,1:length(y))*y;
