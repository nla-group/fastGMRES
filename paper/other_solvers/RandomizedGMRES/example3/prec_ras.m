function b_out=prec_ras(A_in,b_in,n_lap,n_subd)
n=length(b_in);
b_out=zeros(n,1);

n0=round(n/n_subd);
% n1=n-n0*(n_subd-1);

i0=1;
i1=n0;

% For the first subdomain
% compute Ri*b_in
temp=b_in(i0:i1+n_lap,1);

% compute Ri*A*Ri'
subA=A_in(i0:i1+n_lap,i0:i1+n_lap);

% compute inv(Ri*A*Ri')*Ri*b_in
temp=subA\temp;

% multiply Di
temp(end-2*n_lap+1:end)=0.5*temp(end-2*n_lap+1:end);

% multiply Ri'
b_out(i0:i1+n_lap)=temp;

for i=1:n_subd-2

    i0=i*n0+1;
    i1=(i+1)*n0;
    
    % compute Ri*b_in
    temp=b_in(i0-n_lap:i1+n_lap,1);
    
    % compute Ri*A*Ri'
    subA=A_in(i0-n_lap:i1+n_lap,i0-n_lap:i1+n_lap);
    
    % compute inv(Ri*A*Ri')*Ri*b_in
    temp=subA\temp;

    % multiply Di
    temp(1:2*n_lap)=0.5*temp(1:2*n_lap);
    temp(end-2*n_lap+1:end)=0.5*temp(end-2*n_lap+1:end);
    b_out(i0-n_lap:i1+n_lap,1)=b_out(i0-n_lap:i1+n_lap,1)+temp;
    
end

i0=i1+1;
temp=b_in(i0-n_lap:end,1);
subA=A_in(i0-n_lap:end,i0-n_lap:end);
temp=subA\temp;
temp(1:2*n_lap)=0.5*temp(1:2*n_lap);
b_out(i0-n_lap:end,1)=b_out(i0-n_lap:end,1)+temp;
