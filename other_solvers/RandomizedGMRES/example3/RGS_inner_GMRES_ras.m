function x=RGS_inner_GMRES_ras(A,b,m,theta,n_lap,n_subd)
n=length(b);
r0=b;
c=[norm(theta*r0),zeros(1,m)]';

x=zeros(n,1);
if m>1
    [~,H,Z]=RGS_Arnoldi_ras(A,r0,theta,m,n_lap,n_subd);
    y=H\c;
    x=x+Z*y;
else
    x=prec_ras(A,r0,n_lap,n_subd);
end