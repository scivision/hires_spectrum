function T=CovEst_kl(hatT)
if(norm(hatT-hatT')>.001||min(eig(hatT)<0))
    error('Sample covariance must be PSD');
end
n=length(hatT);

cvx_begin
variable t(1,n)
T=toeplitz(t);
minimize(-log_det(T)+trace(T/hatT))
T==semidefinite(n,n);
cvx_end
T=full(T);
