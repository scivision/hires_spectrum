function T=Tmean(That)
% function T=Tmean_centralized(That)
% compute the transportation mean of covariance matrices by solving LMI
% input:  That --- covariance matrices with size nxnxm
% output: T    --- Toeplitz structured covariance

[n,n1,m]=size(That);
I=eye(n);
A=[I -I;-I I];
cvx_begin quiet
variable S(n,n,m)
variable t(1,n)
T=toeplitz(t);
obj=0;
for k=1:m
    [T S(:,:,k);S(:,:,k)' That(:,:,k)]==semidefinite(2*n,2*n);
    obj=obj+trace(A*[T S(:,:,k);S(:,:,k)' That(:,:,k)]);
end
minimize(obj)
cvx_end
