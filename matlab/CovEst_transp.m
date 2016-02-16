function T=CovEst_transp(hatT,A,B)
% function T=CovEst_LMI(hatT,A,B)
% state covariance (Toeplitz covariance) matrix estimation via transportation/Hellinger
% distance minimization

if(norm(hatT-hatT')>.001||min(eig(hatT)<0))
    error('Sample covariance must be PSD');
end
n=length(hatT);

if(nargin==1) %Toeplitz matrix estimation
    
    cvx_begin quiet
    variable t(1,n)
    T=toeplitz(t);
    variable S(n,n);
    minimize(trace(T+hatT-S-S'))
    [T S;S' hatT]==semidefinite(2*n,2*n);
    cvx_end
    T=full(T);
else if(nargin==3) % state co
    [mb,nb]=size(B);    
    cvx_begin quiet
    variable T(n,n)
    variable S(n,n);
    variance H(nb,mb);
    minimize(trace(T+hatT-S-S'))
    [T S;S' hatT]==semidefinite(2*n,2*n);
    T-A*T*A'==B*H+H'*B';
    cvx_end
    T=full(T);
    end
end

