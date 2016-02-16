function [Ts,Tma] = MA_decomp(T,m)
% [Tsignal,Tnoise]=MA_decomposition(T,m)
%
%  T = Tsignal + Tnoise
%
% where Tnoise corresponding to a MA(m) process
% with spectrum Q0+2*Q1*cos(theta)+...+2*Qm*cos(m*theta)
% i.e., with m+1 coefficients.
%
% T,Tsignal,Tnoise are all Toeplitz covariance matrices.
% The routine uses convex optimization tool cvx to solve an LMI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = length(T);
A = compan(eye(1,m+1));
B = eye(m,1);
if(m>0)
    if isreal(T) % real valued covariance matrix
        cvx_begin quiet
        variable c(1,m+1);
        variable P(m,m) symmetric;
        C=c(2:m+1);
        minimize(-c(1))
        P==semidefinite(m,m);
        Tma=toeplitz([c zeros(1,n-m-1)]);
        T-Tma==semidefinite(n,n);
        [P-A'*P*A  C'-A'*P*B; C-B'*P*A c(1)-B'*P*B]==semidefinite(m+1,m+1);
        cvx_end
        Tma=full(Tma);
        Ts=T-Tma;
    else
        cvx_begin quiet
        variable c0(1,1);
        variable Creal(1,m);
        variable Cimag(1,m);
        variable Preal(m,m) symmetric;
        variable Pimag(m,m);
        Tma_real = toeplitz([c0 Creal zeros(1,n-m-1)]);
        Tma_imag = toeplitz([0  Cimag zeros(1,n-m-1)],[0  -Cimag zeros(1,n-m-1)]);

        Mreal = [Preal-A'*Preal*A Creal'-A'*Preal*B; 
                 Creal-B'*Preal*A    c0-B'*Preal*B];
        Mimag = [Pimag-A'*Pimag*A  -Cimag'-A'*Pimag*B;
                 Cimag-B'*Pimag*A  -B'*Pimag*B];
        minimize(-c0);
        [Mreal Mimag; Mimag' Mreal]==semidefinite(2*(m+1),2*(m+1));
        Pimag+Pimag'==zeros(m,m);
        [Preal Pimag;Pimag' Preal] == semidefinite(2*m, 2*m);
        [real(T)-Tma_real  imag(T)-Tma_imag; (imag(T)-Tma_imag)'  real(T)-Tma_real]==semidefinite(2*n,2*n);
        cvx_end
        i=sqrt(-1);
        Tma=full(Tma_real+i*Tma_imag);
        Ts= T-Tma;
    end
else
    Tma=min(eig(T))*eye(n);
    Ts=T-Tma;
end
    
    
    
    