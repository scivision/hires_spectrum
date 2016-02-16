function [P, x, xb]=dlsim_complex(A,B,u)

% [P, x, xb]=dlsim_complex(A,B,u): state cov statistics with vectorial input u
%
% Computes state statistics when a discrete-time input u
% drives the system x_{k+1}=Ax_k+Bu_k
%
% u = [ u(1) u(2) ... ] where u(k) could be a column vector
% and all of A,B,u can be complex-valued
%
% Input:  A,B state-space data normalized so that A*A'+B*B'=I
%             (needed for well-condintioning and for reliable estimation of P)
%
%           u input
%
% Output:   P estimated state covariance
%
%        x,xb state trajectory when the system is driven in forward and reverse sense
%
% Comments:
% In estimating P, a transient segment is erased based on the system timeconstant.
% The choice of ``timeconstant'' affects the estimate of P.

timeconstant=max([floor(-1/log10(max(abs(eig(A)+eps)))) 1]);

[n,m]=size(B);
rA=real(A);  iA=imag(A);  AA=[rA -iA; iA rA];
rB=real(B);  iB=imag(B);  BB=[rB -iB; iB rB];
ru=real(u); ru=detrend(ru')';
iu=imag(u); iu=detrend(iu')';
U=[ru;iu]; fliplrU=fliplr([ru;iu]);

[x,xdelayed]=dlsim(AA,BB,AA,BB,U);           x=x.'; x=x(1:n,:)+1i*x(n+1:2*n,:);
[xb,xbdelayed]=dlsim(AA,BB,AA,BB,fliplrU); xb=xb.'; xb=xb(1:n,:)+1i*xb(n+1:2*n,:);

[n,N]=size(x);
    x=x(:,timeconstant:N); xb=xb(:,timeconstant:N);
    Pf=x*x'/(N+1); Pb=xb*xb'/(N+1); P=(Pf+Pb)/2;
%%% last line of dlsim_complex.m
