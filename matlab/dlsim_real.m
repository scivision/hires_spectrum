function [P, x, xb]=dlsim_real(A,B,u)

% [P, x, xb]=dlsim_real(A,B,u) : state-cov statistics with vectorial input u
%
% Computes state statistics when a discrete-time input u
% drives the system x_{k+1}=Ax_k+Bu_k
%
% u = [ u(1) u(2) ... ] where u(k) could be a column vector
% and all of A,B,u should be real valued
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
% The choice of ``timeconstant'' affect the estimate of P.

timeconstant=max([floor(-1/log10(max(abs(eig(A)+eps)))) 1]);

%u=u-mean(u); % may not be always appropriate, but most of the times it helps <<<
%[n,m]=size(B);
u=u';
[x,xdelayed]=dlsim(A,B,A,B,u);           x=x.';
[xb,xbdelayed]=dlsim(A,B,A,B,flipud(u)); xb=xb.';

[n,N]=size(x);
    x=x(:,timeconstant:N); xb=xb(:,timeconstant:N);
    Pf=x*x'/(N+1); Pb=xb*xb'/(N+1); P=(Pf+Pb)/2;