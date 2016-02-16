function [R,Q,Qh]=CovEst_QQh(Rh)
% function [R,Q, Qh]=CovEst_QQh(Rh)
% [trace(Q+Qh), Q, R, Qh]=toepQQh(Rh)  :  Rh+Qh=R+Q, R,Q,Qh>0, R Toeplitz, and trace(Q+Qh)=minimal
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=length(Rh); 

cvx_begin quiet

variable Q(n,n) symmetric;
variable Qh(n,n) symmetric;
variable r(n,1);
R = toeplitz(r);

minimize(trace(Q+Qh));
Q == semidefinite(n,n);
Qh == semidefinite(n,n);
R == semidefinite(n,n);
Rh+Qh==R+Q;
cvx_end