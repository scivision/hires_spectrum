function [R,Q]=CovEst_Q(Rh)
% function [R,Q]=CovEst_Q(Rh)
% [R, Q]=toepQ(Rh)  :  Rh=R+Q, R,Q>0, R Toeplitz, and trace(Q)=minimal
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=length(Rh);
cvx_begin quiet
variable Q(n,n) symmetric;
variable r(n,1); 
R = toeplitz(r);
minimize(trace(Q))
Q == semidefinite(n,n);
R == semidefinite(n,n);
Rh==R+Q;
cvx_end