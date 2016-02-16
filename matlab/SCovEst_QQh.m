function [R,Q,Qh]=SCovEst_QQh(Rh,A,B)
% function [R,Q,Qh]=SCovEst_QQh(Rh,A,B)
% 
% R+Q=Rh+Qh  (where R>0, R-ARA'=BH+H'B', Q,Qh>0 and tr(Q+Qh) is minimal)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n=length(Rh); [n1,m]=size(B);

cvx_begin quiet
variable Q(n,n) symmetric
variable R(n,n) symmetric
variable Qh(n,n) symmetric
variable H(m,n1)
minimize(trace(Q+Qh))
Q == semidefinite(n,n);
Qh == semidefinite(n,n);
R == semidefinite(n,n);
R-A*R*A' == B*H+H'*B';
R+Q == Rh+Qh;
cvx_end

