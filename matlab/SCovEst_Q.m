function [R,Q]=SCovEst_Q(Rh,A,B)
% function [R,Q]=SCovEst_Q(Rh,A,B)
%
%  R+Q=Rh  (where R>0, R-ARA'=BH+H'B', Q>0 and tr(Q) is minimal)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n=length(Rh); [n1,m]=size(B);

cvx_begin quiet

variable Q(n,n) symmetric
variable R(n,n) symmetric
variable H(m,n1)
minimize(trace(Q))
Q == semidefinite(n,n);
R == semidefinite(n,n);
R+Q==Rh;
R-A*R*A'==B*H+H'*B';
cvx_end
