function [Ahat,bhat,A,b]=cjordan(n,lambdas)

% [Ahat,bhat,A,b]=cjordan(n,lambdas)
% 
% Generates a controllable pair (provided A is cyclic), such that
% A is a Jordan matrix with n(i)-size blocks having eigenvalues lambdas(i).
% The eigenvalues lambdas(i) can be complex, in which case the system
% parameters have complex values.
%
% Ahat,bhat are normalized so that [Ahat bhat] is co-isometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=[];b=[];

for i=1:length(n),
a1=poly(zeros(1,n(i)));
A1=compan(a1); b1=eye(length(A1),1);
A1=A1+lambdas(i)*eye(size(A1));
[rs,cs]=size(A);
[r1,c1]=size(A1);
A=[A zeros(rs,c1); zeros(r1,cs) A1];
b=[b;b1];
end

R=dlyapchol(A,b); invRprime=inv(R)'; Ahat=invRprime*A*R'; bhat=invRprime*b;