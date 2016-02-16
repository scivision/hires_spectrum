function [Ahat,bhat,A,b]=rjordan(n,lambdas)

% [Ahat,bhat,A,b]=rjordan(n,lambdas)
% 
% n=[n1 n2 ...]; lambdas=[lambda1 lambda2 ...]
%
% Then
% A,b : a controllable pair (provided lambda_i \neq lambda_j) with
%       A is a Jordan matrix with n(i)-size blocks having eigenvalues lambdas(i).
%
%       The eigenvalues lambdas(i) can be complex, in which case the system
%       has as eigenvalues their conjugates as well
%
% Ahat,bhat : are normalized so that Ahat*Ahat'+bhat*bhat'=I
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=[];b=[];

for i=1:length(n),
    if isreal(lambdas(i)),
a1=poly(zeros(1,n(i)));
A1=compan(a1); b1=eye(length(A1),1);
A1=A1+lambdas(i)*eye(size(A1));
[rs,cs]=size(A);
[r1,c1]=size(A1);
A=[A zeros(rs,c1); zeros(r1,cs) A1];
b=[b;b1];
    else
        Block=[real(lambdas(i)) imag(lambdas(i));-imag(lambdas(i)) real(lambdas(i))];
        a1=poly(zeros(1,n(i))); A1=compan(a1); %b1=eye(length(A1),1);
        A1r=kron(eye(n(i)),Block)+kron(A1,eye(2)); b1r=ones(length(A1r),1);
        [rs,cs]=size(A);
        [r1r,c1r]=size(A1r);
        A=[A zeros(rs,c1r); zeros(r1r,cs) A1r];
        b=[b;b1r];
    end
end

R=dlyapchol(A,b); invRprime=inv(R)'; Ahat=invRprime*A*R'; bhat=invRprime*b;