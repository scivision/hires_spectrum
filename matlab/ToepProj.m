function T=ToepProj(S)
% function T=ToepProj(S) project a matrix S into Toeplitz matrix
N=length(S);
t=[];
I=eye(N);
for i=1:N
    D=toeplitz(I(:,i),zeros(N,1));
    r=trace(D*S)/sum(sum(D));
    t=[t r];
end
T=toeplitz(t);