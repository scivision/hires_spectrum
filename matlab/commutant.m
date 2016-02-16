function [W,wpoly,ProjP,error]=commutant(P,A,E)

% [W,w,ProjP]=commutant(P,A,E) : solves P=(W*E+E*W')/2 with WA=AW
%
% APPLIES to complex data P,A,E
% Solves the equation P=(W*E+E*W')/2
% for a ``commutant'' W of A.
%
% DATA:
% P,A,E are square matrices of equal size,
% P is Hermitian positive semidefinite
% E satisfies the Lyapunov equation E-A*E*A'=b*b'
%   for some b for which (A,b) is controllable
% A is a stable matrix (i.e., with eigenvalues having modulus <1)
%
% OUTPUT:
% W is the commutant
% w is a polynomial so that W=w(A)
% mineig is the least eigenvalue of (W*E+E*W')/2
% relative_approximation_error is || P-(W*E+E*W')/2 ||_{Frobenius} over
% the Frobenius norm of ProjP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=length(P);

Pvec=cvec(P);
Evec=cvec(E);
pack=Evec;

M=[]; AE=E;
   for i=1:length(A)-1,
   AE=A*AE;
   R=real(AE); I=imag(AE);
   Rs=(R+R')/2; Ras=(R-R')/2; Is=(I+I')/2; Ias=(I-I')/2;
   M1=Rs+j*Ias;
   N1=-Is+j*Ras;
   pack=[pack cvec(M1) cvec(N1)];
end

ri=pack\Pvec;
wpoly=[ri(1) ri(2:2:2*n-2)'+j*ri(3:2:2*n-1)'];

W=polyvalm(fliplr(wpoly),A);
ProjP=(W*E+E*W')/2;

DP=P-ProjP;
error=norm(DP,'fro')/norm(ProjP,'fro');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% This the last line of commutant.m