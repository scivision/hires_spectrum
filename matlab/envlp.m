function [sqrtrho,omega]=envlp(P,Ahat,bhat,omega,noiselevel)

% [sqrtho,omega]=envlp(P,A,B,omega,noiselevel) : spectral envelop
%
% Computes the envelope of spectra consistent with the data
% by using state covariance statistics P for the input-to-state system (A,B).
%
% If 0 =< NOISELEVEL <1 is assigned a value, then a level of
% white noise at the input is assumed which correspond to the
% specified ``noiselevel'' fraction of the smallest sv of the
% the state covariance (default: NOISELEVEL=0).
%
% Values of NOISELEVEL > 0 tend to increase resolution between spectral lines.
%
% System matrices need to be normalized so that:
% A*A'+B*B' = I
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I=eye(size(Ahat));
W=commutant(P,Ahat,I);
P=W+W';

if nargin==5, P=P-min(eig(P))*noiselevel; end
if nargin==3, omega=linspace(0,2*pi,200); end
N=length(omega);
 
Pinvhalf=inv(sqrtm(P));

expomega=exp(1i*omega);
sv=zeros(1,N);
for i=1:N,
sv(i)=norm( Pinvhalf*inv( I-conj(expomega(i))*Ahat ) *bhat );
end

sqrtrho=1./sv;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%