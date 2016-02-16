
function [rsig,theta]=rsigma(A,b,omega)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [rsig,theta]=rsigma(A,b,omega)
%
% Computes sigma for complex data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta=omega(:);
exptheta=exp(j*theta);
I=eye(size(A));
rsig=[];
for i=1:length(omega),
rsig=[rsig; norm( inv( I-conj(exptheta(i))*A ) *b ) ];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Tryphon Georgiou  (revised May 8, 2000)
