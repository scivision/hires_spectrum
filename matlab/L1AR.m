function [omega,mag,phi,a,sigma]=L1AR(y,freqNum,ARorder,freqrange,N)
% function [omega,mag,phi,a,sigma]=L1AR(y,freqNum,ARorder,freqrange,N)
% input
%              y: observation record signal+noise
%        freqNum: number of signals
%        ARorder: order of autoregressive filter (default value 0)
%      freqrange: interested frequency interval (default value [0 pi])
%              N: resolution=2*pi/N (default value 2*length(y))
% output
%          vopt: 
%         omega:  estimated frequencies of signals
%           mag:  magnitudes of signals
%           phi:  estimated phase angles
%             a:  coefficients of AR filter
%         sigma:  whitened noise variance

n=length(y);
if nargin==2
    ARorder = 0;
    freqrange = [0 pi];
    N = 2*n;
else if(nargin==3) 
    freqrange = [0 pi];
    N = 2*n; 
    else if(nargin==4)
        N = 2*n; 
        end
    end
end
%% build weight
[Ba,freq,A]=freqselect(n-ARorder,freqrange,N);
    for i=1:N
        B(:,i)=Ba(:,i)/norm(Ba(:,i));
    end
    W=eye(N);
    Hy=hankel(y(1:n-ARorder),y(n-ARorder:n));
    Hy=fliplr(Hy);
if(ARorder>0)
    H=Hy(:,2:end);
    P=H*inv(H'*H)*H';
    Pperp=eye(size(P))-P;
    Bb=Pperp*B;
    Q=((diag(diag(Bb'*Bb))));
    w=diag(Q);
    w=sqrt(.5*w(1:N/2).^2+.5*w(N/2+1:N).^2);
    W=diag([w; w]);
end
%% optimize 
% differet results could be obtained by changing the coefficient rho 
% larger rho-->sparser vs
rho=1;

cvx_begin quiet
    variable as(ARorder);
    variable vs(N,1);
minimize(rho*sum(abs(W*vs))+.5*(Hy*[1;as]-B*vs)'*(Hy*[1;as]-B*vs)); 
cvx_end
a=[1;as]';
vopt=sqrt(vs(1:N/2).^2+vs(N/2+1:N).^2);
for iter=1:8;
    weight=2./(vopt+0.00000001);
    W=diag([weight;weight]);
    cvx_clear
    cvx_begin quiet
    variable as(ARorder);
    variable vs(N,1);
    minimize(rho*sum(abs(W*vs))+.5*(Hy*[1;as]-B*vs)'*(Hy*[1;as]-B*vs));
    cvx_end
    a=[1;as]';
    vopt=sqrt(vs(1:N/2).^2+vs(N/2+1:N).^2)/sqrt(2);
end
   figure(1);stem(abs(vopt));

%% signal parameters
[Vest,ind]=sort(vopt,'descend');
omega=freq(ind(1:freqNum));
phi=[];
Pj=[cos([0:n-1]'*omega) sin([0:n-1]'*omega)];
c=(Pj'*Pj)\Pj'*y;
for k=1:freqNum
    phi(k)= atan(c(k)/c(k+freqNum));
    mag(k)= norm([c(k) c(k+freqNum)]);
end
[omega, ind]=sort(omega);
phi=phi(ind);
mag=mag(ind);
%% noise level
Bsub=B(:,[ind(1:freqNum) ind(1:freqNum)+N/2]);
vsub=vs([ind(1:freqNum)  ind(1:freqNum)+N/2]);
noise=Hy*[1;as]-Bsub*vsub;
sigma=norm(noise)/sqrt(length(noise));




