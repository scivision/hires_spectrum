function w=music(y,n,m)
%
% The Root MUSIC method for frequency estimation.
%
%  w=music(y,n,m);
%
%      y  ->  the data vector
%      n  ->  the model order
%      m  ->  the order of the covariance matrix in (4.5.14)
%      w  <-  the frequency estimates
%

% Copyright 1996 by R. Moses

y=y(:);
N=length(y);                       % data length

% compute the sample covariance matrix

R=zeros(m,m);
for i = m : N,
   R=R+y(i:-1:i-m+1)*y(i:-1:i-m+1)'/N;
end

% to use the forward-backward approach, uncomment the next line
% R=(R+fliplr(eye(m))*R.'*fliplr(eye(m)))/2;

% get the eigendecomposition of R; use svd because it sorts eigenvalues
[U,D,V]=svd(R);
G=U(:,n+1:m);

GG = G*G';

% find the coefficients of the polynomial in (4.5.16)
 a = zeros(2*m-1,1);
  for j=-(m-1):(m-1)
    a(j+m) = sum( diag(GG,j) );
  end

% find the n roots of the a polynomial that are nearest and inside the unit circle,
ra=roots([a]);
rb=ra(abs(ra)<1);

% pick the n roots that are closest to the unit circle
[dumm,I]=sort(abs(abs(rb)-1));
w=angle(rb(I(1:n)));
