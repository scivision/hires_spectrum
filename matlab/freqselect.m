function [CS, frequencies, A]=freqselect(N,freqrange, n)
% [CS, frequencies]=freqselect(N,freqrange,n)
%
% Inputs --
% N : number of data points in time-series
%     also specifies default frequency resolution 2*pi/N
%
% freqrange=[freqrange(1) freqqrange(2)] range of frequencies
%               in rad/sample in [0,pi]
%
% n : optional... specifies a frequency resolution 2*pi/n different from
%     default 2*pi/N
%
% Outputs --
% frequencies = selected range of frequencies
%               freqrange(1):(2*pi/n):freqrange(2)-pi/n
%
% CS : basis of cosines and sines = [real(A) imag(A)]
%
% A : Nxm matrix of the form
%     exp(1i*2*[0:N-1]'*freqrange)

if nargin==2, n=N+1; end

if 2*n+1<N, n=N+1, end

frequencies=freqrange(1)+pi/n:2*(pi/n):freqrange(2)-pi/n;
A=exp(1i*(0:N-1)'*frequencies);
CS=[real(A) imag(A)];
