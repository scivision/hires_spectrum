"""
 Maximum Entropy Signal Analysis

 [spectrum, halfspectrum, omega]=me(P,A,B,omega)

  computes the maximum entropy spectrum
  based on a state covariance P for the state model (A,B)


 INPUTS: P     state covariance
         A,B   filter G(z)= (I - z^-1 A)^-1 B, A,B possibly complex,
               normalized so that I = AA' + BB'
         omega frequency range as an array of values 0 < omega < 2*pi
               (default: omega=linspace(0,2*pi,200))

 OUTPUTS:
   spectrum  power spectrum for input of maximal entropy
             in case the input is vectorial the spectrum is packed
             horizontally:

             spectrum = [ f11(0) f12(0) etc. f1m(0)    f11(Dtheta) ...
                          f21(0) ............f2m(0)    f21(Dtheta) ...
                          ...
                          fm1(0)             fmm(0     fm1(Dtheta) ... ]

halfspectrum is the square-root of the spectrum
   omega     is the range of frequencies

 Based on:

 [1] T.T. Georgiou, ``Spectral analysis based on the state covariance:
                      state-space formulae''
                      IEEE Trans. on Automatic Control.

 Written by Tryphon Georgiou -- last line (March 23, 2001, revised 5/2005,
 revised March, 2009)
"""
from numpy import linspace

def maxent(P,A,B,omega):
    if nargin==3:
        omega=linspace(0,2*pi,200)

    N=len(omega)

    iP=inv(P)
    e=inv(B.T*iP*B)
    C1=e*B.T*iP
    n,m = B.shape

# Phi=C1*G(z); max entropy spectrum = Phi^{-1} e Phi^{-1}^*

    rA = A.real
    iA = A.img
    AA=[rA -iA; iA rA]
    rB=B.real
    iB=B.imag
    BB=[rB -iB; iB rB]
    rC=C1.real
    iC=B1.imag
    CC1=[rC -iC; iC rC]

    Phi_inv=ss(AA-BB*CC1*AA,BB,-CC1*AA,eye(2*m,2*m),1);
    HH=freqresp(Phi_inv,omega); H=HH(1:m,1:m,:)+1i*HH(m+1:2*m,1:m,:);

    if m==1:
        h=squeeze(H)
        spectrum=real(diag(h*e*h.T)).T
    else:
        print('The spectrum is matricial m*m, where m={}'.format(m))
        spectrum = zeros(m,N*m)
        halfspectrum = zeros(m,N*m)
        [Ue,svdOmega]=svd(e)
        sqrtsvdOmega=sqrt(svdOmega)
        for i in range(N):
            temp=H(:,:,i); %H(:,:,i)=temp*e*temp.T
            spectrum(:,(i-1)*m+1:i*m)=(temp*e*temp.T)
            halfspectrum(:,(i-1)*m+1:i*m)=temp*Ue*sqrtsvdOmega

    return spectrum, halfspectrum, omega