 function [omega_ss, residues_ss, sigmashalf, sigmaotherhalf]=sm(P,Ahat,bhat,m)

% function [omegas,residues]=sm(P,A,B,m)   subspace method generalizing ESPRIT
%
% Signal Estimation via Selected Harmonic Amplification
% Subspace Method 2 in [1] generalizing ESPRIT
%
% Employs an input-to-state filter SYS to generate covariance data.
%
% INPUTS:
%          A,B : state-space data of a filter G(z)= (I - z^-1 A)^-1 B
%                normalized so that A*A'+B*B'=I
%
%            P : state-covariance when the above filter is driven by a stochastic
%                process whose spectrum we seek to analyze
%
%            m : number of anticipated spectral lines
%                (= number of lines produced by the current sub-space analysis)
%
% OUTPUTS:
%     omegas : frequency estimates based on analysis of the signal subspace
%
%   residues : corresponding residues
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Generate a basis for the signal subspace
   [U,svds]=svd(P); Ua=U(:,1:m);
   sigmashalf=diag(svds(1:m,1:m)).^.5; n=length(P); sigmaotherhalf=diag(svds(m+1:n,m+1:n)).^.5;

   phi=[bhat Ahat*Ua]\Ua; %mu=phi(1,:); 
   phi(1,:)=[];    % equation (33) in [1]
   omega_ss=angle(eig(phi)); 
   
   D=diag(exp(1i*omega_ss));
   Gs=axxbc(-Ahat,D',bhat*(diag(D)'));                 % equation (32) in [1]

   Rhalf=Gs\(Ua*diag(sigmashalf));
   residues=diag(Rhalf*Rhalf');

   [residues,order]=sort(residues);

   orderedthetas=omega_ss(order);

   omega_ss=orderedthetas(1:m);      omega_ss=flipud(omega_ss);
   omega_ss=-omega_ss+(omega_ss>0)*2*pi; 
   residues_ss=residues(1:m);        residues_ss=flipud(residues_ss);
   residues_ss=sqrt(abs(residues_ss));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Based on:
% [1] T.T. Georgiou, ``Spectral Estimation via Selective Harmonic Amplification,''
%                      in IEEE Trans. on Automatic Control
%
%%%%%%% Written by Tryphon Georgiou -- last line of sm.m  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
