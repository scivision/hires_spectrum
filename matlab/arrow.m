function arrow(omega,ampl)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% function arrow(omega,ampl)
%%
%% Draws arrows of amplitude specified by ampl
%% and location specified by omega
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin==1, ampl=ones(size(omega)); end
if size(omega)~=size(ampl), disp('amplitudes/frequencies inconsistent'), return, end

hold on

for i=1:length(omega)
 plot([omega(i) omega(i)],[0 ampl(i)],'r');
 plot(omega(i), ampl(i),'r^');
end

%%%%%%%%%%%%%%%%%%%%%%%%% Last line of arrow.m %%%%
%%%%%%%%%%%%%%%%%%%%%%%%% TTG: March 4, 1999 %%%%%%
