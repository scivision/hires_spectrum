clear
%% test signal
Anoise=0.01;
A1=1; 
A2=1;
f1=1000; %[Hz]
f2=1050; %[Hz]

fs = 8000; %[H[]

t = 0:1/fs:50e-3; t=t(:);
Ns = length(t);

df = fs/Ns;
fax= 0:df:fs/2;

disp(['unpadded freq. resolution ',num2str(df)])

y = Anoise*randn(Ns,1) + ...
    A1*exp(1j*(2*pi*f1*t)) + ... %+2*pi*rand)) + ... %random phase
    A2*exp(1j*(2*pi*f2*t));%+2*pi*rand)); %random phase

%% maximum entropy
fmid = mean([f1,f2]);
thetamid = 2*pi * fmid / (fs/2); 
[A,B]=cjordan(5,0.88*exp(thetamid*1j));
R=dlsim_complex(A,B,y');
spectrum=me(R,A,B);
%% plot
figure(1); clf(1)
plot(fax(1:end-1),spectrum)
xlabel('freq [Hz]'), ylabel('relative amplitude')
set(gca,'xlim',[0.9*fmid, 1.1*fmid])

figure(2);clf(2)
plot(t,real(y))
xlabel('time [sec.]'), ylabel('amplitude')