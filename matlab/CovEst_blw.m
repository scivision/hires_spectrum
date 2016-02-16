function T=CovEst_blw(hatT)
% function T=CovEst_blw(hatT)
% maximum likelihood estimation of Toeplitz covariance matrix
% method proposed by Burg, J.P.;   Luenberger, D.G.;   Wenger, D.L.;  
if(norm(hatT-hatT')>.001|| min(eig(hatT)<0))
    error('Sample covariance must be PSD');
end
%sigmamin=1e-5;
sigmamin=1e-6;

n=length(hatT);
I=eye(n);
Q=zeros(n,n,n);
for i=1:n
   Q(:,:,i)=toeplitz(I(i,:));
end

% start from I

Tk=eye(n);
gk=-log(det(Tk))-trace(Tk\hatT);
stop=0;

while(stop~=1)
    
A=zeros(n,n);
c=zeros(n,1);
for i=1:n
    c(i)=trace(Tk\hatT/Tk*Q(:,:,i));
    for j=1:n
        A(i,j)=trace(Tk\Q(:,:,i)/Tk*Q(:,:,j));
    end
end
t=inv(A)*c;
Ttry=toeplitz(t);
shrinkstep=0;

if(min(eig(Ttry))>sigmamin)
    gtry=-log(det(Ttry))-trace(Ttry\hatT);
    if(gtry>gk+.00001)
        gk=gtry;
        Tk=Ttry;
    else
        shrinkstep=1;
    end
else
    shrinkstep=1;
end
 if(shrinkstep==1)
     D=Ttry-Tk;
     for q=1:10
         Ttry=Tk+.5^q*D;
         if(min(eig(Ttry))>sigmamin)
             gtry=-log(det(Ttry))-trace(Ttry\hatT);
             if(gtry>gk+0.00001)
                 gk=gtry;
                 Tk=Ttry;
                 break;
             end
         end
     end
     if(q==10)
         stop=1;
         T=Tk;
     end
 end
end
                 
