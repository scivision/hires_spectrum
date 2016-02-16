
function v=cvec(M)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% M hermitian is stacked into a column vector
%%   keeping only the top triangular part
%%
%% M11 M12 M13 ...
%%     M22 M23 ...
%% 
%% into  [M11;
%%        M22;
%%        ...
%%        Mnn;
%%        real(M12);
%%        real(M13);
%%        ...
%%        imag(M12);
%%        ...]
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[row,col]=size(M);
v=[];
for i=1:col,
v=[v; M(1:i,i)];
end


n=length(M);
   Mreal=real(M);
   Mimag=imag(M);
   Mtopreal=diag(Mreal); Mtopimag=[];
   for i=2:n;
   Mtopreal=[Mtopreal; Mreal(1:i-1,i)];
   Mtopimag=[Mtopimag; Mimag(1:i-1,i)];
   end

v=[Mtopreal; Mtopimag];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% This the last line of cvec.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TTG, March 20, 2000
