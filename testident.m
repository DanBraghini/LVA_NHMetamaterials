% test ident erac
% run FRF before
close all
figure(1)
plot(fv,20*log10(abs(P_SEM_passive(:,1))))
Nv=length(fv);
H=P_SEM_passive(:,1);
H=[real(H(1));H(1:Nv);0;flipud(conj(H(1:Nv)))];
N=length(H);
df=fv(2)-fv(1);
h=real(N*ifft(H));
T=1/(N*df);
%T=1/(2*1500);
t=(0:(N-1))*T;
figure(2)
plot(t,h)
[Ad,Bd,Cd,Dd]=erasiso(T*h');
sd=ss(Ad,Bd,Cd,Dd,T);
sc=d2c(sd);
[A,B,C,D]=ssdata(sc);
s=eig(A);
sort(abs(s))/2/pi