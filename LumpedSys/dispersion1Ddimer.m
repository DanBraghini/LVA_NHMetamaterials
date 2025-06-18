% Controlled 1D lattice dimer
clear all
close all
clc
%%

%Parameters
k1=1; m1=1;
k2=1; m2=2;
rm=m2/m1;
a1=0; a2=0; %How many masses behind feedback comes from for both masses
kc1=-0.1; kc2=-0.1; % control gains

muvec=pi*linspace(-1,1,701);
omega=zeros(2,length(muvec)); 
Modes=zeros(2,2,length(muvec));
M=diag([m1 m2]);

%------------------------------------------------------------------------------------------------------------------%
%------------------Periodic (Infinity) System : Dispersion Relation--------
%----------------------------------------------------------------------------------------------------------------%

%auxiliar variables
a1p=ceil(a1/2); a1pp=ceil((a1+1)/2);
% if a1 is odd, taking a1=1 for instance, a1p=1=a1pp
% if a1 is even, taking a1=2 for instance, a1p=1 and a1pp=2
a2p=floor(a2/2); a2pp=floor((a2+1)/2);
% if a2 is odd, taking a2=1 for instance, a2p=0 and a2pp=1
% if a1 is even, taking a2=2 for instance, a2p=1=a2pp


for ii=1:length(muvec)
    mu=muvec(ii);
    
    K=[ (k1+k2)                  (-k1-k2*exp(-1i*mu))
       (-k1-k2*exp(1i*mu))                      (k1+k2)];
   
    Kc=diag([kc1 kc2])*[(-1)^a1*exp(-1i*mu*a1p)    (-1)^(a1+1)*exp(-1i*mu*a1pp)
                       (-1)^(a2+1)*exp(-1i*mu*a2p) (-1)^a2*exp(-1i*mu*a2pp)];
    
    Kt=K+Kc;
    [V,lambda] = eig(Kt,M);
    omega(:,ii) = sqrt(diag(lambda));
    % Sort complex numbers according to their real part
    % Modes is the mode shapes, corresponding to harmonic excitation in
    % steady state, on the left of structure
    [omega(:,ii),ind] = sort(omega(:,ii), 'ComparisonMethod','real');
    Modes(:,:,ii) = V(:,ind); 
end

% Comparing with passive dimer

muvec_nofeedback=pi*linspace(-1,1,701);
omega_nofeedback=zeros(2,length(muvec_nofeedback)); 
Modes_nofeedback=zeros(2,2,length(muvec_nofeedback));

for j=1:length(muvec_nofeedback)
    K=[k1+k2 -k1-k2*exp(-1*i*muvec_nofeedback(j))
        -k1-k2*exp(1*i*muvec_nofeedback(j)) k1+k2];
    
    [V_nofeedback,lambda_nofeedback]=eig(K,M);
    omega_nofeedback(:,j)=sqrt(diag(lambda_nofeedback));
    [~,ind]=sort(real(omega_nofeedback(:,j)));
    omega_nofeedback(:,j)=omega_nofeedback(ind,j);
    Modes_nofeedback(:,:,j)=V_nofeedback(:,ind); 
end
%---Periodic (Infinit) System : Dispersion Relation------%

figure(1)
plot3(real(omega(1,:)),imag(omega(1,:)),muvec/pi,'r','LineWidth',1.5)
hold on
plot3(real(omega(2,:)),imag(omega(2,:)),muvec/pi,'b','LineWidth',1.5)
plot3(real(omega_nofeedback(1,:)),imag(omega_nofeedback(1,:)),muvec_nofeedback/pi,'k--','LineWidth',1.5)
plot3(real(omega_nofeedback(2,:)),imag(omega_nofeedback(2,:)),muvec_nofeedback/pi,'k--','LineWidth',1.5)
ylabel('\Im(\omega)')
xlabel('\Re({\omega})')
zlabel('\mu/ \pi')
box on
grid on
title('Dispersion Diagram (real part)')
% Dispersion realation
view(0,0) 

figure(3)
plot3(real(omega(1,:)),imag(omega(1,:)),muvec/pi,'r','LineWidth',1.5)
hold on
plot3(real(omega(2,:)),imag(omega(2,:)),muvec/pi,'b','LineWidth',1.5)
plot3(real(omega_nofeedback(1,:)),imag(omega_nofeedback(1,:)),muvec_nofeedback/pi,'k--','LineWidth',1.5)
plot3(real(omega_nofeedback(2,:)),imag(omega_nofeedback(2,:)),muvec_nofeedback/pi,'k--','LineWidth',1.5)
ylabel('\Im(\omega)')
xlabel('\Re({\omega})')
zlabel('\mu/ \pi')
box on
grid on
title('Temporal Atenuation Diagram')
% Imaginary frequency
view(90,0)

figure(4)
plot3(real(omega(1,:)),imag(omega(1,:)),muvec/pi,'r','LineWidth',1.5)
hold on
plot3(real(omega(2,:)),imag(omega(2,:)),muvec/pi,'b','LineWidth',1.5)
ylabel('\Im(\omega)')
xlabel('\Re({\omega})')
zlabel('\mu/ \pi')
box on
grid on
title('Complex Frequency Dispersion Band')

%--------------Finite System: Natural Frequencies------------%

N=300;
K=zeros(N,N);

K(1,1)=k1; K(1,2)=-k1;

mvec=zeros(N,1);
mvec(1:2:N,1)=m1; mvec(2:2:N,1)=m2;
M=diag(mvec);
%even masses (u2, a2)
for ii=2:2:N-1
    K(ii,ii)=k1+k2;
    K(ii,ii-1)=-k1;
    K(ii,ii+1)=-k2;
    if ii>a2+1
        K(ii,ii-a2)=K(ii,ii-a2)+kc2;
        K(ii,ii-a2-1)=K(ii,ii-a2-1)-kc2;
    end
end
%odd masses (u1, a1)
for ii=3:2:N-1
    K(ii,ii)=k1+k2;
    K(ii,ii-1)=-k2;
    K(ii,ii+1)=-k1;
    if ii>a1+1
        K(ii,ii-a1)=K(ii,ii-a1)+kc1;
        K(ii,ii-a1-1)=K(ii,ii-a1-1)-kc1;
    end
end
%last mass
if mod(N,2)==0 %if N is even 
    K(N,N)=k1; K(N,N-1)=-k1;
    K(N,N-a2)=K(N,N-a2)+kc2; K(N,N-a2-1)=K(N,N-a2-1)-kc2;
else
    K(N,N)=k2; K(N,N-1)=-k2;
    K(N,N-a1)=K(N,N-a1)+kc1; K(N,N-a1-1)=K(N,N-a1-1)-kc1;
end

% Normal Modes for active structure (free response)
[V,Lambda]=eig(K,M);
[wn,ind]=sort(sqrt(diag(Lambda)));
Wn=wn*sqrt(m1)/sqrt(k1);
V=V(:,ind);

K_nofeedback=toeplitz([2*k1 -k1 zeros(1,N-2)]);
K_nofeedback(1,1)=k1;
K_nofeedback(N,N)=K_nofeedback(1,1);

% Normal Modes for passive structure (free response)
[V_nofeedback,Lambda]=eig(K_nofeedback,M);
[wn_nofeedback,ind]=sort(sqrt(diag(Lambda)));
Wn_nofeedback=wn_nofeedback*sqrt(m1)/sqrt(k1);
V_nofeedback=V_nofeedback(:,ind);

%  
% figure(2)
% plot3(real(Omega(1,:)),imag(Omega(1,:)),muvec/pi,'r','LineWidth',1.5)
% hold on
% plot3(real(Omega(2,:)),imag(Omega(2,:)),muvec/pi,'b','LineWidth',1.5)
% plot3(real(Omega_nofeedback(1,:)),imag(Omega_nofeedback(1,:)),muvec_nofeedback/pi,'k--','LineWidth',1.5)
% plot3(real(Omega_nofeedback(2,:)),imag(Omega_nofeedback(2,:)),muvec_nofeedback/pi,'k--','LineWidth',1.5)
% ylabel('\Im(\Omega)')
% xlabel('\Re({\Omega})')
% zlabel('\mu/ \pi')
% box on
% grid on
% title('Complex Frequency Plane')
% view(0,90)
% 
% 
% hold on
% plot(real(Wn),imag(Wn),'k.','MarkerSize',14)
% mode1=10;mode2=60;
% 
% figure(5)
% display(Wn(mode1))
% plot(1:1:N, V(:,mode1),1:1:N, V_nofeedback(:,mode1));
% xlabel('Particle Nr.')
% ylabel('û_n(\omega) ');
% legend(['\omega_n = ' num2str(wn(mode1)) 'rad/s'], ['\omega_n = ' num2str(wn_nofeedback(mode1)) 'rad/s']);
% title('First Pass Band Modes')
% figure(6)
% display(Wn(mode2))
% plot(1:1:N, V(:,mode2),1:1:N, V_nofeedback(:,mode2));
% xlabel('Particle Nr.')
% ylabel('û_n(\omega) ');
% legend(['\omega_n = ' num2str(wn(mode2)) 'rad/s'], ['\omega_n = ' num2str(wn_nofeedback(mode2)) 'rad/s']);
% title('Second Pass Band Modes')
% 
% mode3=find((abs(Wn-1.204) < (Wn(2)-Wn(1))/10));
% figure(7)
% display(Wn(mode3))
% plot(1:1:N, V(:,mode3),1:1:N,V_nofeedback(:,mode3));
% xlabel('Particle Nr.')
% ylabel('û_n(\omega) ');
% legend(['\omega_n = ' num2str(wn(mode3)) 'rad/s'], ['\omega_n = ' num2str(wn_nofeedback(mode3)) 'rad/s']);

%---------- Forced Response using Lsim------------------------%
% state vector in R^(N x N) , v = [x(t) x_dot(t)]^T

Ce = 0.*K;
A = [zeros(N,N) eye(N);
    -inv(M)*K   -inv(M)*Ce];
A_passive = [zeros(N,N)               eye(N);
              -inv(M)*K_nofeedback   -inv(M)*Ce];
B = [zeros(N,N);
     inv(M)];
C = [eye(N) zeros(N,N)];
D = 0;

% Defining tone-burst excitation
% Nc is the number of circles from central frequency fc.
% T2 define the envelope frequency f2, which has Nc circles within.
% Nt is the number of entries on time vector t. The heavside function
% defines a window on half a period of the envelope.
% time discretization period dt depends on Tc 

wc=1.5;
fc=wc/2/pi;
Nc=25; 
Tc=1/fc;
T2=Nc*Tc*2; 
f2=1/T2; w2=2*pi*f2;
dt=Tc/60;
Nt=2*9046;
t=0:dt:(Nt-1)*dt;

F=double(sin(wc*t).*sin(w2*t).^2.*heaviside(T2/2-t)); 
df=1/t(end); fm=1/dt;
fvec=0:df:fm/2;

Fw=1/Nt*fft(F);

figure
subplot(1,2,1)
plot(t,F)
title('external force dynamics')
subplot(1,2,2)
plot(2*pi*fvec,abs(Fw(1:Nt/2)))
xlim([0 2])
title('external force spectrum')

U = zeros(N,Nt);
U(round(N/2),:) = F;

ss_dimer = ss(A,B,C,D);
ss_passive = ss(A_passive,B,C,D);
Y=lsim(ss_dimer,U,t,zeros(2*N,1));
Y_passive = lsim(ss_passive,U,t,zeros(2*N,1));

xd=1:1:N;
% transient simulations
k=1;
figure;
plot(xd,Y(1,:))
hold on

for i=2:500:Nt
    plot(xd,Y(i,:)./max(abs(Y(i,:)))+ 2*k,'b')
    k=k+1;
end
title('active system')

figure;
k=1;
plot(xd,Y_passive(1,:))
hold on
for i=2:500:Nt
    plot(xd,Y_passive(i,:)./max(abs(Y_passive(i,:)))+ 2*k,'b')
    k=k+1;
end
title('passive system')

figure
norm = max(max(Y));
for i=1:1000:length(t)
    plot3(xd,t(i)*ones(1,length(xd)),Y(i,:)./norm,'b');
     hold on
end

xlim([0 max(xd)]), ylim([0 max(t)])
xlabel('Length [m]'), ylabel('Time [s]'), zlabel('Displacement u(x,t)[m]')
hold off
title('active system')

norm =max(max(Y_passive));
figure
for i=1:1000:length(t)
    plot3(xd,t(i)*ones(1,length(xd)),Y_passive(i,:)./norm,'b');
     hold on
end

title('passive system')
colormap jet
xlim([0 max(xd)]), ylim([0 max(t)])
xlabel('Length [m]'), ylabel('Time [s]'), zlabel('Displacement u(x,t)[m]')
hold off
 
% 2D ffts
tend=16002; %index of last instant considered for FFts : this must be chosen for each case! Preferebly before reflections to get dispersion curves
U=Y(1:tend,:); 
U_passive=Y_passive(1:tend,:);

Uw=fftshift(fft2(U));
Uw_passive=fftshift(fft2(U_passive));

xsamp=1; Nx=length(xd);
kx=xsamp*(-Nx/2:Nx/2-1)/Nx; kx=2*pi*kx; kx=fliplr(kx);
dt2=t(2)-t(1);
dsamp=1/dt2;
fvec2=dsamp*(-tend/2:tend/2-1)/tend; wvec2=fvec2*2*pi;
figure
pcolor(wvec2,kx/pi,abs(Uw'))
shading interp
xlabel('\omega (rad/s)')
ylabel(' k L_c')
hold on
plot(real(omega(1,:)),muvec/pi,'r');
plot(real(omega(2,:)), muvec/pi,'b');
axis([0 1.8 -1 1])

figure
pcolor(wvec2,kx/pi,abs(Uw_passive'))
shading interp
xlabel('\omega (rad/s)')
ylabel(' k L_c')
hold on
plot(real(omega_nofeedback(1,:)),muvec/pi,'r');
plot(real(omega_nofeedback(2,:)), muvec/pi,'b');
axis([0 1.8 -1 1])

%title('FRF- Receptance ');

% % Direct Method
% Wd=0:.01:2;
% 
% for i=1:length(Wd)
%     D=[(1-Wd(i)^2)*(2-Wd(i)^2*rm-kc2/k1)-1-kc2/k1         -1
%        -(1-kc1/k1)*(1-kc2/k1)                          (1-kc1/k1)*(2-rm*Wd(i)^2-kc2/k1)-1+kc1/k1];
%    D=1/(2-rm*Wd(i)^2-kc2/k1)*D;
%    %D=[1         -1
%    %    0    2-2*Wd(i)^2];
%    %D=inv(D)*[Wd(i)^2+1  -1
%    %              -1     1-2*Wd(i)^2];
%              
%    T=[-D(1,2)^-1*D(1,1)                      D(1,2)^-1
%        D(2,2)*D(1,2)^-1*D(1,1)-D(2,1)   -D(2,2)*D(1,2)^-1];
%    % From Floquet's Theorem
%     if abs(Wd(i))==sqrt((2-kc2/k1)/rm)
%        kL(i)=kL(i-1);
%     else
%    [V,D]=eig(T);
%    Fm1=D(1,1);Fm2=D(2,2);
%    kL(i)=-1i*log(Fm1);
%     end
%     
% end
% 
%         
% figure
% plot(Wd, real(kL)/pi, 'b')
% hold on
% plot(Wd,imag(kL/pi),'--')
% ylabel('kL/ \pi')
% xlabel('\Omega')
% grid on
% title('Dispersion Relation')