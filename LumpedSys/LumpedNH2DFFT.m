%% clear 
clear all 
close all 
clc 
%echo on
%%
% boundary == 0 : free-free (default)
% bundary == 1: fixed-fixed
% boundary == 2 : WABC
boundary=1;
% Laplace variable to define feedback law
s=tf('s');
% flag used for feedback in terms of velocity rather than displacement
%derivative=-1: integral feedback (position)
%derivative=0: position static feedback
%derivative=1: velocity static feedback
%derivative=2: acceleration static feedback
derivative=1;
%%
plotbode=0;
plotforce=1;
%% lumped-parameter system set up
% m = mass 
% k = string elastic coefficient
% range of non-reciprocal coupling, r=0 means next-neighborhoud couple
% e = excitation point for the external 
k=10; m1=1; m2=1;wn=sqrt(k/m1);
%g =[-100:1:-10 -9.9:0.1:0];
gd=-10;
gp=0;%-0.05*k;
gdd=0;gi=0;%-0.1;
r=0;
b=0;%0.01;
e='m';
ndof=1501;
%% Define the feedback law
if derivative==-1 % integral case
    gamma_c=[0 gi];
elseif derivative==0
     gamma_c=[gp 0];
elseif derivative==1
         gamma_c=[gd 0];
elseif derivative==2
         gamma_c=[gdd 0];
end
H=gamma_c(1)+gamma_c(2)/s;
%%
% angular frequency vector limit
flim=3.5e0;
output=function_CL4MSS(ndof,m1,m2,k,b,flim,boundary,e,derivative,gamma_c,H);
%output=function_CL4MSS(ndof,m1,m2,k,b,flim,boundary,'l',derivative,gamma_c,H);
%% verify internal stability of closed loop system
[A,B,~,~]=ssdata(output.CLsys);
[Ap,Bp,~,~]=ssdata(output.OLsys);
% [V, Lambda] = eig(A);
% Lambda = diag(Lambda);
% Vn = C*V;
% so=Lambda;
ndof_n=output.ndof_n;
xd=1:1:ndof_n;
nx_cl=size(A,1);
nxc=nx_cl-output.nx;
%% state-space 4-matrix representation SIMO system for Performance Metric
%C=[eye(ndof_n)  zeros(ndof_n,nx_cl-ndof_n)];
C=[eye(2*(ndof_n)) zeros(2*ndof_n,nx_cl-2*ndof_n)];
D=zeros(size(C,1),size(B,2));
sys_simo=ss(A,B,C,D); 

if plotbode==1
    %flim=10*flim;
    % frequency step, used to compute the performance metric
    df=2e-4;
    % angular frequency vector limit
    wv=2*pi*(df:df:flim);
    % Taking BODE diagrams of first and last DOFs
    %P1=P(1,1);
   C1=[1 zeros(1,ndof_n-1)  zeros(1,ndof_n)];
   CN=[zeros(1,ndof_n-1) 1  zeros(1,ndof_n)];
   sys_1=ss(A,B,C1,0);sys_N=ss(A,B,CN,0);
    % Taking BODE diagrams of first and last DOFs
    [mag,~]=bode(sys_N,wv);qnm(i,:)=mag;
     [mag,~]=bode(sys_1,wv);q1m(i,:)=mag;
    S=2*(trapz(qnm)-trapz(q1m))*df;
    figure
   % plot(fv(1:100:end),20*log10(q10m(1:100:end)),'r*','MarkerSize',8,'LineWidth',1)
    plot(wv,20*log10(q1m),'r','MarkerSize',8,'LineWidth',1)
    hold on
    plot(wv,20*log10(qnm),'k','MarkerSize',8,'LineWidth',1)
    hold off
    ylabel('$ |\hat{q} / \hat{F}| (i \omega)$ [dB]', 'interpreter', 'latex', 'fontsize', 15)
    xlabel('$\omega$ [rad/s]', 'interpreter', 'latex', 'fontsize', 15)
    set(gca,'TickLabelInterpreter','Latex','fontsize',15);
    legend('L','R')
    xlim([0 2*pi*flim])
    box on
    set(gcf, 'Color', 'w');
    set(gca,'TickLabelInterpreter','Latex','fontsize',15);
    grid on
    figure
    plot(wv,q1p,'r','MarkerSize',8,'LineWidth',1)
    hold on
    plot(wv,qnp,'k','MarkerSize',8,'LineWidth',1)
    hold off
    ylabel('$\angle \hat{q} / \hat{F} (i \omega) $', 'interpreter', 'latex', 'fontsize', 15)
    xlabel('$f$ [Hz]', 'interpreter', 'latex', 'fontsize', 15)
    set(gca,'TickLabelInterpreter','Latex','fontsize',15);
    legend('L','R')
    xlim([0 2*pi*flim])
    box on
    set(gcf, 'Color', 'w');
    set(gca,'TickLabelInterpreter','Latex','fontsize',15);
    grid on
    %% state-space 4-matrix representation SISO system for H2 norm
    nw=1;nz=1;
    C=[1 zeros(1,ndof_n-2) -1 zeros(1,size(A,2)-ndof_n)];
    D=zeros(nz,nw);
    sys_siso=ss(A,B,C,D);
    Norm2=norm(sys_siso,2);
end
%% Dispersion
muvec=pi*linspace(-1,1,1000);
        %omega_a1=zeros(length(gp),length(muvec));omega_a2=omega_r;
     for j=1:length(muvec)
            cmu=1-cos(muvec(j));
            emu=1-exp(1i*muvec(j));
            a1=(1-gdd/m1*exp(1i*muvec(j)*r)*emu)*1i;
            a2=2*b/m1*cmu-gd/m1*exp(1i*muvec(j)*r)*emu;
            a3=(gp/m1*exp(1i*muvec(j)*r)*emu-2*k/m1*cmu)*1i;
            a4=gi/m1*exp(1i*muvec(j)*r)*emu;
            omega_r(:,j)=roots([a1 a2 a3 a4]);
    end
    

%% zero state response of the closed loop for input w
    wc=0.5;
    fc=wc/2/pi;
    Nc=5; 
    Tc=1/fc;
    T2=Nc*Tc*2; 
    f2=1/T2; w2=2*pi*f2;
    dt=Tc/100;
    Nt=8*9046;
    t=0:dt:(Nt-1)*dt;
    %tend=round(length(t)/10);
    tend=length(t);
    w=double(sin(wc*t).*sin(w2*t).^2.*heaviside(T2/2-t)); 
    dff=1/t(end); fm=1/dt;
    ffvec=0:dff:fm/2;
    Fw=1/Nt*fft(w);
    tau=t*wn;Omega=2*pi*ffvec/wn;
    ind=find(abs(tau-650)<=1);
    plotend=ind(1);
    %%
    if plotforce==1
        figure
        subplot(1,2,1)
        plot(tau(1:plotend),w(1:plotend))
        xlabel('$\tau$', 'interpreter', 'latex','fontsize', 15)
        %xaxis([0 100])
        box on
        grid on
        title('External Force Dynamics', 'interpreter', 'latex','fontsize', 15)
        subplot(1,2,2)
        plot(Omega,abs(Fw(1:Nt/2)))
        xlabel('$\Omega $', 'interpreter', 'latex','fontsize', 15)
        xlim([0 1])
        title('External Force Spectrum', 'interpreter', 'latex','fontsize', 15)
        box on
        grid on
        set(gcf, 'Color', 'w');
    end
    %%
    
    ni=size(A,1);
    z=lsim(sys_simo,w,t,zeros(ni,1));
    for i=1:length(t)
        E(i)=0;
        for j=1:ndof_n
            E(i)=E(i)+m1*z(i,j+ndof_n)^2/2+k*z(i,j)^2/2;
        end
    end
    Erms=sqrt(trapz(E.^2)*dt);
    %passive system
    ni=size(Ap,1);
    zp=lsim(ss(Ap,Bp,eye(size(Ap)),zeros(size(Ap,1),size(Bp,2))),w,t,zeros(ni,1));
    for i=1:length(t)
        Ep(i)=0;
        for j=1:ndof_n
            Ep(i)=Ep(i)+m1*zp(i,j+ndof_n)^2/2+k*zp(i,j)^2/2;
        end
    end
    Ermsp=sqrt(trapz(Ep.^2)*dt);
    % Energy balance
    Ec=Erms-Ermsp;
%     figure
%         plot(t(1:tend),z(1:tend,1),'b')
%         hold on
%         plot(t(1:tend),z(1:tend,ndof_n),'k')
%         set(gcf, 'Color', 'w');
%         xlim([0,t(tend)])
%         xlabel('$t [s]$','interpreter', 'latex','FontSize',15');
%         legend('$q_1(t) [m]$','$q_N(t) [m]$','interpreter', 'latex','FontSize',15');
%         grid on
        
%%
figure
normm = norm(z(1,:),inf);
c = normm.*ones(1,length(xd));
%c=z(1,:);
scale=1/normm*1e1;
patch([xd NaN],[z(1,:)*scale+tau(1)  NaN],[c NaN],'edgecolor','interp')


hold on
DI=round(plotend/20);
for i=1:DI:plotend
    normm = norm(z(i,:),inf);
    scale=1/normm*1e1;
    c = normm.*ones(1,length(xd));
        patch([xd NaN],[z(i,:)*scale+tau(i)  NaN],[c NaN],'edgecolor','interp')

    colorbar
end
colormap winter
box on
xlabel('$n$','interpreter', 'latex','FontSize',15); 
ylabel('$\tau$','interpreter', 'latex','FontSize',15');
h=colorbar ;
ylabel(h,'$u_n(t)$','interpreter', 'latex','FontSize',15)
set(gca,'TickLabelInterpreter','Latex','fontsize',15);
set(gcf, 'Color', 'w');
%export_fig fig4d.pdf
%title('displacement')
    %%
      figure
     for i=1:1:plotend
            normm = max(abs(z(i,:)));
            zp(i,:)=abs(z(i,:))/normm;
     end
    mesh(xd,tau(1:plotend),z(1:plotend,:));
    set(gcf, 'Color', 'w');
    colormap jet
    h=colorbar ;
    %ylabel(h,'$|p(x,t)|$','interpreter', 'latex','FontSize',15)
    hold on
    xlabel('$n$','interpreter', 'latex','FontSize',15); 
    ylabel('$\tau$','interpreter', 'latex','FontSize',15');
    view([0,90])
    tlim=tau(plotend);
    ylim([0 tlim])
box on
%%
%-------------------------------------------------------------------------%
%    Dispersion Diagram via 2D FFT of the transient response         %
%-------------------------------------------------------------------------%
%tend=round(length(t)/2.5); %index of last instant considered for FFts : this must be chosen for each case! Preferebly before reflections to get dispersion curves
ind=find(abs(tau-600)<=1e-1);
tend=ind(1);
U=z(1:tend,:); 
Uw=fftshift(fft2(U));
% total length of x domain
Nx = length(xd);
%dx = Le_s;
%fxsamp=1/dx;
fxsamp=1;
dkx=fxsamp/Nx;
% Nyquist frequency
kxN=fxsamp/2;
kx=-kxN+dkx/2:dkx:kxN-dkx/2;
%kx=-1/2:1/2-1;
kx=kx*2*pi; 
kx = fliplr(kx);
% sample frequency for t domain
%dt=t(2)-t(1);
ftsamp=1/dt;
% total length of t domain = tend
% resolution of frequency f mapped from t
dff = ftsamp/tend;
% Nyquist frequency
fN = ftsamp/2;

ffvec = -fN+dff/2:dff:fN-dff/2;
wvec=2*pi*ffvec;

%%
%p=0:.1:1;
figure
pcolor(kx/pi,wvec/wn,abs(Uw))
shading interp
ylabel('$\Omega_r$','interpreter', 'latex','FontSize',15')
xlabel('$ \mu/ \pi$','interpreter', 'latex','FontSize',15')
hold on
%plot(real(omega_r(1,:)),muvec/pi,'b');
plot(muvec/pi,real(omega_r(2,:))/wn,'k');
plot(muvec/pi,real(omega_r(3,:))/wn,'--k');
%axis([0 1.8 -1 1])
    box on
    grid on
    set(gcf, 'Color', 'w');
%%
% figure;
% plot(real(omega_r(1,:)), muvec/pi,'b',real(omega_r(2,:)), muvec/pi,'k',real(omega_r(3,:)), muvec/pi,'--k');
% xlabel('$\Re(\omega)$ [rad/s]','interpreter', 'latex','FontSize',15')
% ylabel('$ \mu$','interpreter', 'latex','FontSize',15')
%     box on
%     grid on
%     set(gcf, 'Color', 'w');
% figure;
% plot(imag(omega_r(1,:)), muvec/pi,'b',imag(omega_r(2,:)), muvec/pi,'k',imag(omega_r(3,:)), muvec/pi,'--k');
% xlabel('$\Im(\omega)$ [rad/s]','interpreter', 'latex','FontSize',15')
% ylabel('$ \mu$','interpreter', 'latex','FontSize',15')
%     box on
%     grid on
%     set(gcf, 'Color', 'w');
    %%
 %   omega_n=so.*1i;
 %   [omega_n,~]=sort(omega_n,'ComparisonMethod','real');
%     figure;
% %                 omega1(j)=-gd/2/m1*emu*1i+1/2*sqrt(-gd^2*emu^2/(m1^2)+8*k/m1*cmu);
% %             omega2(j)=-gd/2/m1*emu*1i-1/2*sqrt(-gd^2*emu^2/(m1^2)+8*k/m1*cmu);
% plot(real(omega_r(2,:)),imag(omega_r(2,:)),'k',real(omega_r(3,:)), imag(omega_r(3,:)),'--k',real(omega_r(1,:)),imag(omega_r(1,:)),'b');hold on
%  %   scatter(real(omega_n),imag(omega_n),'kx','LineWidth',1.5)
% xlabel('$\Re(\omega)$ [rad/s]','interpreter', 'latex','FontSize',15')
% ylabel('$\Im(\omega)$ [rad/s]','interpreter', 'latex','FontSize',15')
%     box on
%     grid on
%     set(gcf, 'Color', 'w');
 