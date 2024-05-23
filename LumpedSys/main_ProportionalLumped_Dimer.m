%% clear 
clear all 
close all 
clc 
%echo on
%%
% boundary == 0 : free-free (default)
% bundary == 1: fixed-fixed
% boundary == 2 : WABC
symmetricdimmer=0;%requires ndof=3*ncell-ncell+1 and boundary=0
boundary=1;
eigenmode=0;
dispersion=0;
FRF3D=0;
s=tf('s');
%% lumped-parameter system set up
% m = mass 
% k = string elastic coefficient
% range of non-reciprocal coupling, r=0 means next-neighborhoud couple
% e = excitation point for the external 
k=1; m1=1; m2=1;gp=0;%-33;
r=0;
b=.10;
e='m';
ndof=9;
% angular frequency vector limit
flim=3.5e0;
wv=2*pi*(0.01:0.002:flim);
fv=wv/2/pi;
%%
% D=[s^2*m+2*c*s+2*k  -s*c-k  0  0
%            -s*c-k+g    s^2*m+2*c*s+2*k-g -s*c-k 0 
%              0               -s*c-k+g        s^2*m+2*c*s+2*k-g -s*c-k
%             0                    0                      -s*c-k-g             s^2*m+2*c*s+2*k-g];
% P=simplify(inv(D));

%% periodic (infinity) system with one mass on the unit cell
if dispersion==1
    muvec=pi*linspace(-1,1,701);
    gp=0;
    omega_a=zeros(size(muvec));
    for i=1:701
        omega_a(i)=sqrt(2*k*(1-cos(muvec(i)) ) - (gp/m1)*(1-exp(1i*muvec(i)) )*exp(1i*muvec(i)*r ) );
    end
    wn=sqrt(k/m1);
    %Omega=omega_a./wn;
%     figure(1)
%     plot3(real(omega_a),imag(omega_a),muvec,'r','LineWidth',1.5)
%     hold on
%     ylabel('$\omega_I$ [rad/s]','interpreter', 'latex', 'fontsize', 15)
%     xlabel('$\omega_R$','interpreter', 'latex', 'fontsize', 15)
%     zlabel('$\mu$','interpreter', 'latex', 'fontsize', 15)
%     box on
%     grid on
%     set(gcf, 'Color', 'w');
%     set(gca,'TickLabelInterpreter','Latex','fontsize',20);
    
    figure(2)
    plot(real(omega_a),muvec,'r','LineWidth',1.5)
    hold on
    xlabel('$\omega_R$','interpreter', 'latex', 'fontsize', 15)
    ylabel('$\mu$','interpreter', 'latex', 'fontsize', 15)
    box on
    grid on
    set(gcf, 'Color', 'w');
    set(gca,'TickLabelInterpreter','Latex','fontsize',20);

    figure(4)
    plot(imag(omega_a),muvec,'r','LineWidth',1.5)
    hold on
    xlabel('$\omega_I$','interpreter', 'latex', 'fontsize', 15)
    ylabel('$\mu$','interpreter', 'latex', 'fontsize', 15)
    box on
    grid on
    set(gcf, 'Color', 'w');
    set(gca,'TickLabelInterpreter','Latex','fontsize',20);
    %%
     figure(1)
    plot(real(omega_a),imag(omega_a),'r','LineWidth',1.5)
    hold on
    xlabel('$\omega_R$','interpreter', 'latex', 'fontsize', 15)
    ylabel('$\omega_I$','interpreter', 'latex', 'fontsize', 15)
    box on
    grid on
    set(gcf, 'Color', 'w');
    set(gca,'TickLabelInterpreter','Latex','fontsize',20);

end
%% gains
gi=0;
gd=0;
gamma_c=[gp gi gd];
% quasiperiodicity parameters
theta=0;phi=0;
%[x,~,so_,Vn,sys_wz,nx_cl] = function_LumpedEigenfreq(ndof,m1,m2,k,eta,flim,boundary,e,gamma_c,theta,phi);
[x_wabc,ndof_wabc,so_wabc,Vn_wabc,sys_wz_wabc,nx_cl_wabc] = function_LumpedEigenfreq(100,m1,m2,k,b,flim,2,e,gamma_c,theta,phi);
%% Computing dynamic stiffness matrix
% number of elements
if ndof>1
    mvec=zeros(ndof,1);
    mvec(1:2:ndof,1)=m1; mvec(2:2:ndof,1)=m2;
    Ms=diag(mvec);
     if symmetricdimmer==1
    % ndof=3*ncell-ncell+1 and boundary==1
          Ms(1,1)=m1/2;Ms(end,end)=Ms(1,1);
     end
    D2_=toeplitz([2 -1 zeros(1,ndof-2)]);
    D2_(1,1)=1; D2_(ndof,ndof)=1;
    Ks=k.*D2_;
   %viscous damping model
    Cs=b.*D2_;
    % Feedback Matrix
    Gs=gp.*(diag(-1*ones(1,ndof))+diag(ones(1,ndof-1),-1));
    Gs(1,:)=zeros(1,ndof);
else
disp('Warning: only one mass isnt enough to have feedback' )
end
% structural/histerisis damping model
%fm=flim/2;
%Cs=eta.*Ks/fm/2/pi;
%fixed-fixed
if boundary==1
    Ms(1,:)=[];Ms(:,1)=[];Ms(end,:)=[];Ms(:,end)=[];
    Ks(1,:)=[];Ks(:,1)=[];Ks(end,:)=[];Ks(:,end)=[];
    Cs(1,:)=[];Cs(:,1)=[];Cs(end,:)=[];Cs(:,end)=[];
    Gs(1,:)=[];Gs(:,1)=[];Gs(end,:)=[];Gs(:,end)=[]; 
    ndof_n=ndof-2;
    Gs(1,:)=zeros(1,ndof_n);
% WABC
elseif boundary==2
    Ms(:,1)=Ms(:,1)+Ms(:,end);
    Ms(1,:)=Ms(1,:)+Ms(end,:);
    Ms(:,end)=[];Ms(end,:)=[]; 
    Ks(:,1)=Ks(:,1)+Ks(:,end);  
    Ks(1,:)=Ks(1,:)+Ks(end,:);
    Ks(:,end)=[];Ks(end,:)=[];
    Cs(:,1)=Cs(:,1)+Cs(:,end);  
    Cs(1,:)=Cs(1,:)+Cs(end,:);
    Cs(:,end)=[];Cs(end,:)=[];
    ndof_n=ndof-1;
end
%%
gamma_c=[0*gp 0*gi 0*gd];
[~,~,so_passive,~,sys_wz_passive,~] = function_LumpedEigenfreq(ndof,m1,m2,k,b,flim,boundary,e,gamma_c,theta,phi);
%% Assemblyng dynamic stiffness matrix and solving u(\omega)  
W=length(wv);
Dg=zeros(ndof_n,ndof_n,W);
 %% define where the external load is applied
 % F=1 for FRF
 F=1;
% force in the middle
if e == 'm'
    m = floor(ndof_n/2);
    if mod(ndof_n,2) == 0
         F = [zeros(ndof_n-m-1,1);F;zeros(ndof_n-m,1)];
    else
        F = [zeros(m,1);F;zeros(m,1)];
    end
% force on the left end
elseif e == 'l'
    F = [F; zeros(ndof_n-1,1)];
% force on the right end    
elseif e=='r'
    F = [zeros(ndof_n-1,1);F];
else 
    disp('invalid force position')
    return
end
q=zeros(ndof_n,W);
% This loop runs through the frequency vector
% Computing the FRF through the dynamic stiffness matrix
%     for i=1:W
%        w=wv(i);
%        Dm=-w^2*Ms-1i*w*Cs+Ks+Gs;
%       q(:,i)=Dm\F;
%     end
% q1m=abs(q(1,:));
% qNm=abs(q(end,:));
% q1p=atan2d(imag(q(1,:)),real(q(1,:)));
%qNp=atan2d(imag(q(end,:)),real(q(end,:)));
%end
 % or from the state space 4 matrix representation
 nw=1;nz=ndof_n;
 A=[zeros(ndof_n) eye(ndof_n)
      -Ms\(Ks+Gs)    -Ms\Cs];
B=[zeros(ndof_n,1)
        Ms\F];
C=[eye(ndof_n)  zeros(ndof_n)];
D=zeros(nz,nw);
sys_simo=ss(A,B,C,D);
P=tf(sys_simo);
P1=P(1,1);Pn=P(ndof_n,1);
% BODE diagram of tf for output u(1)
[mag,phase]=bode(P1,wv);
q1m=vec(mag);
q1p=vec(phase);
% BODE diagram of tf for output u(ndof)
[mag,phase]=bode(Pn,wv);
qNm=vec(mag);
qNp=vec(phase);

 nz2=1;
 C=[ 1 zeros(1,ndof_n-2) -1 zeros(1,ndof_n)];
 D=zeros(nz2,nw);
 sys_siso=ss(A,B,C,D);
 %output1= function_stability_c(A,eye(size(A,1)))
% output2=function_H2Norm(A,B,C,D)
%%
figure(3)
plot(fv,20*log10(q1m),'r','MarkerSize',8,'LineWidth',1)
hold on
plot(fv,20*log10(qNm),'k','MarkerSize',8,'LineWidth',1)
hold off
ylabel('$ |\hat{q} / \hat{F}| (i \omega)$ [dB]', 'interpreter', 'latex', 'fontsize', 15)
xlabel('$f$ [Hz]', 'interpreter', 'latex', 'fontsize', 15)
set(gca,'TickLabelInterpreter','Latex','fontsize',15);
legend('L','R')
xlim([0 flim])
box on
set(gcf, 'Color', 'w');
set(gca,'TickLabelInterpreter','Latex','fontsize',15);
grid on
figure(4)
plot(fv,q1p,'r','MarkerSize',8,'LineWidth',1)
hold on
plot(fv,qNp,'k','MarkerSize',8,'LineWidth',1)
hold off
ylabel('$ \hat{q} / \hat{F} (i \omega) $', 'interpreter', 'latex', 'fontsize', 15)
xlabel('$f$ [Hz]', 'interpreter', 'latex', 'fontsize', 15)
set(gca,'TickLabelInterpreter','Latex','fontsize',15);
legend('L','R')
xlim([0 flim])
box on
set(gcf, 'Color', 'w');
set(gca,'TickLabelInterpreter','Latex','fontsize',15);
grid on
%%
if FRF3D==1
    % Computing the FRF through the dynamic stiffness matrix
    for i=1:W
       w=wv(i);
       Dm=-w^2*Ms-1i*w*Cs+Ks+Gs;
      q(:,i)=Dm\F;
    end
    figure(5)
    xs=1:1:ndof;
    if boundary==1
        q=[zeros(1,W); q; zeros(1,W)];
    elseif boundary==2
        xs(:,end)=[];
    end
    mesh(fv,xs,20*log10(abs(q)))
    %contour(fv(10:end)/1000,xs,20*log10(abs(u_M_SEM)))
    %view([0 -90])
    shading INTERP
    colormap jet
    h=colorbar;
    ylabel(h,'$| \hat{u}(x) |$ [dB]','interpreter', 'latex','FontSize', 15)
    %w_real_TM_4PB(end)/2/pi/1000 =  919.9540
    ylabel('$x$ [m]', 'interpreter', 'latex', 'fontsize', 15)
    xlabel('$\Re(f)$ [kHz]', 'interpreter', 'latex', 'fontsize', 15)
    %zlabel('$| \hat{u}(x) | [dB]$', 'interpreter', 'latex', 'fontsize', 15)
    set(gca,'TickLabelInterpreter','Latex','fontsize',15);
    box on
    set(gcf, 'Color', 'w');
end
%% s=-i \omega
%((-i\omega)^2 M -i\omega C +K+G)q(\omega) = 0
if b ~=0
    [V,so]=polyeig(Ks+Gs,Cs,Ms);
else
    [V,L]=eig(Ks+Gs,-Ms);
    so=[sqrt(diag(L));-sqrt(diag(L))];
end
omega_n=so.*1i;
[omega_n,ind]=sort(omega_n,'ComparisonMethod','real');
%[V,L]=eig(Ks,-Ms);
%so_h=[sqrt(diag(L));-sqrt(diag(L))];
%%
%     figure(6)
%     scatter(real(so),imag(so),'kx','LineWidth',1.5)
%     hold on
%     scatter(real(so_h),imag(so_h),'bx','LineWidth',1.5)
%     ylim([0 8])
%     grid on
%     title('Eigenvalues', 'interpreter', 'latex', 'fontsize', 15)
%     ylabel('$\Im ( s )$', 'interpreter', 'latex', 'fontsize', 15)
%     xlabel('$\Re ( s )$', 'interpreter', 'latex', 'fontsize', 15')
%     set(gcf, 'Color', 'w');
%     set(gca,'TickLabelInterpreter','Latex','fontsize',20);
%     box on
    %%
if dispersion==1
    figure(7)
    scatter(real(so),imag(so),'kx','LineWidth',1.5)
    hold on
    plot(real(-1i.*omega_a),imag(-1i.*omega_a),'r','LineWidth',1.5)
    plot(real(1i.*omega_a),imag(1i.*omega_a),'r','LineWidth',1.5)

  %  scatter(real(so_h),imag(so_h),'bx','LineWidth',1.5)
    scatter(real(so_wabc),imag(so_wabc),'ro')
    scatter(real(so_passive),imag(so_passive),'bx','LineWidth',1.5)
    ylim([0 11])
    grid on
    %title('Eigenvalues', 'interpreter', 'latex', 'fontsize', 15)
    ylabel('$\Im ( s )$', 'interpreter', 'latex', 'fontsize', 15)
    xlabel('$\Re ( s )$', 'interpreter', 'latex', 'fontsize', 15')
    set(gcf, 'Color', 'w');
    set(gca,'TickLabelInterpreter','Latex','fontsize',20);
    box on
end
%% Eigenmodes animation
x=1:1:ndof_n;
if eigenmode==1
   % ylim([-20 0])
    wnv=-1i*so;
        [fi_imag,fi_real]=ginput(1);
        dfi=1/2;
        j=1;
        % searching the eigenfrequency by its real and imaginary parts
        ind_real=find(abs(real(wnv)-fi_real)<=dfi);
        ind_imag=find(abs(imag(wnv)-fi_imag)<=dfi);
        ind=intersect(ind_real,ind_imag);
        dfi=dfi/10;
        while length(ind) ~= 1 && j <=100
            ind_real=find(abs(real(wnv)-fi_real)<=dfi);
            ind_imag=find(abs(imag(wnv)-fi_imag)<=dfi);
            ind=intersect(ind_real,ind_imag);
            % updates accuracy until the only eigenfrequency find is the selected
            if isempty(ind)
                dfi=2*dfi;
            else
                dfi=dfi/10;
            end
    
            j=j+1;
        end
        % norm1 = max(abs(real(Vn(:,mode1))));
        disp('ind [rad/s]')
        wn=wnv(ind);
        display(wn)
       
     
        %modo=sin(x)+1j*cos(x)/2;
        mode=V(:,ind);
        %wn=-1i*so(end);
        norm = max(abs(mode));
        % numero de frames
        Nf=100;
        % numero de ciclos
        nc=1; 
        teta=0:2*nc*pi/(Nf-1):2*nc*pi; 
        figure
        % l=length(x);
        % x = [x x(l).*ones(1,l-1)+x(2:l)];
        for i=1:Nf
            y=abs(mode).*cos(angle(mode)+teta(i));
            scatter(x(1:ndof_n),y'/norm)
        %    y=[y;y];
        %    plot(x(1:end-1),y'/norm)
            axis([min(x),max(x), -1,1]);
            xlabel('$n$', 'interpreter', 'latex', 'fontsize', 15)
            ylabel('$ \hat{g}(\omega_n) $', 'interpreter', 'latex', 'fontsize', 15)     
            set(gca,'TickLabelInterpreter','Latex','fontsize',15);
            box on
            set(gcf, 'Color', 'w');
            legend(['\omega_n = ' num2str(wn,'%0.3f') 'rad/s']);
            drawnow
           
        end
end
%% zero state response of the closed loop for input w
% wc=0.5;
% fc=wc/2/pi;
% Nc=30; 
% Tc=1/fc;
% T2=Nc*Tc*2; 
% f2=1/T2; w2=2*pi*f2;
% dt=Tc/100;
% Nt=8*9046;
% t=0:dt:(Nt-1)*dt;
% tend=round(length(t)/10);
% 
% w=double(sin(wc*t).*sin(w2*t).^2.*heaviside(T2/2-t)); 
% df=1/t(end); fm=1/dt;
% fvec=0:df:fm/2;
% 
% Fw=1/Nt*fft(w);
% % 
% figure
% subplot(1,2,1)
% plot(t(1:tend),w(1:tend))
% xlabel('$t$ [s]', 'interpreter', 'latex','fontsize', 15)
% %xaxis([0 100])
% box on
% grid on
% title('External Force Dynamics', 'interpreter', 'latex','fontsize', 15)
% subplot(1,2,2)
% plot(2*pi*fvec,abs(Fw(1:Nt/2)))
% xlabel('$\omega $ [rad/s]', 'interpreter', 'latex','fontsize', 15)
% xlim([0 2])
% title('External Force Spectrum', 'interpreter', 'latex','fontsize', 15)
% box on
% grid on
% set(gcf, 'Color', 'w');
% %% zero state response of the closed loop for input w
% 
% ni=size(A,1);
% z=lsim(sys_wz,w,t,zeros(ni,1));
% %%
% figure
% Nplots=100;
% DI=round(tend/Nplots);
% %i=tend;
% for i=1:10:tend
%     scatter3(x(1:size(z,2)),(t(i)*ones(1,size(z,2))),z(i,:),'b');
%      hold on
% end
% xlabel('$n$', 'interpreter', 'latex', 'fontsize', 15)
% ylabel('$t$ [s]', 'interpreter', 'latex', 'fontsize', 15)
% zlabel('$q_n(t)$ [m]', 'interpreter', 'latex', 'fontsize', 15)
% hold off
% box on
% grid on
% set(gcf, 'Color', 'w');
% %%
% figure
% %tend=length(t);
% plot(t(1:tend),z(1:tend,1),'b')
% hold on
% plot(t(1:tend),z(1:tend,end),'k')
% set(gcf, 'Color', 'w');
% xlim([0,t(tend)])
% xlabel('$t [s]$','interpreter', 'latex','FontSize',15');
% legend('$q_1(t) [m]$','$q_N(t) [m]$','interpreter', 'latex','FontSize',15');
% grid on

%% Stabilizable dynamic output controller synthetis
%options.hinf=1e+20;
%output = function_Hinf_Dynamic_OutputFeedback(Ass,B1ss,B2ss,C1ss,C2ss,D11ss,D12ss,D21ss)
%% Internal stalibility test
%  output.feas=0 = > internally unstable
% output.feas=1 = > internally stable
% output = function_stability_c(Ass)
% output = function_Hinf_Optimal_Controller_c(Ass,B1ss,B1ss,C1ss,D11ss,D12ss)
% %%
% function_stabilizability_c(Ass,B2ss)
%  %%
%  % sample frequency for x domain
% fxsamp=1;% The distance between each mass is 0.5 for the dimer and 1 for monoatomic model 
% %(unitary distance is assumed as the length of the cell)
% % total length of x domain
% Nx = length(x);
% % resolution of spacial frequency (kx) mapped from x
% dkx = fxsamp/Nx;
% % Nysquist wavenumber
% kxN = fxsamp/2;
% % if Nx is odd
% if mod(Nx,2) ~= 0
% %    kx=-kxN+dkx/2:dkx:kxN-dkx/2;
%     kx = (1/dkx)*(-(Nx-1)/2:(Nx-1)/2)/Nx;
% else
%     kx = -kxN:dkx:kxN-dkx;
% end
% kx=2*pi*kx;
% p=0:.1:1;
% for j=1:size(Vn,2)
%     Vn_til(:,j)=1/Nx*fft(Vn(:,j));
% end
% 
% scatter(imag(wn),real(wn),'k')
% hold on
% normm = max(max(abs(Vn_til)));
% contour(wn,kx/pi,abs(Vn_til)./normm,p)
% set(gca,'TickLabelInterpreter','Latex','fontsize',15);
% box on
% grid on
% set(gcf, 'Color', 'w');
% colormap('jet')
% shading interp
% %ylim([0 2])
%echo off