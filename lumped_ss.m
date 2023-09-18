%% clear 
clear all 
%close all 
clc 
%% lumped-parameter system set up
k=1; m1=1; m2=m1;
%k2=1; m2=2;
%r1=0; r2=0; %range of non-reciprocal coupling, r=0 means next-neighborhoud
%couple
r=0;
eta=0;
e='m';
% boundary == 0 : free-free (default)
% bundary == 1: fixed-fixed
% boundary == 2 : PBC
boundary=2;
% angular frequency vector limit
flim=3.5e3;
wv=2*pi*(2:2:flim);
fv=wv/2/pi;
Nd=100*2*pi*flim;
% gains
kp=-1;
ki=kp;
kd=kp;
gamma_c=[1*kp 0*ki 0*kd];
s=tf('s');
% qusiperiodicity parameters
theta=0;phi=0;
%% FEM
ncell=10;
ne_cell=1;
output = function_buildFEM_xl(m1,m2,k, eta ,ne_cell,ncell,flim,'boundary',boundary);
%M=output.M;C=output.C;K=output.K;
Ms=output.Ms;Cs=output.Cs;Ks=output.Ks;
x=output.x;ndof=output.ndof;%% State Space of interconnection (S,H) 
%% S
% number of inputs=outputs
ny=ncell*ne_cell-r-1;
output=function_buildSS_xl(Ms,Cs,Ks,ndof,e,ny);
Ass=output.Ass;nx=size(Ass,1);
% PS: (Ms,Cs,Ks) were previously multiplied by B*rho,
% so we need to multiply matrixes F (inside B1ss) and T(inside B2ss) also 
B1ss=output.B1ss;nw=size(B1ss,2);
B2ss=output.B2ss;nu=size(B2ss,2);
C1ss=output.C1ss;nz=size(C1ss,1);
C2ss=output.C2ss;ny=size(C2ss,1);
D11ss=output.D11ss;
D12ss=output.D12ss;
D21ss=output.D21ss;

%% H 
% Feedback law in pressure by volume velocity
H_v=gamma_c(1)+gamma_c(2)/s+gamma_c(3)*Nd/(1+Nd/s);

%% setting analog filters
% RC filter
% wf=1e3*2*pi;
% Hf = s/wf/(1+s/wf);
% no filter
Hf=1;
%butterworth filter
% w1=1800*2*pi;
% w2=2000*2*pi;
% % lp= 'low pass', bp = 'band pass'
% type='hp';
% out=function_filter(-3,-10,w1,w2,type);
% Hf=out.Hf;
%Hf=Filtro_ver2_completo_v1;
% final feedback transfer function
H=H_v*Hf;
%%  closed loop (S,H)
if  norm(gamma_c,1)~=0
%     figure
%     bode(H)
%     title('bode diagram of the feedback law')
    [num,dem]=tfdata(H,'v');
    if size(num,1) > size(dem,1)
        disp('improper feedback law\n');
    end
    %% controlable state space realization of H
    [Ac1,Bc1,Cc1,Dc1]=tf2ss(num,dem);
    nxc=size(Ac1,1);
    %nu=ny
    Acss=zeros(ny*nxc,ny*nxc);
    Bcss=zeros(ny*nxc,ny);
    Ccss=zeros(ny,ny*nxc);
    Dcss=zeros(ny,ny);
    %% building state space realization of H (Acss,Bcss,Ccss,Dcss)
    % by blocks, since each cell has a decoupled feedback law
    for k=1:ny
        lines=(k-1)*nxc+1:k*nxc;
        Acss(lines,lines)=Ac1;
        Bcss(lines,k)=Bc1;
        Ccss(k,lines)=Cc1.*cos(theta*k+phi);
        Dcss(k,k)=Dc1.*cos(theta*k+phi);
    end
    %Acss=Ac1*eye(ncell);   Bcss=Bc1*eye(ncell);
    %Ccss=Cc1*eye(ncell);  
    %Dcss=Dc1*eye(ny);
    %% closed loop (S,H) interconnection
    Ass_cl=[Ass+B2ss*Dcss*C2ss B2ss*Ccss
            Bcss*C2ss           Acss];
    Bss_cl=[B1ss+B2ss*Dcss*D21ss
                Bcss*D21ss];
    Css_cl=[C1ss+D12ss*Dcss*C2ss     D12ss*Ccss];
    Dss_cl=D11ss+D12ss*Dcss*D21ss;
    sys_wz = ss(Ass_cl,Bss_cl,Css_cl,Dss_cl);
    Css=Css_cl;
    nx_cl=size(Ass_cl,1);
else 
    % only logical case other then the previous is gamma_c==0, meaning open
    % loop system, u =0
    sys_wz = ss(Ass,B1ss,C1ss,D11ss);
    Css=C1ss;
end
%% verify internal stability of closed loop system
if  norm(gamma_c,1)~=0
    [V, Lambda] = eig(Ass_cl);
else
    [V, Lambda] = eig(Ass);
end
Lambda = diag(Lambda);
Vn = Css*V;
wn = -1i*Lambda;
[wn,ind] = sort(wn,'ComparisonMethod','real');
Vn = Vn(:,ind); 

%%
 figure
scatter(imag(wn),real(wn),'k')
grid on
title('Eigenfrequencies', 'interpreter', 'latex', 'fontsize', 15)
xlabel('$\Im \{ \omega_n \}$', 'interpreter', 'latex', 'fontsize', 15)
ylabel('$\Re \{ \omega_n \}$', 'interpreter', 'latex', 'fontsize', 15')
set(gcf, 'Color', 'w');
box on
%% zero state response of the closed loop for input w
wc=0.5;
fc=wc/2/pi;
Nc=30; 
Tc=1/fc;
T2=Nc*Tc*2; 
f2=1/T2; w2=2*pi*f2;
dt=Tc/100;
Nt=8*9046;
t=0:dt:(Nt-1)*dt;

F=double(sin(wc*t).*sin(w2*t).^2.*heaviside(T2/2-t)); 
df=1/t(end); fm=1/dt;
fvec=0:df:fm/2;

Fw=1/Nt*fft(F);
% 
% figure
% subplot(1,2,1)
% plot(t,F)
% title('external force dynamics')
% subplot(1,2,2)
% plot(2*pi*fvec,abs(Fw(1:Nt/2)))
% xlim([0 2])
% title('external force spectrum')
%% zero state response of the closed loop for input w
if norm(gamma_c,1)~=0
    ni=nx_cl;
else
    ni=nx;
end
z=lsim(sys_wz,F,t,zeros(ni,1));
%%
figure
Nplots=100;
tend=length(t)/10;
DI=round(tend/Nplots);
for i=1:DI:tend
    plot3(x(1:size(z,2)),t(i)*ones(1,length(x)),z(i,:),'b');
     hold on
end
xlabel('$x$ [m]', 'interpreter', 'latex', 'fontsize', 15)
ylabel('$t$ [s]', 'interpreter', 'latex', 'fontsize', 15)
zlabel('$u (x,t)$ [m]', 'interpreter', 'latex', 'fontsize', 15)
hold off
box on
grid on
set(gcf, 'Color', 'w');
%%
figure
for i=1:1:tend
    normm = max(abs(z(i,:)));
    z(i,:)=abs(z(i,:))/normm;
end
mesh(t(100:tend).*1e3,x(1:size(z,2)),z(100:tend,:).');
set(gcf, 'Color', 'w');
colormap jet
h=colorbar ;
ylabel(h,'$||P_n(t)||_{\infty}$','interpreter', 'latex','FontSize',15)
hold on
ylabel('$x [m]$','interpreter', 'latex','FontSize',15); 
xlabel('$t [ms]$','interpreter', 'latex','FontSize',15');
zlabel('$|p| [Pa]$','interpreter', 'latex','FontSize',15);
view([0,90])
tlim=t(end).*1e3;
%xlim([0 tlim])
box on
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