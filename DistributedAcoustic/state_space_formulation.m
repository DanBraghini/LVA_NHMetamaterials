%% clear
clear all
clc
%close all
%% Created on : 17/05/2021 by Danilo Braghini
%  Last update: 29/08/2021 by Danilo Braghini
%% macros
s=tf('s');
%feedback = input('enter the type of feddback\n(string: d = derivative, i = integrative and p = proportional)\n','s');
feedback='pid';
%dmodel = input('enter the selected damping model\n(v = viscous, s = structural)')
dmodel='s';
e='m';
% angular frequency vector limit
flim=3.5e3;
Nd=100*2*pi*flim;
%PS: both structural damping model and derivative feedback aproximation
%depends on flim
%%  acoustic metamaterial set up
%PC properties
% acoustic metamaterial
if dmodel == 's'
    eta=.01;
elseif dmodel == 'v'
   eta= 100;
elseif dmodel == 'n' 
    eta=0;
else
   disp('invalid damping model\n choose structural (s),viscous (v) or no damping (n).') ;
end
Lc = 50e-2;
r1 = 2e-2;
A1 = pi*r1^2;
r2 = r1;
A2 = A1;
xs = Lc/4;
xi = xs+Lc/2;
rho =1.225;
c=343;
B=rho*c^2;
L1= xs;
L2 = xi- xs;
kp=1e-5;
ki=1.5e-3;
kd=2e-10;
%% Feedback control law and sensor locality paramater 
if feedback=='pid'
    gamma_c=[1*kp 0*ki 0*kd];
elseif feedback=='l-l'
    % lead-lag
    z=2*pi*1000;
    p=2*pi*1;
    gamma_c=1e-4*ki;
elseif feedback=='new'
    % double pole
    gamma_c=2e3*ki;
    p1=2*pi*1e3;
    p2=p1;
end
a=0;

%% FEM
% ncell:number of cells that make the structure (multiple of 3 to be divided between the 3 segments)
% number of elements per cell
ne_cell=7*3;
%number of cells on the structure
ncell=18;
% boundary == 0 : free-free (default)
% bundary == 1: fixed-fixed
output = function_buildFEM(rho, L1, L2, Lc, A1, A2, c, eta,ne_cell,ncell,flim,xs,xi,'dmodel',dmodel);
M=output.M;C=output.C;K=output.K;
Ms=output.Ms;Cs=output.Cs;Ks=output.Ks;
x=output.x;ndof=output.ndof;
n1=output.n1;n2=output.n2;
%% State Space of interconnection (S,H) 
%% S
output=function_buildSS(Ms,Cs,Ks,ndof,ncell,n1,n2,e,a);
Ass=output.Ass;nx=size(Ass,1);
% PS: (Ms,Cs,Ks) were previously multiplied by B*rho,
% so we need to multiply matrixes F (inside B1ss) and T(inside B2ss) also 
B1ss=(B*rho).*output.B1ss;nw=size(B1ss,2);
B2ss=(B*rho).*output.B2ss;nu=size(B2ss,2);
C1ss=output.C1ss;nz=size(C1ss,1);
C2ss=output.C2ss;ny=size(C2ss,1);
D11ss=output.D11ss;
D12ss=output.D12ss;
D21ss=output.D21ss;
%% H 
% Feedback law in pressure by volume velocity
if feedback=='pid'
    H_v=gamma_c(1)+gamma_c(2)/s+gamma_c(3)*Nd/(1+Nd/s);
elseif feedback=='l-l'
    H_v=gamma_c*(s+z)/(s+p);
 elseif feedback=='new'
    H_v=gamma_c/((s+p1)*(s+p2));
end
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
Hf=Filtro_ver2_completo_v1;
% final feedback transfer function
H_v=H_v*Hf;
%% Correction to volume acceleration
H=H_v*Nd/(1+Nd/s);
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
    Acss=zeros((ncell-a)*nxc,(ncell-a)*nxc);
    Bcss=zeros((ncell-a)*nxc,ncell-a);
    Ccss=zeros(ncell-a,(ncell-a)*nxc);
    %% building state space realization of H (Acss,Bcss,Ccss,Dcss)
    % by blocks, since each cell has a decoupled feedback law
    for k=1:ncell-a
        lines=(k-1)*nxc+1:k*nxc;
        Acss(lines,lines)=Ac1;
        Bcss(lines,k)=Bc1;
        Ccss(k,lines)=Cc1;
    end
    %Acss=Ac1*eye(ncell);   Bcss=Bc1*eye(ncell);
    %Ccss=Cc1*eye(ncell);  
    Dcss=Dc1*eye(ncell-a);
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
% figure
% scatter(real(Lambda),imag(Lambda),'+r')
% grid on
% title('Eiganvalues of state space active model')
% xlabel('$\Re \{ s \}$', 'interpreter', 'latex', 'fontsize', 15)
% ylabel('$\Im \{ s \}$', 'interpreter', 'latex', 'fontsize', 15')
%%
l_max=max(real(Lambda))
figure
scatter(real(wn),imag(wn),'k')
grid on
title('Eigenfrequencies', 'interpreter', 'latex', 'fontsize', 15)
xlabel('$\Re \{ \omega_n \}$', 'interpreter', 'latex', 'fontsize', 15)
ylabel('$\Im \{ \omega_n \}$', 'interpreter', 'latex', 'fontsize', 15')
set(gcf, 'Color', 'w');
box on

%% Transient Response: applied force
% Nc is the number of circles from central frequency fc.
% T2 define the envelope frequency f2, which has Nc circles within.
% Nt is the number of entries on time vector t. The heavside function
% defines a window on half a period of the envelope.
% time discretization period dt depends on N, in such manner that
% the total time of the simulation T can be modified if nedeed.
Nt=4*4096*2;
% T=0.01; % total analysis time
Ap=1;
%fc=250;
fc=2500;
wc = 2*pi*fc;
Nc=15; 
Tc=1/fc;
T2=Nc*Tc*2; 
f2=1/T2; w2=2*pi*f2;
dt=T2/2/4096;
%dt = T2/Nt;
t=0:dt:(Nt-1)*dt;
w=Ap*sin(wc*t).*sin(w2*t).^2.*heaviside(T2/2-t); 

%% zero state response of the closed loop for input w
if norm(gamma_c,1)~=0
    ni=nx_cl;
else
    ni=nx;
end
z=lsim(sys_wz,w,t,zeros(ni,1));
%%
Nplots=100;
tend=length(t);
DI=round(tend/Nplots);
figure
for i=1:DI:tend
    plot3(x,(t(i).*10^3)*ones(1,length(x)),z(i,:),'b');
    hold on
end
xlabel('Length [m]'), ylabel('Time [ms]'), zlabel('Pressure p(x,t)[Pa]')
hold off
%ylim([0 100])
box on
set(gcf, 'Color', 'w');
%% 
% figure
% for i=1:1:tend
%     normm = max(abs(z(i,:)));
%     z(i,:)=abs(z(i,:))/normm;
% % end
% mesh(t(100:tend).*1e3,x,z(100:tend,:).');
% set(gcfDI, 'Color', 'w');
% colormap jet
% h=colorbar ;
% ylabel(h,'$||P_n(t)||_{\infty}$','interpreter', 'latex','FontSize',15)
% hold on
% ylabel('$x [m]$','interpreter', 'latex','FontSize',15); 
% xlabel('$t [ms]$','interpreter', 'latex','FontSize',15');
% zlabel('$|p| [Pa]$','interpreter', 'latex','FontSize',15);
% view([0,90])
% tlim=t(end).*1e3;
% xlim([0 tlim])
% box on