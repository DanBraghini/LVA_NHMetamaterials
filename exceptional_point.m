%% clear 
clear all 
close all 
clc 
%% Macros
%feedback = input('enter the type of feddback\n(string: 
%0 = integral, 1=derivative, 2 = proportional, 3 = double integral)\n','s');
plotfilter=0;
filter=0;
feedback=0;%PC properties
if feedback==0
    string_par='\gamma_I';
    kp=1.5e-2;
elseif feedback==1
    string_par='\gamma_D';
    kp=2e-9;
elseif feedback==2
    string_par='\gamma_P';
    kp=1e-7;
elseif feedback==3
     string_par='\gamma_{II}';
    kp=1000000e-7;
else
    disp('Please choose other kind of feedback or improve this code.\n')
end
    %dmodel = input('enter the selected damping model\n(v = viscous, s = structural)')
dmodel='s';
% filter type
ideal_filter=0;
% angular frequency vector
flim=1.5e3;%in Hz
wv=2*pi*(2:2:flim);
%wv=2*pi*(-flim:2:flim);
fv=wv/2/pi;
Nd=100*2*pi*flim;
s=tf('s');
e='m';
% locality parameter
a=0;
% boundary == 0 : closed-closed (default)
% boundary == 1 : open-open
% boundary == 2 : periodic(approaches infinity system as ncell gets higher)
boundary=0;
% ne_ cell: number of finite elements per cell
%(multiple of 3 to be divided between the 3 segments)
ne_cell=3*3;
% ncell:number of cells that make the structure
ncell=18;
modulation = 'alpha';
% qusiperiodicity parameters
if modulation == 'phiii'
    theta=1/3;phi=(0:1:ncell*2*pi)./ncell;
    parm=phi;
elseif modulation == 'theta'
     phi=0;%theta=0:.1:2*pi;
     theta=(0:1:ncell*2*pi)./ncell;
     parm=theta;
elseif modulation == 'alpha'
     phi=0;theta=0;
     kp=linspace(-kp,kp, 100);
     parm=kp;
     kp=kp';
end

%% acoustic metamaterial set up
if dmodel == 's'
    eta=0.001;
else
   eta= 100;
end
eta=0;
Lc = 50e-2;
rho=1.225;
c=343;
B=rho*c^2;
r1 = 2e-2;
A1 = pi*r1^2;
r2 = r1;
A2 = A1;
xs = Lc/4;
xi = xs+Lc/2;
L1= xs;
L2 = xi- xs;
%% setting analog filters
% no filter
if filter == 0
    Hf=1;
elseif filter == 1% butterworth filter
    out=function_filter(-3,-10,filterinfo.w1,filterinfo.w2,filterinfo.type);
    Hf=out.Hf;
elseif filter == 2% RC filter 
    wf=filterinfo.wf;
    Hf = s/wf/(1+s/wf);
elseif filter == 3
%------ filter built for the experiment---------------
     Hf=Filtro_ver2_completo_v1;
end
if plotfilter==1
    s1=1i*wv;
    FRFf=zeros(1,length(wv));
    for j=1:length(wv)
         [num,dem]=tfdata(Hf,'v');
         num=real(num);
         dem=real(dem);
         nn=size(num,2);
         nd=size(dem,2);
         Ns=0;Ds=0;
         for i=1:nn 
            Ns=Ns+num(i)*s1(j)^(nn-i);
         end
         for i=1:nd
            Ds=Ds+dem(i)*s1(j)^(nd-i);
         end
         FRFf(j)=Ns/Ds;
    end
    figure
    plot(wv/2/pi/1000,20*log10(abs(FRFf)))
     xlabel('$f$ [kHz]', 'interpreter', 'latex', 'fontsize', 15)
     ylabel('$[dB]$', 'interpreter', 'latex', 'fontsize', 15)
     figure
    plot(wv/2/pi/1000,180*angle(FRFf)/pi)
     xlabel('$f$ [kHz]', 'interpreter', 'latex', 'fontsize', 15)
     ylabel('$degrees$', 'interpreter', 'latex', 'fontsize', 15)
    box on
    set(gcf, 'Color', 'w');
    grid on
end

%% setting feedback law
 aux=zeros(size(kp)) ;
if feedback==0
    gamma_c=[kp aux aux aux  ];
elseif feedback==1
     gamma_c=[ aux kp aux aux ];
elseif feedback==2
     gamma_c=[aux aux kp aux];
elseif feedback==3
     gamma_c=[aux aux aux kp];
end
% Feedback law in pressure by volume velocity
Nd=100*2*pi*flim;
H=gamma_c(:,1)/s+...
    gamma_c(:,2)*Nd/(1+Nd/s)...
    +gamma_c(:,3)...
     +gamma_c(:,4)/s/s;
H_v=H.*Hf;
%% Correction to volume acceleration
H=H_v.*(Nd/(1+Nd/s));
%% FEM
output = function_buildFEM(rho, L1, L2, Lc, A1, A2, c, eta,ne_cell,ncell,flim,xs,xi...
    ,'dmodel',dmodel,'boundary',boundary);
%M=output.M;C=output.C;K=output.K;
Ms=output.Ms;Cs=output.Cs;Ks=output.Ks;
x=output.x;ndof=output.ndof;
n1=output.n1;n2=output.n2;
%% State Space of interconnection (S,H) 
%% S
output=function_buildSS(Ms,Cs,Ks,ndof,ncell,n1,n2,e,a);
Ass=output.Ass;nx=size(Ass,1);
% PS: (Ms,Cs,Ks) were previously multiplied by A*rho on FEM, following 
% the formulation from Neurodinne&Atalla, section 4.2. The initial 
% dimension of external inputs to the PDE was [B/A*(volume acceleration)].
% So we need to multiply matrixes F (inside B1ss) and T(inside B2ss) by
% (B*rho) so that unitary generalized forces correspond to volume
% accelerationB1ss=(B*rho).*output.B1ss;nw=size(B1ss,2);
B1ss=(B*rho).*output.B1ss;nw=size(B1ss,2);
B2ss=(B*rho).*output.B2ss;nu=size(B2ss,2);
C1ss=output.C1ss;nz=size(C1ss,1);
C2ss=output.C2ss;ny=size(C2ss,1);
D11ss=output.D11ss;
D12ss=output.D12ss;
D21ss=output.D21ss;
%%
np=length(parm);
wo=cell(size(parm));
Vo=wo;fo=wo;
%wo=zeros(nxcl,np);
%Vo=zeros(length(x),nxcl,np);
for i=1:np
   if modulation == 'phiii'
       phim=phi(i);thetam=theta;Hm=H;
   elseif modulation == 'theta'
       phim=phi;thetam=theta(i);Hm=H;
   else
       phim=0;thetam=0;Hm=H(i);
   end
%    if norm(Hm) == 0
%         continue
%    end
    
   output = function_eigenspectrum_varyingp(Hm, Ass,B2ss,C1ss,C2ss,D12ss,thetam,phim);
   % storing eigenfrequencies in rad/s
   wo{i}=output.wo;
   % storing eigenfrequencies in kHz
   fo{i}=wo{i}./(2*pi*1000);
   % storing eigenmodes
   Vo{i}=output.Vo;
end
%% 
nx_cl=output.nx_cl;
%% effective electric gains
% Transducers: 
% T1 = sensor, T2 = actuator
T1 = function_getMicTransd;
T2 =  function_getLoudSpeakerTransd;
%T1=1;T2=T1;
%H_v=T2.*H_e*T1;
parm=parm./(T1*T2);

%% plots

figure
for i=1:np
    plot(parm(i).*ones(1,nx_cl),real(fo{i}),'.k')
    hold on
end
grid on
title('Modulated spectrum', 'interpreter', 'latex', 'fontsize', 15)
ylabel('$\Re(f) [kHz]$', 'interpreter', 'latex', 'fontsize', 15)
xlabel(['$' string_par '$'], 'interpreter', 'latex', 'fontsize', 15')
set(gcf, 'Color', 'w');
box on
ylim([0 flim/1000])
xlim([parm(1) parm(end)])
%xlim([])
%% 
figure
for i=1:np
    plot(parm(i).*ones(1,nx_cl),imag(fo{i}),'.k')
    hold on
end
grid on
title('Modulated spectrum', 'interpreter', 'latex', 'fontsize', 15)
ylabel('$\Im(f) [kHz]$', 'interpreter', 'latex', 'fontsize', 15)
xlabel(['$' string_par '$'], 'interpreter', 'latex', 'fontsize', 15')
set(gcf, 'Color', 'w');
box on
ylim([-0.25 0.25])
xlim([parm(1) parm(end)])
%%
% 
% figure
% for i=1:np
%     plot(imag(wo(:,i))/2/pi/1000,real(wo(:,i))/2/pi/1000,'.k')
%     hold on
% end
% %plot(real(wo(1:nx,:)),phi/pi)
% grid on
% title('Modulated spectrum', 'interpreter', 'latex', 'fontsize', 15)
% xlabel('$\Im \{ \omega \}$', 'interpreter', 'latex', 'fontsize', 15)
% ylabel('$\Re  \{ \omega \}$', 'interpreter', 'latex', 'fontsize', 15')
% set(gcf, 'Color', 'w');
% box on
%%
 figure
for i=1:np
    plot3(parm(i).*ones(1,nx_cl),imag(fo{i}),real(fo{i}),'.k')
    hold on
end
%%
grid on
title('Modulated spectrum', 'interpreter', 'latex', 'fontsize', 15)
xlabel(['$' string_par '$'], 'interpreter', 'latex', 'fontsize', 15')
ylabel('$\Im(f) [kHz]$', 'interpreter', 'latex', 'fontsize', 15)
zlabel('$\Re(f) [kHz]$', 'interpreter', 'latex', 'fontsize', 15')
zlim([0 flim/1000])
ylim([-.5 .5])
set(gcf, 'Color', 'w');
box on
%%
    v1=2e-10/T1/T2;v2=2;v3=flim/1000;
    fh1 = fill3([v1 v1 v1 v1],[-v2  -v2  v2 v2],[0 v3 v3 0],'r');
    set(fh1,'LineWidth',1.0);
    set(fh1,'EdgeColor','none');
    set(fh1,'facealpha',.2);