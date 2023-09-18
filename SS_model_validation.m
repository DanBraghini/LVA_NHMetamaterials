%% Created on : 2707/2021 by Danilo Braghini
clear all
clc
%close all
%% macros
% excitation point for FRF simulation
e='m';
flim=1e3;
wv=2*pi*(2:2:flim);
fv=wv/2/pi;
nf=length(fv);
%PS: both structural damping model and derivative feedback aproximation
%depends on flim
%dmodel = input('enter the selected damping model\n(v = viscous, s = structural)')
dmodel='s';
s=tf('s');
%% Crystal Set up
% acoustic metamaterial
if dmodel == 's'
    eta=0.001;
else
   eta= 10;
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
%% Feedback control law and sensor locality paramater 
kp=1e-7;
ki=1.5e-2;
kd=2e-10;
gamma_c= [0*kp ki 0*kd];
a=0;
%% 1)Adding PID to FEM matrix, using the standard physical model
% FEM
% ncell:number of cells that make the structure 
% ne: number of elements per cell (multiple of 3 to be divided between the 3 segments)
ne=3*3;
%number of cells on the structure
ncell=18;
outputFEM = function_buildFEM(rho, L1, L2, Lc, A1, A2, c, eta,ne,ncell,flim,xs,xi,'dmodel',dmodel);
M=outputFEM.M;C=outputFEM.C;K=outputFEM.K;
x=outputFEM.x;ndof=outputFEM.ndof;
n1=outputFEM.n1;n2=outputFEM.n2;
% Inserting local feedback interactions 
gamma_c_ph=(B*rho).*gamma_c;
%gamma_c_ph=gamma_c;
M(n2,n1) = M(n2,n1)-gamma_c_ph(3);
K(n2,n1) = K(n2,n1)-gamma_c_ph(2);
C(n2,n1) = C(n2,n1)-gamma_c_ph(1);
% joining cells to make a structure
Ms=zeros(ndof,ndof);Cs=Ms;Ks=Ms;

for i=1:ne:ndof-ne
    Ms(i:i+ne,i:i+ne) = Ms(i:i+ne,i:i+ne) + M;
    Cs(i:i+ne,i:i+ne) = Cs(i:i+ne,i:i+ne) + C;
    Ks(i:i+ne,i:i+ne) = Ks(i:i+ne,i:i+ne) + K;
end

%% 2) Using SS model
% State Space of (S,H) interconnection 
% S
Ms_2=outputFEM.Ms;Cs_2=outputFEM.Cs;Ks_2=outputFEM.Ks;
output=function_buildSS(Ms_2,Cs_2,Ks_2,ndof,ncell,n1,n2,e,0);
Ass=output.Ass;nx=size(Ass,1);
% PS: In this discretization, source is -B/A times volume velocity.
% Also, the matrixes (Ms,Cs,Ks) were previously multiplied by A*rho
B1ss=(B*rho).*output.B1ss;nw=size(B1ss,2);
B2ss=(B*rho).*output.B2ss;nu=size(B2ss,2);

C1ss=output.C1ss;nz=size(C1ss,1);
C2ss=output.C2ss;ny=size(C2ss,1);
D11ss=output.D11ss;
D12ss=output.D12ss;
D21ss=output.D21ss;
% H 
Nd=100*2*pi*flim;
% Feedback law in volume velocity
Hpid_v=gamma_c(1)+gamma_c(2)/s+gamma_c(3)*Nd/(1+Nd/s);
% Correction to volume acceleration
% Filter
%Hf=-10000*s/(s^2+10312.5*s+3125000);
%Hi=-0.702/(1.551*s+1);
%Hf = s^2 /(s^2 + 14140*s + 100000000);
Hf=1;
% final feedback transfer function
H=Hf*Hpid_v;
H=H*Nd/(1+Nd/s);
% figure
% bode(H_v);
% title('bode diagram of the feedback law')
[num,dem]=tfdata(H,'v');
[Ac1,Bc1,Cc1,Dc1]=tf2ss(num,dem);
nxc=size(Ac1,1);
Acss=zeros(ncell*nxc,ncell*nxc);
Bcss=zeros(ncell*nxc,ncell);
Ccss=zeros(ncell,ncell*nxc);
% building state space representation of H (Acss,Bcss,Ccss,Dcss)
% by blocks, since each cell has a decoupled controller
for k=1:ncell
    lines=(k-1)*nxc+1:k*nxc;
    Acss(lines,lines)=Ac1;
    Bcss(lines,k)=Bc1;
    Ccss(k,lines)=Cc1;
end
%Acss=Ac1*eye(ncell);   Bcss=Bc1*eye(ncell);
%Ccss=Cc1*eye(ncell);  
Dcss=Dc1*eye(ncell);
if norm(gamma_c,1)~=0
    %closed loop (S,H) interconnection
    Ass_cl=[Ass+B2ss*Dcss*C2ss B2ss*Ccss
            Bcss*C2ss           Acss];
    Bss_cl=[B1ss+B2ss*Dcss*D21ss
                Bcss*D21ss];
    Css_cl=[C1ss+D12ss*Dcss*C2ss     D12ss*Ccss];
    Dss_cl=D11ss+D12ss*Dcss*D21ss;
    nx_cl=size(Ass_cl,1);
    sys_wz = ss(Ass_cl,Bss_cl,Css_cl,Dss_cl);
else
    %open-loop system, u=0
    sys_wz = ss(Ass,B1ss,C1ss,D11ss);
end
% In case ki~=0 and kp==0 and kd==0, H output is alike a static compensator L
%     L = gamma_c(2).*eye(ncell);
%     %closed loop (S,L) interconnection
%     Ass_cl=Ass+B2ss*L*C2ss;
%     Bss_cl=B1ss+B2ss*L*D21ss;
%     Css_cl=C1ss+D12ss*L*C2ss;
%     Dss_cl=D11ss+D12ss*L*D21ss;
%     sys_wz = ss(Ass_cl,Bss_cl,Css_cl,Dss_cl);
%     nx_cl
%% Comparing models via Harmonic Forced Response 

%-------------------------------------------------------------------------%
%                            FRF via FEM                                      % 
%-------------------------------------------------------------------------%
u_FEM = zeros(ndof,nf);
% External force in unitary volume acceleration 
%F=B*rho;
% frequency loop
for i = 1:nf
    w = wv(i);
    % Dinamic stiffness matrix via FEM
    Dg = -w^2*Ms + 1i*w*Cs + Ks;
    % External force in constant volume velocity 
    % PS: In this discretization, source is -B/A times volume velocity.
    % Also, the matrixes (Ms,Cs,Ks) were previously multiplied by A*rho
    F=B*rho*1i*w;
    %F=1i*w;
    % force in the middle
    if e == 'm'
        m = floor(ndof/2);
        if mod(ndof,2) == 0
             Fv = [zeros(ndof-m-1,1);F;zeros(ndof-m,1)];
        else
            Fv = [zeros(m,1);F;zeros(m,1)];
        end
    % force on the left end
    elseif e == 'l'
        Fv = [F; zeros(ndof-1,1)];
    % force on the right end    
    elseif e=='r'
        Fv = [zeros(ndof-1,1);F];
    else 
        disp('invalid force position')
        return
    end
    u_FEM(:,i) = Dg\Fv;
end 
%%
%-------------------------------------------------------------------------%
%                 Numerical FRF via SS                                    % 
%-------------------------------------------------------------------------%

nxa=size(Ass_cl,1);
Hiw_cl=zeros(ndof,nf);Hiw_cla=Hiw_cl;
for j=1:nf
    w=wv(j);
    if norm(gamma_c,1)~=0
        Hiw_cla(:,j)=Css_cl*inv(1i*w*eye(nxa)-Ass_cl)*Bss_cl+Dss_cl;
    else
         %open-loop system, u=0
        Hiw_cla(:,j)=C1ss*inv(1i*w*eye(nx)-Ass)*B1ss+D11ss;
    end
    % Change FRF From pressure by volume acceleration 
    % to pressure by volume velocity
    Hiw_cl(:,j)=Hiw_cla(:,j)*1i*w;
end
%%
%-------------------------------------------------------------------------%
%                  Analytical FRF via SEM                                 % 
%-------------------------------------------------------------------------%
feedback='i';%PC properties
%dmodel = input('enter the selected damping model\n(v = viscous, s = structural)')
ideal_filter=0;
filter=0;
plotfilter=0;
ncell=18;
e='m';
%         Boundary conditions opitions:
%            boundary == 0 : free-free (default)
%            boundary == 1 : fixed-fixed
%            boundary == 2 : periodic(infinity system)
boundary=0;
% B=rho*c^2;
% F=B*rho;
F=1;
% gains
% gamma_c=(B*rho).*gamma_c;
H_pv=[1; 1/s; s];
output= function_FRF_Ac_v2(dmodel, rho, L1, L2, A1, A2, c, eta, gamma_c, H_pv, ideal_filter, wv,ncell,e, boundary,F);
%% Plots
figure
scatter(fv,20*log10(abs(u_FEM(1,:))),'r')
hold on
P_SEM=output.ps_v;
scatter(fv,20*log10(abs(Hiw_cl(1,:))),'r*')
scatter(fv,20*log10(abs(P_SEM(1,:))),'rs')
scatter(fv,20*log10(abs(u_FEM(end,:))),'k')
scatter(fv,20*log10(abs(Hiw_cl(end,:))),'k*')
scatter(fv,20*log10(abs(P_SEM(end,:))),'ks')

xlabel('Frequency [Hz]')
ylabel('FRF [dB]')
legend('physical-L','SS-L','SEM-L','physical-R', 'SS-R','SEM-R')
%legend('physical-L','SEM-L','physical-R','SEM-R')
set(gcf, 'Color', 'w');
box on
title('Pressure by volume velocity')
%title('Pressure by volume acceleration')
%export_fig FRF_p_V.pdf.pdf
% %%
% figure
% plot(fv,20*log10(abs(Hiw_cl(end,:))),'k','MarkerSize',8,'LineWidth',1);hold on;
% plot(fv,20*log10(abs(Hiw_cl(1,:))),'r','MarkerSize',8,'LineWidth',1.2)
% xlabel('Frequency [kHz]')
% ylabel('FRF [dB]')
% legend('L','R')
% title('Numerical FRF via SS-formulation')
% set(gcf, 'Color', 'w');
% %% Transient Response: applied force
% % Nc is the number of circles from central frequency fc.
% % T2 define the envelope frequency f2, which has Nc circles within.
% % Nt is the number of entries on time vector t. The heavside function
% % defines a window on half a period of the envelope.
% % time discretization period dt depends on N, in such manner that
% % the total time of the simulation T can be modified if nedeed.
% Nt=4096*2;
% % T=0.01; % total analysis time
% Ap=1;
% fc=1618;
% wc = 2*pi*fc;
% Nc=15; 
% Tc=1/fc;
% T2=Nc*Tc*2; 
% f2=1/T2; w2=2*pi*f2;
% dt=T2/2/4096;
% %dt=T2/Nt;
% t=0:dt:(Nt-1)*dt;
% w=Ap*sin(wc*t).*sin(w2*t).^2.*heaviside(T2/2-t); 
% 
% %% Transient Response: Simulation by integration of physical model
% 
% Ass=[zeros(ndof,ndof) eye(ndof,ndof)
% -inv(Ms)*Ks -inv(Ms)*Cs];
% 
% if e == 'm'
%     m = floor(ndof/2);
%     if mod(ndof,2) == 0
%          Fm = [zeros(ndof-m-1,1);1;zeros(ndof-m,1)];
%     else
%         Fm = [zeros(m,1);1;zeros(m,1)];
%     end
% % excitation on the left end
% elseif e == 'l'
%     Fm = [1; zeros(ndof-1,1)];
% % excitation on the right end    
% elseif e=='r'
%     Fm = [zeros(ndof-1,1);1];
% else 
%     disp('invalid force position')
%     return
% end
% 
% Bss=[zeros(ndof,1)
%      (B*rho).*inv(Ms)*Fm];
% 
% Cqsi=[eye(ndof,ndof) zeros(ndof,ndof)];
% D=0;
% sys_qsi = ss(Ass,Bss,Cqsi,D);
% 
% qsi=lsim(sys_qsi,w,t,zeros(2*ndof,1));
% %% 
% Nplots=100;
% DI=round(length(t)/Nplots);
% tend=length(t);
% figure
% for i=1:DI:tend
%     plot3(x,(t(i).*10^3)*ones(1,length(x)),qsi(i,:),'b');
%      hold on
% end
% xlabel('Length [m]'), ylabel('Time [ms]'), zlabel('Pressure p(x,t)[Pa]')
% hold off
% box on
% set(gcf, 'Color', 'w');
% 
% %% zero state response of the closed loop of SS model for input w
% if norm(gamma_c,1)~=0
%     ni=nx_cl;
% else
%     ni=nx;
% end
% z=lsim(sys_wz,w,t,zeros(ni,1));
% %%
% Nplots=100;
% DI=round(length(t)/Nplots);
% tend=length(t);
% figure
% for i=1:DI:tend
%     plot3(x,(t(i).*10^3)*ones(1,length(x)),z(i,:),'b');
%      hold on
% end
% xlabel('Length [m]'), ylabel('Time [ms]'), zlabel('Pressure p(x,t)[Pa]')
% hold off
% box on
% set(gcf, 'Color', 'w');