% Created on : 05/01/2022
% ERA identification algorithm credited to Arruda 4/1992

%% clear 
clear all 
%close all 
clc 
%% Macros
global dmodel rho L1 L2 Lc A1 A2 c eta gamma_c H_pv kp ki kd ideal_filter kL
%feedback = input('enter the type of feddback\n(string: 
feedback='i';%PC properties
%kp=1e-2;
ki=0*1.5e-2;
%ki=-1.5e-3;
%kd=2e-10;
%dmodel = input('enter the selected damping model\n(v = viscous, s = structural)')
dmodel='v';
ideal_filter=0;
% Laplace variable
s=tf('s');
% angular frequency vector
flim=1.5e3;%in kHz
% parameter used for time derivative aproximation
Nd=100*2*pi*flim;
wv=2*pi*(2:1/4:flim);
%wv=2*pi*(-flim:2:flim);
fv=wv/2/pi;
% pb =input('enter the number of pass bands / bulk bands\n(any integer number between 1 and 10)')
pb=4;
% a=coupling range (in cells)
a=0;
% filter = 0 : no filter
% filter = 1 : butterworth filter
% filter = 2 : RC filter
% filter = 3 : filter built for the experiment
filter=0;
plotfilter=0;
function_metasetup(filter,plotfilter,feedback);
% labelplot(1): wavenumber x real frequency
% labelplot(2): wavenumber x imaginary frequency
% labelplot(3): complex frequency plane
% labelplot(4): 3D plot
labelplot=[0 0 1 0];
%Boundary conditions of the finite structure:
%            boundary == 0 : free-free (default)
%            boundary == 1 : fixed-fixed
%            boundary == 2 : periodic(infinity system)
boundary=0;
% number of cells on the finite structure
ncell=18;
%FEM
ne_cell=7*3;
xs = Lc/4;
xi = xs+Lc/2;
% excitation point for zero-state response
e='m';

%% initial conditions set up for w(k) method
m=function_initialconditions_setup(pb,fv);
% pb=1;
%  df=fv(2)-fv(1);
%  m1=find(abs(fv-150) <= df/2);
%  m2=find(abs(fv+150) <= df/2);
%  m3=find(abs(fv-550) <= df/2);
%  m4=find(abs(fv+550) <= df/2);
%  m=[m1 m2 m3 m4];
m_passive=m;
%% Inverse method via transfer matrix and SEM formulation
kLv=-pb*pi:8*pi/400:pb*pi;
N=length(kLv); 
fun= @caracteristic_equation_T_Acoustic;
options = optimoptions('fsolve','Display','off');
options.StepTolerance =1e-21;
options.OptimalityTolerance = 1e-21;
options.MaxIterations = 100;
w_TM=zeros(pb,N);w_TM_passive=w_TM;

% This loop runs through the vector kLv.
% PS: remember to adjust the initial conditions depending on the PC. You
% will need to adjust wv range also via flim.
        
gamma_c_aux =gamma_c;
for n =1:N
    kL=kLv(n);
    % for non-local feedback, introduce spatial delay of a cells
    gamma_c = gamma_c_aux.*exp(1i*a*kL);
    for i=1:pb
        w_TM(i,n) = fsolve(fun,wv(m(i)),options);
    end

end
gamma_c = gamma_c_aux;
%%Hermitian counterpart(passive system)
gamma_c_aux =gamma_c;
for n =1:N
    kL=kLv(n);
    gamma_c=[0 0 0];
    for i=1:pb
        w_TM_passive(i,n) = fsolve(fun,wv(m_passive(i)),options);
    end

end
gamma_c = gamma_c_aux;


% Detach frequency and damping factor for each bulk band and converting to kHz:
for i=1:pb
    w_real_TM{i} = real(w_TM(i,:))/2/pi/1000;
    w_imag_TM{i} = imag(w_TM(i,:))/2/pi/1000;
    w_real_TM_passive{i} = real(w_TM_passive(i,:))/2/pi/1000;
    w_imag_TM_passive{i} = imag(w_TM_passive(i,:))/2/pi/1000;
end

%% FEM
output = function_buildFEM(rho, L1, L2, Lc, A1, A2, c, eta,ne_cell,ncell,flim,xs,xi...
    ,'dmodel',dmodel,'boundary',boundary);
M=output.M;C=output.C;K=output.K;
% Passive structure matrices
Ms=output.Ms;Cs=output.Cs;Ks=output.Ks;
x=output.x;ndof=output.ndof;
n1=output.n1;n2=output.n2;np=output.np;
%% State Space of interconnection (S,H) 
%% S
output=function_buildSS(Ms,Cs,Ks,ndof,ncell,n1,n2,e,a);
Ass=output.Ass;nx=size(Ass,1);
% PS: (Ms,Cs,Ks) were previously multiplied by B*rho,
% so we need to multiply matrixes F (inside B1ss) and T(inside B2ss) also 
B=rho*c^2;
B1ss=(B*rho).*output.B1ss;nw=size(B1ss,2);
B2ss=(B*rho).*output.B2ss;nu=size(B2ss,2);
C1ss=output.C1ss;nz=size(C1ss,1);
C2ss=output.C2ss;ny=size(C2ss,1);
D11ss=output.D11ss;
D12ss=output.D12ss;
D21ss=output.D21ss;
%% H 
% Feedback law in pressure by volume velocity
H_v=gamma_c*H_pv;
% Correction to volume acceleration
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
wn_FEM = -1i*Lambda;
[wn_FEM,ind] = sort(wn_FEM,'ComparisonMethod','real');
Vn = Vn(:,ind); 

%% Dispersion Diagram Aproximation via FEM
% number of internal nodes ni
ni=np-2;
% PS: In this discretization, source is -B/A times volume velocity.
% Also, the matrixes (Ms,Cs,Ks) was multiplied by A*rho
gamma_c_ph=(B*rho).*gamma_c;
[w_real_FEM,w_imag_FEM]= function_dipersion_FEM(ni,pb,M,C,K,gamma_c_ph,a,dmodel,kLv);
%%
%-------------------------------------------------------------------------%
%                            FRF via FEM      
% External force in unitary volume acceleration 
% PS: In this discretization, source is -B/A times volume velocity.
% Also, the matrixes (Ms,Cs,Ks) were previously multiplied by A*rho
F=B*rho;
P_FEM= function_FRF_FEM(F,wv,M,C,K,e,gamma_c_ph,ndof,ne)
%% Identifying eigenmodes via SEM (FRF)
e='m';
F=1;
output=function_FRF_Ac(wv,ncell,e,boundary,F);
%% FRF pressure by volume velocity
P_SEM_passive=output.ps_passive.';
P_SEM=output.ps_v.';
figure
plot(fv/1000,20*log10(abs(P_SEM(:,1))),'r','MarkerSize',8,'LineWidth',1)
hold on
%plot(fv/1000,20*log10(abs(P_SEM_passive(:,1))),'r--','MarkerSize',8,'LineWidth',1)
plot(fv/1000,20*log10(abs(P_SEM(:,end))),'k','MarkerSize',8,'LineWidth',1)
%plot(fv/1000,20*log10(abs(P_SEM_passive(:,end))),'k--','MarkerSize',8,'LineWidth',1)
ylabel('$ p / G $ [dB]', 'interpreter', 'latex', 'fontsize', 15)
xlabel('$\Re(f)$ [kHz]', 'interpreter', 'latex', 'fontsize', 15)
set(gca,'TickLabelInterpreter','Latex','fontsize',15);
%legend('left side','Passive','right side')
%xlim([0 1.3])
box on
set(gcf, 'Color', 'w');
grid on
%%
Nv=length(fv);
H=P_SEM_passive(:,1);
H=[real(H(1));H(1:Nv);0;flipud(conj(H(1:Nv)))];
N=length(H);
df=fv(2)-fv(1);
h=real(N*ifft(H));
ts=1/(N*df);
%T=1/(2*1500);
t=(0:(N-1))*ts;
figure
plot(t,h)
title('Impulse response ', 'interpreter', 'latex', 'fontsize', 15)
box on
set(gcf, 'Color', 'w');
grid on
[Ad,Bd,Cd,Dd]=erasiso(ts*h');
% Ex: 
% No. de colunas da matriz de Hankel:200 (ordem do sistema estimado inicialmente)
% No. e de vezes que h e sobreposta na matriz de Hankel: 1000 
% No. e de raizes (picos observados na FRF): 80 (ordem depois do truncamento
% seguindo a decomposição SVD)
sd=ss(Ad,Bd,Cd,Dd,ts);
sc=d2c(sd);
[A,B,C,D]=ssdata(sc);
s=eig(A);
wn_SEM_ident=-1i*s;
%wn_SEM_ident=sort(abs(s))/2/pi;

%%
%-------------------------------------------------------------------------%
%    Unfolded Dispersion Diagram Plots by SEM  with Eigenmodes  by FEM   
%                               and SEM                                   
%-------------------------------------------------------------------------%
% labelplot(1): wavenumber x real frequency
% labelplot(2): wavenumber x imaginary frequency
% labelplot(3): complex frequency plane
% labelplot(4): 3D plot
function_unfolded_dispersion_diagrams(pb,kLv,w_real_TM,w_imag_TM,w_real_TM_passive,w_imag_TM_passive...
   ,'setplot',labelplot,'periodic',1);
zlim([-1 1])
hold on
scatter(imag(wn_FEM/2/pi/1000),real(wn_FEM/2/pi/1000),'k')
scatter(imag(wn_SEM_ident/2/pi/1000),real(wn_SEM_ident/2/pi/1000),'*k')
 ylim([0 0.5])
%  xlim([-0.3 0.3]);
