% Created on : 04/02/2022 by Danilo Braghini
%% clear 
clear all 
%close all 
clc 
%% Macros
global ne rho1 rho2 A1 A2 kL E1 E2 Iz1 Iz2 R s01 s02 H_fu ideal_filter  eta
%feedback = input('enter the type of feddback\n(string: 
feedback='p';%PC properties
kp=1e4;
%ki=1.5e-2;
%ki=-1.5e-3;
%kd=2e-10;
ideal_filter=0;
% angular frequency vector
flim=1e3;%in kHz
Nd=100*2*pi*flim;
plotlog=0;
if plotlog==1
    fv=logspace(-1, 4, 500);
    wv=fv*2*pi;
else
    wv=2*pi*(1:2:flim);
    fv=wv/2/pi;
end
% pb =input('enter the number of pass bands / bulk bands\n(any integer number between 1 and 10)')
pb=1;
% a=coupling range (in cells)
a=0;
% labelplot(1): wavenumber x real frequency
% labelplot(2): wavenumber x imaginary frequency
% labelplot(3): complex frequency plane
% labelplot(4): 3D plot
labelplot=[1 1 1 1];
periodic=0;
passive=0;
colors=0;
% Ring(closed curved beam) type metamaterial set up
%R                           radius [m]
%bh                          width [m]
%hh                          height [m]
%rhoh                        mass density [kg/m3]
%E                           elastic modulus [Pa]
%s0c                          circumferential length of the unit cell
% Spectral Element Mesh %%%%%%%%%%%
%Nuc=4;
Nuc = 10;                       % number of substructures or unit cells
dtheta = 2*pi/Nuc;              % angular length of the unit cell
%eta=0.001;
eta=1e-8;
% Section 1
h1=0.001;
%h1=0.002;
b1=0.010;
%b1=0.15;
R=0.050;
%R=0.25;
s0c=R*dtheta; 
%s01=s0c;
% number of SEM elements per cell;
ne=3;
% curve length of each segment
if ne==1
    s01=s0c;
else
    s01=(s0c/2)/2;
end
E1=5*10^9;
%E1=220e9;
rho1=1500;
%rho1=7200;
A1=b1*h1;
Iz1=b1*h1^3/12;
% Section 2
h2=h1;
b2=b1;
s02=2*s01;
E2=E1;
rho2=rho1;
A2=b2*h2;
Iz2=b2*h2^3/12;

%% setting feedback law
% gains
s=tf('s');
if feedback=='p'
    gamma_c=[kp 0 0];
    %gamma_c=[0 0 0];
elseif feedback=='i'
    gamma_c=[0 ki 0];
elseif feedback=='d'
     gamma_c=[0 0 kd];
elseif feedback=='pi'
    gamma_c=[kp ki 0];
end
PID=[1; 1/s; Nd/(1+Nd/s)];
H_fu=zeros(3,3).*s;
H_fu(1,1)=gamma_c*PID;
%% Dispersion Relation of passive system (gamma_c = 0)
%use this section to define the initial values of fsolve optimization
[kL_sem_PB,kL_sem_SB] = function_SEM_EBRing(wv);
%%
if plotlog==0
    figure
    plot(kL_sem_PB/pi,fv/1000,'.');
    hold on
    ylabel('$f$ [kHz]', 'interpreter', 'latex', 'fontsize', 15)
    xlabel('$\Re(k L_c/ \pi)$', 'interpreter', 'latex', 'fontsize', 15)
    set(gca,'TickLabelInterpreter','Latex','fontsize',15);
    box on
    grid on
    set(gcf, 'Color', 'w');
    figure
    plot(kL_sem_SB/pi,fv/1000,'.');
    hold on
    ylabel('$f$ [kHz]', 'interpreter', 'latex', 'fontsize', 15)
    xlabel('$\Im(k L_c/ \pi)$', 'interpreter', 'latex', 'fontsize', 15)
    set(gca,'TickLabelInterpreter','Latex','fontsize',15);
    box on
    grid on
    set(gcf, 'Color', 'w');
elseif plotlog==1
    figure
    semilogy(kL_sem_PB/pi,fv/1000,'.');
    hold on
    ylabel('$f$ [kHz]', 'interpreter', 'latex', 'fontsize', 15)
    xlabel('$\Re(k L_c/ \pi)$', 'interpreter', 'latex', 'fontsize', 15)
    set(gca,'TickLabelInterpreter','Latex','fontsize',15);
    box on
    grid on
    set(gcf, 'Color', 'w');
    figure
    semilogy(kL_sem_SB/pi,fv/1000,'.');
    hold on
    ylabel('$f$ [kHz]', 'interpreter', 'latex', 'fontsize', 15)
    xlabel('$\Im(k L_c/ \pi)$', 'interpreter', 'latex', 'fontsize', 15)
    set(gca,'TickLabelInterpreter','Latex','fontsize',15);
    box on
    grid on
    set(gcf, 'Color', 'w');
end
%%
% figure
% plot(fv/1000,kL_sem_PB/pi,'k+',fv/1000,kL_sem_SB/pi,'r+');
% xlabel('$\Re(f)$ [kHz]', 'interpreter', 'latex', 'fontsize', 15)
% ylabel('$k L_c/ \pi$', 'interpreter', 'latex', 'fontsize', 15)
% set(gca,'TickLabelInterpreter','Latex','fontsize',15);
% box on
% grid on
% set(gcf, 'Color', 'w');
 %%
% figure
% plot(fv,kL_sem_PB,'k+');
% xlabel('$\Re(f)$ [kHz]', 'interpreter', 'latex', 'fontsize', 15)
% ylabel('$k L_c/ \pi$', 'interpreter', 'latex', 'fontsize', 15)
% set(gca,'TickLabelInterpreter','Latex','fontsize',15);
% box on
% grid on
% % set(gcf, 'Color', 'w');
%% initial conditions set up for w(k) method
%m=function_initialconditions_setup(pb,fv);
%   pb=1;
df=fv(2)-fv(1);
%m1=find(abs(fv-4000) <= df/2);
%m2=find(abs(fv-8000) <= df/2);
%m=[m1 m2];
m= find(abs(fv-780) <= df/2);
%m_passive=m;
%% Inverse method via transfer matrix and SEM formulation
kLv=-pb*pi:pb*2*pi/100:pb*pi;
N=length(kLv); 
fun= @caracteristic_equation_T_EBRing;
options = optimoptions('fsolve','Display','off');
options.StepTolerance =1e-21;
options.OptimalityTolerance = 1e-21;
options.MaxIterations = 100;
w_TM=zeros(pb,N);w_TM_passive=w_TM;
% This loop runs through the vector kLv.
% PS: remember to adjust the initial conditions depending on the PC. You
% will need to adjust wv range also via flim.

% for i=1:pb
% gamma_c_aux =gamma_c;
%     wik=wv(m);
%     for k =1:N
%         kL=kLv(k);
%         % for non-local feedback, introduce spatial delay of a cells
%         gamma_c = gamma_c_aux.*exp(1i*a*kL);
%         w_TM(i,k) = fsolve(fun,wik,options);
% %         if kL==0
% %             wik=wv(m);
% %         else
%             wik=real(w_TM(i,k));
%         %end
%     end
% gamma_c = gamma_c_aux;
% end

gamma_c_aux =gamma_c;
for n =1:N
    kL=kLv(n);
    gamma_c = gamma_c_aux.*exp(1i*a*kL);
    for i=1:pb
        w_TM(i,n) = fsolve(fun,wv(m(i)),options);
    end

end
gamma_c = gamma_c_aux;

% Detach frequency and damping factor for each bulk band and converting to kHz:
for i=1:pb
    w_real_TM{i} = real(w_TM(i,:))/2/pi/1000;
    w_imag_TM{i} = imag(w_TM(i,:))/2/pi/1000;
%     w_real_TM_passive{i} = real(w_TM_passive(i,:))/2/pi/1000;
%     w_imag_TM_passive{i} = imag(w_TM_passive(i,:))/2/pi/1000;
end
%% Plots
%-------------------------------------------------------------------------%
%          Unfolded Dispersion Diagram Plots by SEM        %
%-------------------------------------------------------------------------%

function_unfolded_dispersion_diagrams(pb,kLv,w_real_TM,w_imag_TM,w_real_TM,w_imag_TM...
    ,'setplot',labelplot,'periodic',periodic,'plotpassive',passive,'plotcolors',colors);
%zlim([-1 1])
