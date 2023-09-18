% Created on : 22/02/2022 by Danilo Braghini
%% clear 
clear all 
%close all 
clc 
%% Macros
global ne rho1 rho2 L1 L2 A1 A2 kL E1 E2 Iz1 Iz2  gamma_c H_fu ideal_filter  eta
%feedback = input('enter the type of feddback\n(string: 
feedback='p';%PC properties
kp=1e6;
%kp=-1e-2;
ki=1.5e-2;
%ki=-1.5e-3;
kd=2e-10;
%dmodel = input('enter the selected damping model\n(v = viscous, s = structural)')
dmodel='n';
ideal_filter=0;
% angular frequency vector
flim=10e3;%in kHz
Nd=100*2*pi*flim;
wv=2*pi*(2:2:flim);
%wv=2*pi*(-flim:2:flim);
fv=wv/2/pi;
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
passive=1;
colors=0;
% beam type metamaterial set up
% if dmodel == 's'
%     eta=0.001;
% elseif dmodel == 'v'
%    eta= 100;
% else
%     eta=0;
% end
%%
eta=1e-8;
ne=3;
% Section 1
h1=0.001;
b1=0.01;
E1=5e9;
rho1=1500;
eta1=0.01;
E1=E1*(1+eta1*1i);
A1=b1*h1;
Iz1=b1*h1^3/12;

% Section 2
h2=h1;
%h2=0.01;
b2=b1;
E2=E1;
rho2=rho1;
%eta2=0.01;
A2=b2*h2;
Iz2=b2*h2^3/12;
Lc = (2*pi)/10*0.05;
xs = Lc/4;
xi = xs+Lc/2;
L1= xs;
L2 = xi- xs;

%% setting feedback law
% gains
s=tf('s');
if feedback=='p'
    gamma_c=[kp 0 0];
elseif feedback=='i'
    gamma_c=[0 ki 0];
elseif feedback=='d'
     gamma_c=[0 0 kd];
elseif feedback=='pi'
    gamma_c=[kp ki 0];
end
H_fu=[1; 1/s; Nd/(1+Nd/s)];
%% Dispersion Relation
%
%-------------------------------------------------------------------------%
%                               SEM   Model                               %  
%-------------------------------------------------------------------------%
%use this section to define the initial values of fsolve optimization
[kL_sem_PB,kL_sem_SB] = function_SEM_EBBeam(wv);
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



%% initial conditions set up for w(k) method
%m=function_initialconditions_setup(pb,fv);
%   pb=1;
df=fv(2)-fv(1);
%   m=find(abs(fv-1200) <= df/2);
%   m2=find(abs(fv-550) <= df/2);
% m=[m1 m2];
m= find(abs(fv-3000) <= df/2);
m_passive=m;
%% Inverse method via transfer matrix and SEM formulation
kLv=-pb*pi:8*pi/400:pb*pi;
N=length(kLv); 
fun= @caracteristic_equation_T_EBBeam;
options = optimoptions('fsolve','Display','off');
options.StepTolerance =1e-21;
options.OptimalityTolerance = 1e-21;
options.MaxIterations = 100;
w_TM=zeros(pb,N);w_TM_passive=w_TM;
clear w_TM
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
%Hermitian counterpart(passive system)
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
%% Plots
%-------------------------------------------------------------------------%
%          Unfolded Dispersion Diagram Plots by SEM        %
%-------------------------------------------------------------------------%

function_unfolded_dispersion_diagrams(pb,kLv,w_real_TM,w_imag_TM,w_real_TM_passive,w_imag_TM_passive...
    ,'setplot',labelplot,'periodic',periodic,'plotpassive',passive,'plotcolors',colors);
%zlim([-1 1])
