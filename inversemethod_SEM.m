% Created on : 15/02/2021 by Danilo Braghini
%% clear 
clear all 
%close all 
clc 
%% Macros
global rho L1 L2 A1 A2 kL c eta gamma_c dmodel ideal_filter 
% angular frequency vector
flim=3.5e3;%in kHz
wv=2*pi*(2:2:flim);
fv=wv/2/pi;
dmodel='s';
% pb =input('enter the number of pass bands / bulk bands\n(any integer number between 1 and 10)')
pb=4;
ideal_filter=0;
gamma_c=0;
%% acoustic metamaterial set up
eta=0.001;
rho=1.225;
c=343;
Lc = 50e-2;
r1 = 2e-2;
%r1 = 3e-2;
A1 = pi*r1^2;
r2 = r1;
A2 = 2*A1;
xs = Lc/4;
xi = xs+Lc/2;
L1= xs;
L2 = xi- xs;
%% Dispersion Relation
%
%-------------------------------------------------------------------------%
%                               SEM   Model                               %  
%-------------------------------------------------------------------------%
%use this section to define the initial values of fsolve optimization
[kL_sem_PB,kL_sem_SB] = function_SEM_Acoustic(wv);
figure
plot(real(fv)/1000,abs(kL_sem_PB)/pi,'k--',real(fv)/1000,-abs(kL_sem_SB)/pi,'r--');grid on
hold on
xlabel('$\Re(f)$ [kHz]', 'interpreter', 'latex', 'fontsize', 15)
ylabel('$k L_c/ \pi$', 'interpreter', 'latex', 'fontsize', 15)


%% initial conditions set up for w(k) method
% define m as a vector with the initial conditions, frequencies within each
% one of the bulk bands
m=[100 500 850 1100];
%m=function_initialconditions_setup(pb,fv);
%% Inverse method via transfer matrix and SEM formulation
kLv=-pb*pi:8*pi/400:pb*pi;
N=length(kLv); 
fun= @caracteristic_equation_T_Acoustic;
options = optimoptions('fsolve','Display','off');
options.StepTolerance =1e-21;
options.OptimalityTolerance = 1e-21;
options.MaxIterations = 100;
w_TM=zeros(pb,N);
% This loop runs through the vector kLv.
% PS: remember to adjust the initial conditions depending on the PC. You
% will need to adjust wv range also.
        
%Hermitian counterpart(passive system)
gamma_c_aux =gamma_c;
for n =1:N
    kL=kLv(n);
    gamma_c=0;
     for i=1:pb
        w_TM(i,n) = fsolve(fun,wv(m(i)),options);
    end

end
% Detach frequency and damping factor for each bulk band and converting to kHz:
for i=1:pb
    w_real_TM{i} = real(w_TM(i,:))/2/pi/1000;
    w_imag_TM{i} = imag(w_TM(i,:))/2/pi/1000;
end

%% Plots
%-------------------------------------------------------------------------%
%          Unfolded Dispersion Diagram Plots by SEM        %
%-------------------------------------------------------------------------%
% labelplot(1): wavenumber x real frequency
% labelplot(2): wavenumber x imaginary frequency
% labelplot(3): complex frequency plane
% labelplot(4): 3D plot
labelplot=[1 0 0 0];
function_unfolded_dispersion_diagrams(pb,kLv,w_real_TM,w_imag_TM,w_real_TM,w_imag_TM...
    ,'setplot',labelplot,'periodic',1);
%zlim([-1 1])
