clc
%close all
clear
%%
feedback='p';
%Number of eigenvalues or number of modes
Ntrunc=10;
%Number of Fourier coefficients
% Nfou=2*Neig+1;
%Plane waves evaluated
Npw=-Ntrunc:Ntrunc;
%Wavenumber discretization
deltak=100;
%Inclusion area
D=0.04;
%Total length
L=50e-2;
%Speed of sound
c=343;
%Varying area through x-direction
x=linspace(-L/2,L/2,2*Ntrunc+1);
S=pi*D^2/4*ones(size(x));

%Range of m and mbar parameters: g=2*pi*m/a and gbar=2*pi*mbar/a
[m,mbar]=meshgrid(Npw,Npw);
%g and gbar difference
g_gbar=(2*pi/L)*(m-mbar);
Scoef=fft(S)/length(x);
A=toeplitz(Scoef);

%Evaluation of K(omega)
%Matrix G for matricial evaluation
G=2*pi/L*m;
%Wavenumber vector
pb=4;
k=-pb*pi/L:2*pi/L/(deltak):pb*pi/L;
wreal=zeros(Ntrunc,length(k));
wimag=zeros(Ntrunc,length(k));
%The matrix A corresponds to the "mass matrix" in the eigenproblem
M=A/c^2;
x2=3*L/4;
x1=L/4;
% n=0;
rho=1.225;
BB=rho*c^2;
kp=1e-6;
ki=1.5e-2;
kd=2e-10;
if feedback=='d'
    kfb=kd;
elseif feedback=='i'
    kfb=ki;
else
    kfb=kp;
end
gamma=rho*kfb/(pi*D^2/4);
%gamma=1.5e-2*(pi*D^2/4*L)/rho;
a=0;
for i=1:length(k)
   
    %"Stiffness matrix" in the eigenproblem.
    K=(k(i)*ones(length(Npw))+G).*(k(i)*ones(length(Npw))+G.').*A;

    B=exp(1j*x2*G);
    UM=exp(-1j*x1*(k(i)+G(1,:)));
    C=gamma/L*exp(1j*k(i)*(a*L+x2))*sum(B.*A,2)*UM;

    %The eigenvalues for each wavenumber
    if feedback == 'd'
           omega2=eig(K,M-C);
    elseif feedback == 'i'
        omega2=eig(K-C,M);
    elseif feedback == 'p'
         omega2=polyeig(K,C,M);
    end

    %Ordering the eigenvalues
    omega2=sort(omega2,'ComparisonMethod','real');
    %Frequencies for the wavenumber
    wreal(:,i)=real(sqrt(omega2(1:Ntrunc)));
    wimag(:,i)=imag(sqrt(omega2(1:Ntrunc)));
    if k(i)==0
        wimag(:,i)=zeros(Ntrunc,1);
    end
end
%% Build structeres to plot function 
for i=1:Ntrunc
    w_real_PWE{i} = wreal(i,:)/2/pi/1000;
    w_imag_PWE{i} = wimag(i,:)/2/pi/1000;
end
kLv=k.*L;
%% Plots
%-------------------------------------------------------------------------%
%          Unfolded Dispersion Diagram Plots by SEM        %
%-------------------------------------------------------------------------%
% labelplot(1): wavenumber x real frequency
% labelplot(2): wavenumber x imaginary frequency
% labelplot(3): complex frequency plane
labelplot=[1 1 1 0];
function_unfolded_dispersion_diagrams(pb,kLv, w_real_PWE,w_imag_PWE,w_real_PWE,w_imag_PWE,'setplot',labelplot);

