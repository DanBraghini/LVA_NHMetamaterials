%%
clc
close all
clear
%% Number of eigenvalues or number of modes
Ntrunc=4;
%Number of Fourier coefficients
% Nfou=2*Neig+1;
%Plane waves evaluated
Npw=-Ntrunc:Ntrunc;
%Wavenumber discretization
deltak=50;
%Inclusion area
D=0.04;
%Total length
L=50e-3;
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

%% Evaluation of K(omega)
%Matrix G for matricial evaluation
G=2*pi/L*m;
%Wavenumber vector
k=-Ntrunc*pi/L:2*pi/L/(deltak):Ntrunc*pi/L;
freqreal=zeros(Ntrunc,length(k));
freqimag=zeros(Ntrunc,length(k));
%The matrix A corresponds to the "mass matrix" in the eigenproblem
M=A/c^2;
x1=3*L/4;
x2=L/4;
% n=0;
gamma=-.5e1;
a=0;rho=1.225;
p=2*pi*2100;
z=2*pi*1000;
% w0=linspace(0,5000*2*pi,2*Ntrunc+1)';
[K0,C0]=pwematrices(0,G,A,Npw,L,x1,x2,a,gamma);
[KL,CL]=pwematrices(pi/L,G,A,Npw,L,x1,x2,a,gamma);
w0=(sqrt(sort(eig(K0,M)))+sqrt(sort(eig(KL,M))))/2;
%%
options = optimoptions('fsolve','Display','off');
options.StepTolerance =1e-20;
options.OptimalityTolerance = 1e-20;
options.MaxIterations = 100;
xn=1e-5;
%%
for i=1:length(k)
   
    %"Stiffness matrix" in the eigenproblem.
%     K=(k(i)*ones(length(Npw))+G).*(k(i)*ones(length(Npw))+G.').*A;
% 
%     B=exp(-1j*x1*G);
%     UM=exp(1j*x2*(k(i)+G(1,:)));
%     C=gamma/L*exp(-1j*k(i)*(a*L+x1))*sum(B.*A,2)*UM;
    [K,C]=pwematrices(k(i),G,A,Npw,L,x1,x2,a,gamma);
    %The eigenvalues for each wavenumber
%     omega2=eig(K,M-C);%D
%     omega2=eig(K-C,M);%I
%     omega2=polyeig(K,C,M);%P
    fun=@(w)(det(K-w^2*M-.5*(1j*w+z)/(1j*w+p)*C));
%     fun=@(w)(det(K-w^2*M-(xn)/(xn/1j*w+1)*C));
%     fun=@(w)(det(K-C-w^2*M));
    omega2=zeros(2*Ntrunc+1,1);
%     w0=sqrt(sort(eig(K,M)));
    for j=1:2*Ntrunc+1
        omega2(j)=(fsolve(fun,w0(j),options))^2;
    end
    %Ordering the eigenvalues
    omega2=sort(omega2,'ComparisonMethod','real');
    %Frequencies for the wavenumber
    freqreal(:,i)=real(sqrt(omega2(1:Ntrunc)))/2/pi;
    freqimag(:,i)=imag(sqrt(omega2(1:Ntrunc)))/2/pi;
    
end

%% Build structeres to plot function 
for i=1:Ntrunc
    f_real_PWE{i} = freqreal(i,:);
    f_imag_PWE{i} = freqimag(i,:);
end
kLv=k.*L;
%% Plots
%-------------------------------------------------------------------------%
%          Unfolded Dispersion Diagram Plots by SEM        %
%-------------------------------------------------------------------------%
% labelplot(1): wavenumber x real frequency
% labelplot(2): wavenumber x imaginary frequency
% labelplot(3): complex frequency plane
labelplot=[1 1 1];
function_unfolded_dispersion_diagrams(Ntrunc,kLv, f_real_PWE,f_imag_PWE,f_real_PWE,f_imag_PWE,labelplot);

%%
function [K,C]=pwematrices(k,G,A,Npw,L,x1,x2,a,gamma)

K=(k*ones(length(Npw))+G).*(k*ones(length(Npw))+G.').*A;

B=exp(-1j*x1*G);
UM=exp(1j*x2*(k+G(1,:)));
C=gamma/L*exp(-1j*k*(a*L+x1))*sum(B.*A,2)*UM;

end
