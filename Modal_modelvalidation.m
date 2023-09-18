% Elastic Structurer example. Homogenious rod under free-free boundary
% conditions
% modal (lumped) model validation against FEM (lumped) model
close all
clear all
clc
%% set up
%polymer
E=4.35*10^9; 
rho=1180;
% no damping
eta=0;
Lc = 50e-2;
r = 2e-2;
A = pi*r^2;
xs = Lc/4;
xi = xs+Lc/2;
L1= xs;
L2 = xi- xs;
%% FEM
% ncell:number of cells that make the structure (multiple of 3 to be divided between the 3 segments)
% number of elements per cell
ne_cell=25*3;
%number of cells on the structure
ncell=1;
% number of elements
ne_1 = ne_cell/3; ne_2 = ne_1;
%ne==ne_cell
ne = ne_2 + 2*ne_1;
% np: number of nodes on the mesh per unity cell
ndof = ne + 1;
% lengths of each element
Le_2=L2/ne_2; Le_1 = L1/ne_1;
% indexing grid
n1 = ne_1 +1; n2 = ne_1 + ne_2 +1;
Ms=zeros(ndof,ndof);Ks=Ms;Cs=Ms;
Me1=(rho*A*Le_1/6)*[ 2 1
                    1 2];
Me2=(rho*A*Le_2/6)*[ 2 1
                    1 2];
Ke1=(E*A/Le_1)*[1 -1
                -1 1];
Ke2=(E*A/Le_2)*[1 -1
                -1 1];
Ce1=(eta*A*Le_1/6)*[ 2 1
                    1 2] ;
Ce2=(eta*A*Le_2/6)*[ 2 1
                    1 2] ;
for i=1:ne_1
    Ms(i:i+1,i:i+1) = Ms(i:i+1,i:i+1)+ Me1;
    Ks(i:i+1,i:i+1) = Ks(i:i+1,i:i+1)+ Ke1;
    Cs(i:i+1,i:i+1) = Cs(i:i+1,i:i+1)+ Ce1;
end
for i=n1:ne_1+ne_2
    Ms(i:i+1,i:i+1) = Ms(i:i+1,i:i+1)+ Me2;
    Ks(i:i+1,i:i+1) = Ks(i:i+1,i:i+1)+ Ke2;
    Cs(i:i+1,i:i+1) = Cs(i:i+1,i:i+1)+ Ce2;
end
for i=n2:ne
    Ms(i:i+1,i:i+1) = Ms(i:i+1,i:i+1)+ Me1;
    Ks(i:i+1,i:i+1) = Ks(i:i+1,i:i+1)+ Ke1;
    Cs(i:i+1,i:i+1) = Cs(i:i+1,i:i+1)+ Ce1;
end
Cs=0.*Ks;
% Spatial vector x 
xA1=0:Le_1:xs;
xB1=xs+Le_2:Le_2:xi;
xA2=xi+Le_1:Le_1:Lc;

xcell = [xA1 xB1 xA2];
x = xcell;
for i = 1:ncell-1
    xaux = xcell + i*Lc + Le_1;
    x = [x xaux];
    x(end) = [];
end

%% modal superposition
% number of modes considered
N=100; 
% number of sensors/actuators
M=1;
% modal mass and stiffness vectors
Mj=rho*Lc/2/E; Kj=zeros(N,1);
% B and W matrixes
B=zeros(N,M);W=zeros(N,N);
% C matrix
Cm=(eta/rho).*eye(N);
% space vector s with sensors positions
if M==1
    s=xs;
else
    for j=1:M
        s=linspace(0,Lc,M);
    end
end
for i =1:N
    %Kj(i)=(i*pi)^2/2/Lc;
    W(i,i)=E*(i*pi)^2/(rho*Lc^2);
    for j=1:M
        B(i,j) = cos(i*pi*s(j));
    end
    
end

%% External Force
% Defining tone-burst excitation
% Nc is the number of circles from central frequency fc.
% T2 define the envelope frequency f2, which has Nc circles within.
% Nt is the number of entries on time vector t. The heavside function
% defines a window on half a period of the envelope.
% time discretization period dt depends on N, in such manner that
% the total time of the simulation T can be modified if nedeed.
Nt=2*4096;
% T=0.01; % total analysis time
Ap=1;
fc=7000;
wc = 2*pi*fc;
Nc=15; 
Tc=1/fc;
T2=Nc*Tc*2; 
f2=1/T2; w2=2*pi*f2;
dt=Tc/60;
% dt = T/N;
t=0:dt:(Nt-1)*dt;
% t=0:dt:(Nt-1)*dt;
% criar funçao heaviside mais confiável
w=Ap*sin(wc*t).*sin(w2*t).^2.*heaviside(T2/2-t); 
df=1/t(end); 
fm=1/dt;
 fvec=0:df:fm/2;
% fN=N*df/2;
% f=0:df:fN-1;
Fw=1/Nt*fft(w);

figure
subplot(1,2,1)
plot(t,w)
title('external force dynamics')
subplot(1,2,2)
plot(fvec/1000,abs(Fw(1:Nt/2)))
xlim([0 10])
title('external force spectrum')
%% Integration of the system in time using state space formulation
%% FEM
A_FEM=[zeros(ndof,ndof) eye(ndof)
-inv(Ms)*Ks -inv(Ms)*Cs];
% excitation at s=0
F = [1; zeros(ndof-1,1)];
B_FEM = [zeros(ndof,1)
         inv(Ms)*F];

C_FEM=[eye(ndof,ndof) zeros(ndof,ndof)];
D_FEM=zeros(ndof,1);
sys_FEM = ss(A_FEM,B_FEM,C_FEM,D_FEM);
%zero input response
x_FEM=lsim(sys_FEM,w',t,zeros(2*ndof,1));
%% modal
A_modal=[zeros(N,N) eye(N)
          -W        zeros(N,N)];
% excitation at s=0
F = [2/rho/A/Lc; zeros(N-1,1)];     
B_modal= [zeros(N,M);
          F];
C_modal= [B' zeros(M,N)];
D_modal= zeros(M,M);
sys_modal = ss(A_modal,B_modal,C_modal,D_modal);
%zero input response
x_modal=lsim(sys_modal,w',t,zeros(2*N,1));
%% sensor displacement
figure
plot(t,x_FEM(:,n1).*1e9)
hold on
plot(t,x_modal.*1e9)
grid on
ylabel('displacement [nm]')
xlabel('time [s]')
legend('FEM','modal')

