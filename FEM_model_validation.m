%% clear 
clear all 
%close all 
clc 
%% 
% Created on : 15/02/2021 by Danilo Braghini
%feedback = input('enter the type of feddback\n(d = derivative, i = integrative and p = proportional)\n','s');
feedback='pid';
%dmodel = input('enter the selected damping model\n(v = viscous, s = structural)')
dmodel='v';
%%
% angular frequency vector
flim=1.5e3;
wv=2*pi*(2:2:flim);
fv = wv./(2*pi);
nf=length(wv);
%% acoustic metamaterial set up
%PC properties
if dmodel == 's'
    eta=0.01;
else
   eta= 10;
end
%eta=0;
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
% Feedback control law and sensor locality paramater 
kp=1e-7;
ki=1.5e-2;
kd=2e-10;
gamma_c=[0*kp ki 0*kd];
a=0;
%%
% inputs:
% ne: number of finite elements per unity cell
% np: number of nodes on the mesh per unity cell
% n1: index for the last node of segment 1 on the unity cell
% n2: index for the last node of segment 2 on the unity cell
% ncell: number of cells that make the structure

%FEM
% number of elements for sensor (S) and actuator (A) PZT materials
ne_1 = 5; ne_2 = ne_1;
%ne_1 = 5; ne_2 = ne_1;
ne = ne_2 + 2*ne_1;
np = ne + 1;
% lengths of each element
Le_2=L2/ne_2; Le_1 = L1/ne_1;
% indexing grid
n1 = ne_1 +1; n2 = ne_1 + ne_2 +1;
M=zeros(np,np);K=M;

% mass matrixes
Me1 = rho*A1*Le_1/6*[ 2 1
                      1 2] ;
Me2 = rho*A2*Le_2/6*[ 2 1
                     1 2];  
% viscous damping matrixes   
if dmodel=='v'
    C=M;
    Ce1=eta*A1*Le_1/6*[2 1
                      1 2];                     
    Ce2=eta*A2*Le_2/6*[2 1
                       1 2];  
end
% stiffness matrixes
Ke1 =A1*B/Le_1*[1 -1
                -1 1 ];
Ke2 =A2*B/Le_2*[1 -1
                -1 1 ];   

% assembling the unity cell
% this is a loop through elements, not nodes. (each i is one element in the one segment, each element adds a 2x2 matrix) 
for i=1:ne_1
    M(i:i+1,i:i+1) = M(i:i+1,i:i+1)+ Me1;
    if dmodel=='v'
       C(i:i+1,i:i+1) = C(i:i+1,i:i+1)+ Ce1;
    end
    K(i:i+1,i:i+1) = K(i:i+1,i:i+1)+ Ke1;
end
for i=ne_1+1:ne_1+ne_2
    M(i:i+1,i:i+1) = M(i:i+1,i:i+1)+ Me2;
    if dmodel=='v'
        C(i:i+1,i:i+1) = C(i:i+1,i:i+1)+ Ce2;
    end
    K(i:i+1,i:i+1) = K(i:i+1,i:i+1)+ Ke2;
end
for i=n2:ne
    M(i:i+1,i:i+1) = M(i:i+1,i:i+1)+ Me1;
    if dmodel=='v'
        C(i:i+1,i:i+1) = C(i:i+1,i:i+1)+ Ce1;
    end
    K(i:i+1,i:i+1) = K(i:i+1,i:i+1)+ Ke1;
end

% structural/histerisis damping model
if dmodel == 's'
    fm=flim/2;
    C=eta.*K/fm/2/pi;
    %C=0.*K;
    %K=K*(1+1i*eta);
end

gamma_c_fem=B*rho.*gamma_c;
% Inserting local feedback interactions 
M(n2,n1) = M(n2,n1)-gamma_c_fem(3);
K(n2,n1) = K(n2,n1)-gamma_c_fem(2);
C(n2,n1) = C(n2,n1)-gamma_c_fem(1);

% joining cells to make a structure
ncell = 18;
ndof = ncell*ne+1;

Ms=zeros(ndof,ndof);Cs=Ms;Ks=Ms;

for i=1:ne:ndof-ne
    Ms(i:i+ne,i:i+ne) = Ms(i:i+ne,i:i+ne) + M;
    Cs(i:i+ne,i:i+ne) = Cs(i:i+ne,i:i+ne) + C;
    Ks(i:i+ne,i:i+ne) = Ks(i:i+ne,i:i+ne) + K;
end

% Spatial vector x for FEM model 
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

%% Harmonic Forced Response 
% External force in unitary volume acceleration
% F=B*rho;
%-------------------------------------------------------------------------%
%                            FRF via FEM                                      % 
%-------------------------------------------------------------------------%
u_L_FEM = zeros(ndof,nf);u_R_FEM=u_L_FEM;
% frequency loop
for i = 1:nf
    w = wv(i);
    % Dinamic stiffness matrix via FEM
    Dg = -w^2*Ms + 1i*w*Cs + Ks;
    % External force in unitary volume velocity 
    % PS: (Ms,Ks,Cs) were previously multiplied by B*rho
    F=B*rho*1i*w;
    % force on left end
    F_L = [F; zeros(ndof-1,1)];
    % force on right end
    F_R = [zeros(ndof-1,1);F];
    aux = Dg\[F_L F_R];
    u_L_FEM(:,i) = aux(:,1);
    u_R_FEM(:,i) = aux(:,2);
end 
%%
%-------------------------------------------------------------------------%
%                            FRF via SEM                                  % 
%-------------------------------------------------------------------------%
% edit e to change where the force is applied 
e='l';
output = function_FRF_Ac(A1,A2,L1,L2,Lc,rho,c,feedback,dmodel,gamma_c,wv,ncell,e,eta);
P_SEM = output.ps_v; 

%%
%-------------------------------------------------------------------------%
%                     Comparing FEM and SEM models                        % 
%-------------------------------------------------------------------------%
figure
plot(fv/1000,20*log10(abs(P_SEM(end,:))),'k','MarkerSize',8,'LineWidth',1);hold on;
plot(fv/1000,20*log10(abs(u_L_FEM(end,:))),'r','MarkerSize',8,'LineWidth',1.2)
xlabel('Frequency [kHz]')
ylabel('FRF [dB]')
legend('SEM','FEM')
title('Left excitation and right measure')
set(gcf, 'Color', 'w');