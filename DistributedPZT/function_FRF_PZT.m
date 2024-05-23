function [u_L,u_M,u_R,u_passive,x,kDg] = function_FRF_PZT(varargin)
% by Danilo Braghini 16.08.2020
% modified at 08/2022
% Using the spectral element dynamic stiffness matrix for Piezoelectric
% Materials, this function return Forced Response in frequency domain
%% parametric inputs
param=struct(varargin{:});
E_p=param.Young_Modulus_PZT;
rho_p=param.density_PZT;
L_s=param.sensor_length;
L_a=param.actuator_length;
A_s=param.sensor_cross_area;
A_a=param.actuator_cross_area;
Gamma_c=param.feedback_matrix_actuator;
K_a=param.passive_matrix_actuator;
B_s=param.coeficient_sensor;
e_p=param.piezoelectric_constant;
alpha_p=param.dielectric_constant;
Lc=param.cell_length;
wv=param.frequency_vector;
F=param.impulse_amplitude;
ncell=param.number_cells;
a=param.non_locality;
%%
N=length(wv);
% Dynamic Stiffness Matrix  PZT rod with applied feedback: 
% Gamma_c
% Dynamic Stiffness Matrix  PZT rod with open-circuit (sensor)
mz = 3;
% number of spectral elements on each cell
ne = mz - 1;
Les = Lc/ne;
x = 0:Les:ncell*Lc;
ndof = ncell*ne+1;
% left excitation
% Fv(1) = F
m = round(ndof/2);
F_L = [F; zeros(ndof-1,1)]; % force on left end 
F_M = [zeros(ndof-m,1);F;zeros(ndof-m,1)]; % force in the middle
F_R = [zeros(ndof-1,1);F]; % force on right end
% This loop runs through the frequency vector
u_L=zeros(ndof,N);u_M=u_L;u_R=u_L;u_passive=u_L;
K_s = B_s.*[1 -1;-1 1];

% This loop runs through the frequency vector
 for n=1:N 
    w=wv(n);
    
    k_local_p = w*sqrt(rho_p/(E_p+e_p^2/alpha_p));
    % Dynamic Stiffness Matrix sensor PZT
    qsi_s = 1i*(E_p+e_p^2/alpha_p)*A_s*k_local_p;
    K_ps = (-1i*qsi_s/sin(k_local_p*L_s)).*[cos(k_local_p*L_s) -1; -1 cos(k_local_p*L_s)];   
    
    % Dynamic Stiffness Matrix actuator PZT
    qsi_a = 1i*(E_p+e_p^2/alpha_p)*A_a*k_local_p;
    K_pa = (-1i*qsi_a/sin(k_local_p*L_a)).*[cos(k_local_p*L_a) -1; -1 cos(k_local_p*L_a)];   

    % Passive part of Dynamic Stiffness Matrix  for actuator PZT 
    K_a_sc = K_pa + K_a;
 
    % combining Sensor PZT + Actuator PZT

    % Defining auxiliary variables
    D =  K_ps(2,2) + K_s(2,2) + K_a_sc(1,1);
    C =  K_ps(2,1) + K_s(2,1);
     
K_cell=[ K_ps(1,1) + K_s(1,1)   K_ps(1,2) + K_s(1,2)                      0                
              C                      D                          K_a_sc(1,2)
          0                         K_a_sc(2,1)                 K_a_sc(2,2)];

    % Defining auxiliary variables (K_g = 0)
    D =  K_ps(2,2) + K_s(2,2) + K_a_sc(1,1);
    C =  K_ps(2,1) + K_s(2,1) ;
     
K_cell_passive=[ K_ps(1,1) + K_s(1,1)   K_ps(1,2) + K_s(1,2)            0                
               C                      D                   K_a_sc(1,2)
               0                 K_a_sc(2,1)      K_a_sc(2,2)];
      
%   forced response using global dynamic stiffness matrix

    Dg=zeros(ndof,ndof);
    for i=(1:ne:ndof-ne)
        Dg(i:i+ne,i:i+ne)=Dg(i:i+ne,i:i+ne)+K_cell;
    end
    
 Gamma_a =   [0           0               0
              Gamma_c(1,1)  Gamma_c(1,2)  0
              Gamma_c(2,1)  Gamma_c(2,2)  0];
        
for i = a*ne+1:ne:ndof-ne
    lin = i:i+ne;
    col = i-a*ne:i-(a-1)*ne;
    Dg(lin,col) = Dg(lin,col) +  Gamma_a;
end

% for i = 1:ne:ndof-(a+1)*ne
%     Dg(a*ne+i:(a+1)*ne+i, i:ne+i) = Dg(a*ne+i:(a+1)*ne+i, i:ne+i) + Gamma_a;
% end

    %u(:,n) = Dg\Fv;
    u = Dg\[F_L F_M F_R];
    u_L(:,n) = u(:,1);
    u_M(:,n) = u(:,2);
    u_R(:,n) = u(:,3);
    kDg(n) =cond(Dg);
    Dg_passive=zeros(ndof,ndof);
    for i=(1:ne:ndof-ne)
        Dg_passive(i:i+ne,i:i+ne)=Dg_passive(i:i+ne,i:i+ne)+K_cell_passive;
    end
    u_passive(:,n) = Dg_passive\F_M; 
    
 end
end
