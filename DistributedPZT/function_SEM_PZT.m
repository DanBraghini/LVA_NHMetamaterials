function [kL_sem_PB,kL_sem_SB] = function_SEM_PZT(Gamma_c,K_a,B_s,L_s,L_a,E_p,A_s,A_a,rho_p,e_p,alpha_p,wv);
% by Danilo Braghini 16.08.2020
% Using the spectral element dynamic stiffness matrix for Piezoelectric
% Rods, this function implements the direct method to calculate both passage and stop
% bands of the dispersion relation
% inputs:
%  Gamma_c = Dynamic Stiffness Matrix  PZT rod with applied feedback: 
%  B_s = Coeficient related to sensor PZT
%  L_s = Sensor's length
%  L_a = Acuator's length
%  E_p = Young's modulos of piezoelectric material
%  A_s = Sensor's ross section area
%  A_a = Actuator's ross section area
%  rho_p = Volumetric density of piezoelectric material
%  e_p = piezoelectric constant of piezoelectric material
%  alpha = dielectric constant of piezoelectric material
%  wv = angular frequency vector
% outputs (simulated with SEM formulation):
% kL_sem_PB = wavenumber x cell's length related to the passage bands
% kL_sem_SB = wavenumber x cell's length related to the stop bands

N=length(wv);
% Dynamic Stiffness Matrix  for sensor PZT 
 K_s = B_s.*[1 -1;-1 1];

 
 kL_sem_PB = zeros(N,1);kL_sem_SB=kL_sem_PB;   
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
    K_passive_PZT = K_pa + K_a;
 
    % combining Sensor PZT + Actuator PZT

    % Defining auxiliary variables
    D =  K_ps(2,2) + K_s(2,2) + K_passive_PZT(1,1) + Gamma_c(1,2);
    C =  K_ps(2,1) + K_s(2,1) +  Gamma_c(1,1);
     
    K_cell=[ K_ps(1,1) + K_s(1,1)   K_ps(1,2) + K_s(1,2)                      0                
                  C                      D                          K_passive_PZT(1,2)
              Gamma_c(2,1)         K_passive_PZT(2,1)+Gamma_c(2,2)  K_passive_PZT(2,2)];
              
% condensation
    np = length(K_cell);
%     ni = np -2;

    Kcc = K_cell(2:np-1,2:np-1);
    Klu = K_cell(1,1);
    Krb = K_cell(np,np);
    Kuu = K_cell(1,2:np-1);
    Klc = K_cell(2:np-1,1);
    Kru = K_cell(1,np); %=0
    Krc = K_cell(2:np-1,np);
    Klb = K_cell(np,1); % = 0
    Kbb = K_cell(np,2:np-1);
    Z1 = Kcc\Klc;
    Z2 = Kcc\Krc;

    Kr=[Klu-Kuu*Z1              Kru-Kuu*Z2
         Klb-Kbb*Z1            Krb-Kbb*Z2];

    Tc=[-Kr(1,1)/Kr(1,2)                                  1/Kr(1,2)
          -Kr(2,1)+Kr(2,2)*Kr(1,1)/Kr(1,2)       -Kr(2,2)/Kr(1,2)];
     
    [~,Lambda]=eig(Tc);
    % Floquet multiplier 1
    Fm1 = Lambda(1,1);
    mu = -1i*log(Fm1);
    %normalized wavenumber kL
    kL_sem_PB(n)=real(mu);
    kL_sem_SB(n)=imag(mu);
 end
 
end
