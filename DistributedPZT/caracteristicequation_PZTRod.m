function f=caracteristicequation_PZTRod(w,varargin)
%by Danilo Braghini updated at 08/2022. Function used to compute the
%caracterisc equation of the transfer matrix for the cell of a PZT rod with
%feedback
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
    B1=param.coeficient_sensor;
    e_p=param.piezoelectric_constant;
    alpha_p=param.dielectric_constant;
    kL_global=param.wavenumber;
 %%
    %c_s=sqrt((E_p+e_p^2/alpha_p)/rho_p);
    %ci_a=sqrt((E_p+e_p^2/alpha_p)/rho_p);
    k_local_p = w*sqrt(rho_p/(E_p+e_p^2/alpha_p));
    %% Dynamic Stiffness Matrix sensor PZT
    qsi_s = 1i*(E_p+e_p^2/alpha_p)*A_s*k_local_p;
    K_ps = (-1i*qsi_s/sin(k_local_p*L_s)).*[cos(k_local_p*L_s) -1; -1 cos(k_local_p*L_s)];   
    
    K_s = B1.*[1 -1;-1 1];
    %% Dynamic Stiffness Matrix actuator PZT
    qsi_a = 1i*(E_p+e_p^2/alpha_p)*A_a*k_local_p;
    K_pa = (-1i*qsi_a/sin(k_local_p*L_a)).*[cos(k_local_p*L_a) -1; -1 cos(k_local_p*L_a)];   
    % Passive part of Dynamic Stiffness Matrix  for actuator PZT 
    K_passive_PZT = K_pa + K_a;
    %% combining Sensor PZT + Actuator PZT
    % Defining auxiliary variables
    D =  K_ps(2,2) + K_s(2,2) + K_passive_PZT(1,1) + Gamma_c(1,2);
    C =  K_ps(2,1) + K_s(2,1) +  Gamma_c(1,1);
     
    K_cell=[ K_ps(1,1) + K_s(1,1)   K_ps(1,2) + K_s(1,2)                      0                
              C                      D                          K_passive_PZT(1,2)
          Gamma_c(2,1)         K_passive_PZT(2,1)+Gamma_c(2,2)  K_passive_PZT(2,2)];
    %% condensation
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
    %% Transfer matriz
    Tc=[-Kr(1,1)/Kr(1,2)                                  1/Kr(1,2)
          -Kr(2,1)+Kr(2,2)*Kr(1,1)/Kr(1,2)       -Kr(2,2)/Kr(1,2)];
    %% caracteristic equation for given wavenumber
    f=det(Tc-exp(-1i*kL_global).*eye(2));

end


