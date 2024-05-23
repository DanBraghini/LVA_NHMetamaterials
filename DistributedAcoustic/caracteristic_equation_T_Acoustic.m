%by Danilo Braghini 16.08.2020
function f=caracteristic_equation_T_Acoustic(w)
% Function used to compute the dispersion relation indirectly w(k).
% This function writes down the caracteristic equation vof the transfer matrix
% using SEM model for active acoustic duct(homogenious duct) with local
% feedback.
% As a function of frequency w, corresponding to a fixed wavenumber k,i.e, w(k) the equation is transcendental, since the
% matrix has trigonometric or hiperbolic functions (exponentials). Thus,
% this non-linear can be calling the function fsolve on the main code.
% Inputs( given as global variables):
% dmodel : string with the model chosen.
%     's': structural;
%     'v': viscousrho;
%     'n': none;
% L1, L2: lengths of the segments of the unit cell L1-L2-L1;
% A1, A2: cross sectional areas of the unit cell A1-A2-A1;
% kL: value of wavenumber times length. kL \in R | [-pb*pi, pb*pi];
% for pb Brillouin Zones;
% c: sound velocity within the dute;
% eta: damping parameter;
% gamma_c: feedback gain. [k_p k_i k_d] for PID;
% H_pv : transfer fucntion in pressure by volume velocity;
% ideal_filter: 1 if ideal filter is choosen. 0 otherwise.
global  dmodel rho L1 L2 A1 A2 kL c eta gamma_c H_pa ideal_filter

 B=rho*c^2;
%% select the damping model
B=rho*c^2;
% viscous damping
if dmodel == 'v'
    s=tf('s');
    H_p=H_pa/s;
    k_local=sqrt((w^2*rho - 1i*w*eta)/B);

    k11=1+exp(-2*1i*k_local*L1);
        k12= -2*exp(-1i*k_local*L1);
        k21=k12;
        k22=k11;
    Kd1=[ k11   k12
               k21   k22 ].*(A1/(rho*c*(1-exp(-2*1i*k_local*L1))));

    k11=1+exp(-2*1i*k_local*L2);
        k12= -2*exp(-1i*k_local*L2);
        k21=k12;
        k22=k11;

    Kd2=[ k11   k12
              k21   k22 ].*(A2/(rho*c*(1-exp(-2*1i*k_local*L2))));
% structural damping (histeresis)
elseif dmodel == 's'  || dmodel == 'n'
    H_p=H_pa;
    k_local= w/c;
    k_local=k_local*(1-1i*eta/2);
    % acceleration
    Am1=[0                                      -rho/A1;
          A1*k_local^2/rho       0];
    %velo
    %Am1=[0                                         -1j*k_local*sqrt(B*rho)/A1
    %        -1j*k_local*A1/sqrt(rho*B)          0];
    T=expm(Am1*L1);
    Kd1=-[T(1,1)/T(1,2)                                -1/T(1,2)
          T(2,1)-T(2,2)*T(1,1)/T(1,2)         T(2,2)/T(1,2)];
    % acceleration
     Am2=[0                          -rho/A2;
     A2*k_local^2/rho         0];
    %velocity
    %Am2=[0                                         -1j*k_local*sqrt(B*rho)/A2
    %        -1j*k_local*A2/sqrt(rho*B)          0];
     T=expm(Am2*L2);
     Kd2=-[T(1,1)/T(1,2)                                -1/T(1,2)
          T(2,1)-T(2,2)*T(1,1)/T(1,2)         T(2,2)/T(1,2)];
else 
    disp('invalidy damping model')
    return
end

%% build dynamic stiffness matrix on pressure by volume velocity
% k11=-1-exp(-2*1i*k_local*L1);
% k12= 2*exp(-1i*k_local*L1);
% k21=k12;
% k22=k11;
% Kd1=[ k11   k12
%       k21   k22 ].*(A1/(rho*c*(exp(-2*1i*k_local*L1)-1)));
% 
%               
% k11=-1-exp(-2*1i*k_local*L2);
% k12= 2*exp(-1i*k_local*L2);
% k21=k12;
% k22=k11;
% 
% Kd2=[ k11   k12
%      k21   k22 ].*(A2/(rho*c*(exp(-2*1i*k_local*L2)-1)));
% auxiliary variables
 % feedback gain is given in terms of applied volume velocity
 C1=Kd1(2,2)+Kd2(1,1);
%% adding feedback law
gain=evalfr(H_p,1i*w);
if ideal_filter == 1
    fase=0;
    gain=gain*function_idealfilter(w,2*pi*0.5e3,2*pi*2e3,fase);
end
%*function_hwindow(w,2*pi*0.5e3,2*pi*1e3);
C2=Kd2(2,1)-gamma_c*gain;
 %% build cell matrix
 % combining L1 and L2 and L1
Kc=[Kd1(1,1)   Kd1(1,2)       0                 0
     Kd1(2,1)    C1          Kd2(1,2)           0
        0        C2           C1                Kd1(1,2)                    
        0        0          Kd1(2,1)          Kd1(2,2)];

%% condensation
D=inv(Kc(2:3,2:3));
Kll = Kc(1,1);
Krr = Kc(4,4);
Kuu = Kc(1,2:3);
Klc = Kc(2:3,1);
Klr = Kc(1,4); %=0
Krc = Kc(2:3,4);
Krl = Kc(4,1); % = 0
Kbb = Kc(4,2:3);


Kr=[Kll-Kuu*D*Klc              Klr-Kuu*D*Krc
      Krl-Kbb*D*Klc             Krr-Kbb*D*Krc];
%% Transfer matrix
Tc=[-Kr(1,1)/Kr(1,2)                                  1/Kr(1,2)
      -Kr(2,1)+Kr(2,2)*Kr(1,1)/Kr(1,2)       -Kr(2,2)/Kr(1,2)];
 

% Kc11 = K1(1,1) - K1(1,2)*C1/(D1*D3);   
% Kc12 = K1(1,2)*K2(1,2)*K3(1,2)/(D1*D2*D3);
% Kc21 = C1*C2*K3(2,1)/(D1*D2*D3);
% Kc22 = K3(2,2) - K3(2,1)*K3(1,2)/(D2*D3);
% 
% Tc(1,1) = -Kc22/Kc12;
% Tc(1,2) = 1/Kc12;
% Tc(2,1) = Kc11*Kc22/Kc12 -Kc21;    
% Tc(2,2) = -Kc22/Kc12;

%% caracteristic equation
f=det(Tc-exp(-1i*kL).*eye(2));

end


