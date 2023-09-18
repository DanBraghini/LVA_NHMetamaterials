% UNICAMP - FEM - DMC
% Laboratório de Vibroacústica
%
% Danilo Beli 
% Priscilla Brandao Silva
% Prof. Dr. José Roberto de França Arruda
%
% 29/08/2013
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Curved Beam
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [D_elem_flexiblering2,gamma] = D_elem_flexiblering2(R,b,h,rho,E,s0,omega)

% R = raio [m]
% b = largura da seção [m]
% h = altura da seção [m]
% rho = densidade do material [kg/m^3]
% E = módulo de elasticidade do material [Pa]
% s0 = comprimento angular do elemento [rad]
% omega = velocidade de propagação da onda [rad/s]
%
im = sqrt(-1); % número imaginário
%
A = b*h; % área da seção transversal [m^2]
I = b*h^3/12; % segundo momento de área [m^4]
%
% Equação característica / Relação de dispersão
p6 = -E*A*E*I*R^2;
p5 = 0;
p4 = 2*E*A*E*I+rho*A*E*I*R^2*omega^2;
p3 = 0;
p2 = -E*I*E*A/R^2+E*I*rho*A*omega^2+E*A*rho*A*R^2*omega^2;
p1 = 0;
p0 =+E*A*rho*A*omega^2-rho^2*A^2*R^2*omega^4;
%
eq_caracteristica = [p6 p5 p4 p3 p2 p1 p0];
%
% Raízes da equação característica
rr = roots(eq_caracteristica);
%
% Separar os números de onda - positivos e negativos
ar = 1; br = 1;
r_prop_positivo=[];
for cont2 = 1:1:length(rr)
    if imag(rr(cont2))< 0 
        r_prop_positivo(ar,1) = rr(cont2); ar=ar+1; % números de onda que se propagam no sentido positivo
    else
        r_prop_negativo(br,1) = rr(cont2); br=br+1; % números de onda que se propagam no sentido negativo
    end
end
%%
% a1=p4/p6;
% a2=p2/p6;
% a3=p0/p6;
% Q= (3*a2-a1^2)/9;
% R=(9*a1*a2 -27*a3-2*a1^3)/54;
% S1=(R+(Q^3+R^2)^(1/2))^(1/3);
% S2=(R-(Q^3+R^2)^(1/2))^(1/3);
% r1=S1+S2-a1/3;
% r2=-(S1+S2)/2-a1/3+im*sqrt(3)*(S1-S2)/2;
% r3=-(S1+S2)/2-a1/3-im*sqrt(3)*(S1-S2)/2;
% 
% gamma1= sqrt(r1);gamma2= -gamma1; 
% gamma3= sqrt(r2);gamma4= -gamma3;
% gamma5= sqrt(r3);gamma6= -gamma5;
% gamma = [gamma1;gamma2;gamma3;gamma4;gamma5;gamma6];
%%
% Ordenando os números de onda -
% Válido para velocidade de rotação (Omega0) > 0

r_pos = sort(r_prop_positivo,'ComparisonMethod','real'); 
r_neg = sort(r_prop_negativo,'ComparisonMethod','real');
%
gamma = [r_pos;r_neg];
%
% Razão de amplitude
alpha = im*(E*I.*(gamma.^3)+E*A.*gamma)./ ...
                ((-E*I/R).*(gamma.^2)-E*A*R.*(gamma.^2)+rho*A*R*omega^2);
%
% Matriz que relaciona os deslocamentos com as constantes: U = G*C
G = [        alpha(1)                                       alpha(2)                                    alpha(3)                      (alpha(4))*exp(im*(gamma(4))*s0)      (alpha(5))*exp(im*(gamma(5))*s0)         (alpha(6))*exp(im*(gamma(6))*s0)
                1                                               1                                           1                               exp(im*(gamma(4))*s0)                  exp(im*(gamma(5))*s0)                   exp(im*(gamma(6))*s0)
            -im*gamma(1)                                    -im*gamma(2)                                -im*gamma(3)                -im*(gamma(4))*exp(im*(gamma(4))*s0)    -im*(gamma(5))*exp(im*(gamma(5))*s0)    -im*(gamma(6))*exp(im*(gamma(6))*s0)
     alpha(1)*exp(-im*(gamma(1))*s0)                alpha(2)*exp(-im*(gamma(2))*s0)          alpha(3)*exp(-im*(gamma(3))*s0)                      alpha(4)                              alpha(5)                                   alpha(6)
        exp(-im*(gamma(1))*s0)                           exp(-im*(gamma(2))*s0)                    exp(-im*(gamma(3))*s0)                            1                                     1                                           1
     -im*(gamma(1))*exp(-im*(gamma(1))*s0)      -im*(gamma(2))*exp(-im*(gamma(2))*s0)       -im*(gamma(3))*exp(-im*(gamma(3))*s0)              -im*gamma(4)                            -im*gamma(5)                              -im*gamma(6)];
%
% Matriz da razão de amplitude 
Alpha = [alpha(1) 0 0 0 0 0
         0 alpha(2) 0 0 0 0
         0 0 alpha(3) 0 0 0
         0 0 0 alpha(4) 0 0
         0 0 0 0 alpha(5) 0
         0 0 0 0 0 alpha(6)];

% Matriz que relaciona as forças com as constantes: F = H*C
H = [-E*A*  (deltaNt_deltas(gamma,s0,0)*Alpha              +(1/R)*Nt_s(gamma,s0,0))         
     -E*I*   ((1/R)*delta2Nt_deltas2(gamma,s0,0)*Alpha      -delta3Nt_deltas3(gamma,s0,0))
     E*I*  ((1/R)*deltaNt_deltas(gamma,s0,0)*Alpha             -delta2Nt_deltas2(gamma,s0,0))
     E*A*   (deltaNt_deltas(gamma,s0,s0)*Alpha             +(1/R)*Nt_s(gamma,s0,s0))        
     E*I*  ((1/R)*delta2Nt_deltas2(gamma,s0,s0)*Alpha     -delta3Nt_deltas3(gamma,s0,s0))
     -E*I*   ((1/R)*deltaNt_deltas(gamma,s0,s0)*Alpha       -delta2Nt_deltas2(gamma,s0,s0))];
%
 % Matriz Espectral de Rigidez Dinâmica: F = H*(inv(G))*U -> K = H*(inv(G))
 warning off all;
D_elem_flexiblering2 = H*(inv(G));
 warning on all;