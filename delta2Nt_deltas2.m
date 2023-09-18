% UNICAMP - FEM - DMC
% Laborat�rio de Vibroac�stica

% Danilo Beli 
% Prof. Dr. Jos� Roberto de Fran�a Arruda

% 29/08/2013

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derivada Segunda de N transposta em fun��o de s - delta2Nt_deltas2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [delta2Nt_deltas2] = delta2Nt_deltas2(gamma,s0,s)
% gama = matriz dos n�meros de onda [1/m]
% s = posi��o do inicio do elemento [m]
% s0 = posi��o do fim do elemento [m]  

im = sqrt(-1); % n�mero imagin�rio
               
delta2Nt_deltas2 = [((-im*gamma(1))^2)*exp(-im*gamma(1)*s), ...
                    ((-im*gamma(2))^2)*exp(-im*gamma(2)*s), ...
                    ((-im*gamma(3))^2)*exp(-im*gamma(3)*s), ...
                    ((-im*gamma(4))^2)*exp(-im*gamma(4)*(s-s0)), ...
                    ((-im*gamma(5))^2)*exp(-im*gamma(5)*(s-s0)), ...
                    ((-im*gamma(6))^2)*exp(-im*gamma(6)*(s-s0))];