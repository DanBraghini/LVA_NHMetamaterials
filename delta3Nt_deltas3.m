% UNICAMP - FEM - DMC
% Laborat�rio de Vibroac�stica

% Danilo Beli 
% Prof. Dr. Jos� Roberto de Fran�a Arruda

% 29/08/2013

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derivada Terceira de N transposta em fun��o de s - delta3Nt_deltas3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [delta3Nt_deltas3] = delta3Nt_deltas3(gamma,s0,s)
% gama = matriz dos n�meros de onda [1/m]
% s = posi��o do inicio do elemento [m]
% s0 = posi��o do fim do elemento [m]  

im = sqrt(-1); % n�mero imagin�rio

delta3Nt_deltas3 = [((-im*gamma(1))^3)*exp(-im*gamma(1)*s), ...
                    ((-im*gamma(2))^3)*exp(-im*gamma(2)*s), ...
                    ((-im*gamma(3))^3)*exp(-im*gamma(3)*s), ...
                    ((-im*gamma(4))^3)*exp(-im*gamma(4)*(s-s0)), ...
                    ((-im*gamma(5))^3)*exp(-im*gamma(5)*(s-s0)), ...
                    ((-im*gamma(6))^3)*exp(-im*gamma(6)*(s-s0))];