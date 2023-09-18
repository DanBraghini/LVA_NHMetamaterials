% UNICAMP - FEM - DMC
% Laboratório de Vibroacústica

% Danilo Beli 
% Prof. Dr. José Roberto de França Arruda

% 29/08/2013

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derivada Terceira de N transposta em função de s - delta3Nt_deltas3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [delta3Nt_deltas3] = delta3Nt_deltas3(gamma,s0,s)
% gama = matriz dos números de onda [1/m]
% s = posição do inicio do elemento [m]
% s0 = posição do fim do elemento [m]  

im = sqrt(-1); % número imaginário

delta3Nt_deltas3 = [((-im*gamma(1))^3)*exp(-im*gamma(1)*s), ...
                    ((-im*gamma(2))^3)*exp(-im*gamma(2)*s), ...
                    ((-im*gamma(3))^3)*exp(-im*gamma(3)*s), ...
                    ((-im*gamma(4))^3)*exp(-im*gamma(4)*(s-s0)), ...
                    ((-im*gamma(5))^3)*exp(-im*gamma(5)*(s-s0)), ...
                    ((-im*gamma(6))^3)*exp(-im*gamma(6)*(s-s0))];