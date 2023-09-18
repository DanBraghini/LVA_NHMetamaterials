% UNICAMP - FEM - DMC
% Laboratório de Vibroacústica

% Danilo Beli 
% Prof. Dr. José Roberto de França Arruda

% 29/08/2013

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derivada de N transposta em função de s - deltaNt_deltas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [deltaNt_deltas] = deltaNt_deltas(gamma,s0,s)
% gama = matriz dos números de onda [1/m]
% s = posição do inicio do elemento [m]
% s0 = posição do fim do elemento [m]  

im = sqrt(-1); % número imaginário

deltaNt_deltas = [((-im*gamma(1))^1)*exp(-im*gamma(1)*s), ...
                  ((-im*gamma(2))^1)*exp(-im*gamma(2)*s), ...
                  ((-im*gamma(3))^1)*exp(-im*gamma(3)*s), ...
                  ((-im*gamma(4))^1)*exp(-im*gamma(4)*(s-s0)), ...
                  ((-im*gamma(5))^1)*exp(-im*gamma(5)*(s-s0)), ...
                  ((-im*gamma(6))^1)*exp(-im*gamma(6)*(s-s0))];