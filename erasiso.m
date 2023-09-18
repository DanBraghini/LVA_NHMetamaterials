function [Aera,Bera,Cera,Dera]=erasiso(h)
%------------------------------------------------------------
%  Eigensystem Realization Algorithm
%  Estimacao de modelo de estados discreto
%  pelo metodo ERA de Juang e Pappa da NASA
%  Implementacao SISO Arruda 4/1992
% Agradecimentos J. Camino 06/2017
%  Entrada IRF's h vetor (1 x nt), nt >> (s+r)
%  Saída modelo discreto A,B,C,D   
%  Adaptação SIMO 01/2022
%------------------------------------------------------------

% Montar matriz de Hankel H(1):
[ndof, nt]=size(h);
fprintf('\n No. de time samples  nt = %2.0f ',nt);
fprintf('\n No. de DOFs medidos  ndof = %2.0f ',ndof);
s=input('\n No. de colunas da matriz de Hankel ( s+r+2 < nt ): ');
r=input('\n No. r de vezes que h e sobreposta na matriz de Hankel ( s+r+2 < nt ): ');
% Montar H1 e H2 a partir de h(0:nt-1)
%H=hankel(h(2:(r+s+2)));H1=H(1:r,1:s);H2=H(1:r,2:s+1);
H1=hankel(h(:,2:r+1) , h(r+1:r+s));
H2=hankel(h(3:r+2) , h(r+2:r+s+1));
% SVD de H1:
[U S V]=svd(H1);
% Selecionar SVs nao-nulos (desprezados):
figure(2)
plot(10*log10(abs(diag(S))))
title(' Escolha o No. de raizes com os SVs: ')
ylabel(' Singular Values ')
xlabel(' No. de raizes ')
nsvd=input(' No. de raizes: ');
% Truncar U,S e V com nsvd escolhido (nsvd=no de variáveis de estados, ordem do sistema):
St=S(1:nsvd,1:nsvd);
Ut=U(:,1:nsvd);
Vt=V(:,1:nsvd);
% Calcular system matrices Aera,Bera,Cera, Dera:
St2=sqrt(St);Sti2=inv(St2);
Aera=Sti2*Ut'*H2*Vt*Sti2;
Ob=Ut*St2;
Co=St2*Vt';
Cera=Ob(1:ndof,:); % (ndof x nsvd) 
Bera=Co(:,1); % (nsvd x 1) one input
Dera=h(1);