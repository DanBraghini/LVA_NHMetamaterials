
% geracao da matriz h para aplicacao do metodo POLYREFERENCE
% NOTA: Para rodar ITD2 faca Ni = 1
% sistema de 3 DOF com freqs. naturais da ordem de 100 rad/s
% Arruda 28/11/89

fprintf(' ATENCAO!! memoria sera apagada antes de iniciar!\n aperte qqer tecla p/ continuar...');
pause

% limpando a memoria
clear
% ordem do sistema
N=3;
% no. de inputs
Ni=input('No. de referencias=');
if Ni==1,
   Ne=input('No. do pto. de excitacao=');
end
% no. de pontos no tempo para cada linha de h
M=input('No. de pontos do vetor de tempo M=');
% no. de pontos de medida da resposta
Np=3;
% resolucao no tempo
Delta=input('Delta t (0.01) =');
% porcentagem de ruido
Porc=input('Porcentagem de ruido=');
Porc=Porc/100.;

% parametros do sistema
m=[2 0 0
   0 1 0
   0 0 1];
k=[10  -5    0
   -5  15  -10
    0 -10   20];
k=1000*k;
c=[ 25 -10   0
   -10  20 -10
     0 -10  15];

% autovalores e autovetores
O=zeros(size(k));
A=[c,m;m,O];
B=[k,O;O,-m];
% -B Psi_r = A Psi_r s_r
[Psi,s]=eig(-B,A);
nx=size(A,1);
%normalizacao
norm=conj(Psi')*A*Psi;
Psi=Psi*inv(sqrt(norm));
for ii=1:nx
  exps(ii)=exp(s(ii,ii)*Delta);
end

% monta matriz h
h=zeros(Ni*Np,M);
if Ni>1
    for ii=1:Np
      for j=1:M
         for r=1:nx
           h(Ni*(ii-1)+1:Ni*ii,j)=h(Ni*(ii-1)+1:Ni*ii,j)+real(Psi(ii,r)*exps(r)^(j-1)*Psi(1:Ni,r));
         end
      end
      ii
    end
elseif Ni==1
    for ii=1:Np
      for j=1:M
         for r=1:2*N
           h(ii,j)=h(ii,j)+real(Psi(ii,r)*exps(r)^(j-1)*Psi(Ne,r));
         end
      end
    end
end


% poluicao de h com ruido
hrms=zeros(1,Ni*Np);
for i=1:Ni*Np
   for j=1:M
     hrms(i)=hrms(i)+h(i,j)^2;
   end
   hrms(i)=sqrt(hrms(i)/M);
   h(i,:)=h(i,:)+Porc*hrms(i)*randn(1,M);
end

% calculo das FRFs 
H=fft(h.');
H=H(1:M/2,:);
w=2*pi*((1:M/2)-1)/(M*Delta);