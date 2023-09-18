function [Aera,Bera,Cera,Dera]=erapoly(h)
%------------------------------------------------------------
%  Eigensystem Realization Algorithm
%  Estimacao de parametros modais pelo metodo ERA de
%  Juang e Pappa, NASA Implementacao Arruda 4/1992 rev. 10/97
%   NOTA: Rodar antes "gerah" para gerar matriz com IRF's h.
%  erapoly.m
%------------------------------------------------------------
%         Resultados em: EFreq, EDamp, EPsi
%   IMPLEMENTACAO POLIREFERENCIA
%   fazendo  Y(k) (mxnr) com  nr = No. de referencias
%                              m = No. de amostras temporais
%------------------------------------------------------------

% Montar matriz de Hankel H(1):
H1=0;
H2=0;
[ntot nt]=size(h);
fprintf('\n No. de time samples  nt = %2.0f ',nt);
fprintf('\n No. total de linhas de h ntot = %2.0f  ',ntot);
% inputs
nr=input('No. de referências: ');
% measures:
%ndof=ntot/nr;
ndof=input('No. de DOFs medidos: ');
if ntot ~= nr*ndof,
  fprintf('\n No. de refs e de DOFs não compatíveis !!! \n');
  break;
end

% "arrumando" h para o ERA Polireferência:

hr=[];
for j=1:nt,
  hj=[];
  for i=1:ndof,
    hj=[hj
        h((i-1)*nr+1:i*nr,j)'];
  end
  hr=[hr hj];
end


s=input('No. s para s*nr colunas da matriz de Hankel ( s+r+2 < nt ): ');
r=input('No. r de vezes que h e sobreposta na matriz de Hankel: ');

s=s-1;
for j=1:s+1
  im=1;
  ir=0;
  for i=1:r*ndof
    H1(i,(j-1)*nr+1:j*nr)= hr(im,(j+ir-1)*nr+1:(j+ir)*nr);
    im=im+1;
    if im > ndof
       im=1;
       ir=ir+1;
    end
  end
end

% SVD de H1:
[U, S, V]=svd(H1);

% Selecionar SVs nao-nulos:
plot(diag(S))
title(' Escolha o No. de raizes com os SVs: ')
ylabel(' Singular Values ')
xlabel(' No. de raizes ')

nsvd=input(' No. de raizes: ');

% Truncar U,S e V com nsvd escolhido:
St=S(1:nsvd,1:nsvd);
Ut=U(:,1:nsvd);
Vt=V(:,1:nsvd);

% Montar H(2):
for j=1:s+1
  im=1;
  ir=0;
  for i=1:r*ndof
    H2(i,(j-1)*nr+1:j*nr)= hr(im,(j+ir)*nr+1:(j+ir+1)*nr);
    im=im+1;
    if im > ndof
       im=1;
       ir=ir+1;
    end
  end
end

% Calcular system realization matrices:
St2=sqrt(St);
Sti2=inv(St2);
Enr=[eye(nr) zeros(nr,s*nr)]';
Endof=[eye(ndof) zeros(ndof,(r-1)*ndof)];

Aera=Sti2*Ut'*H2*Vt*Sti2;
Bera=St2*Vt'*Enr;
Cera=Endof*Ut*St2;


% % Calcular eigensolutions para As:
% ss=0;
% Freq=0;
% Damp=0;
% [Fi,Lamb]=eig(Aera);
% for i=1:nsvd
%   ss(i)=log(Lamb(i,i))/Delta;
%   wd(i)=imag(ss(i));
%   wr(i)=abs(ss(i));
%   Freq(i)=wr(i)/(2*pi);
%   Damp(i)=real(ss(i))/wr(i);
% end
% 
% % Transformando Psi para coordenadas fisicas:
% Ct=Ut*St2;
% Psi=Ct*Fi;
% 
% % calculo do  MCF:
% jm=sqrt(-1);
% Q=inv(Fi)*St2*Vt';
% Q=Q(:,1:nr:(s+1)*nr);
% Qe=Q(:,1);
% for i=1:s,
%    esj=Lamb^i;
%    Qe=[Qe esj*Q(:,1)];
% end
% clear esj;
% for i=1:nsvd,
%    MCF(i)=(Qe(i,:)*Q(i,:)')^2/((Qe(i,:)*Qe(i,:)')*(Q(i,:)*Q(i,:)'));
%    MCF(i)=sqrt(abs(MCF(i)));
% end
% 
% 
% % eliminando os pares conjugados para impressao dos resultados
% ij=0;
% nm=0;
% for i=1:nsvd,
%    if wd(i)<0.,
%      ij=ij+1;
%      nm(ij)=i;
%    end
% end
% % nm(i) ; i=1,ij. guarda localizacao doa autovalores p/ impressao
% 
% % normalizacao dos modos :
% for i=1:ij,
%   Psi(:,nm(i))=Psi(:,nm(i))/Psi(1,nm(i));
% end
% 
% % resultados finais para comparacao entre metodos :
% EPsi=zeros(r*ndof,ij);EDamp=zeros(ij);EFreq=EDamp;EMCF=EDamp;
% for i=1:ij,
%   EFreq(i)=abs(Freq(nm(i)));
%   EDamp(i)=abs(Damp(nm(i)));
%   EPsi(:,i)=Psi(:,nm(i));
%   EMCF(i)=MCF(nm(i));
% end
% 
% % ordenamento dos modos:
% for il=ij:-1:2,
%   for i=ij-il+2:ij,
%     if EFreq(ij-il+1) > EFreq(i),
%         Buf=EFreq(ij-il+1);
%         Bud=EDamp(ij-il+1);
%         Bup=EPsi(:,ij-il+1);
%         Bum=EMCF(ij-il+1);
%         EFreq(ij-il+1)=EFreq(i);
%         EDamp(ij-il+1)=EDamp(i);
%         EPsi(:,ij-il+1)=EPsi(:,i);
%         EMCF(ij-il+1)=EMCF(i);
%         EFreq(i)=Buf;
%         EDamp(i)=Bud;
%         EPsi(:,i)=Bup;
%         EMCF(i)=Bum;
%     end
%   end
% end
% 
% 
% % reimpressao de resultados de polyref e polymode
% % impressao dos resultados
% fprintf('\n modo no.    freq. natural  amortecimento    MCF \n\n');
% for i=1:ij,
%   p1=EFreq(i);
%   p2=EDamp(i);
%   p3=EMCF(i);
%   fprintf('   %2.0f  ',i);
%   fprintf('      %8.3g Hz       %5.3f         %5.3f\n',p1,p2,p3);
% end
% 
% 
% 
% 
% % impressao dos resultados
% im=1;
% while im>0,
%   im=input('entre o no. do modo desejado ou 0 para sair:');
%   if im>0,
%     fprintf('\n modo no.: %2.0f  frequencia: %8.3f Hz\n',im,abs(EFreq(im)));
%     fprintf('  GDL no.         Amplitude do modo\n');
%     for i=1:ndof,
%       fprintf('    %2.0f         %8.3g %8.3g i\n',i,real(EPsi(i,im)),imag(EPsi(i,im)));
%     end
%   end
% end