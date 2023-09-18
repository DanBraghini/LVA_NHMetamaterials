function [kL_sem_PB,kL_sem_SB] = function_SEM_Acoustic(L1,L2,rho,A1,A2,c,wv,eta,varargin)
%%
if nargin > 8
    if nargin == 9
        options = varargin{1};
    else
        options = struct(varargin{:});
    end
else
    options = [];
end
if ~isfield(options,'feedback')
    options.feedback = 'pid';
end
if ~isfield(options,'gain')
   if options.feedback == 'pid'
        options.gain = [0 0 0];
   else
       options.gain=0;
   end
end

if ~isfield(options,'dmodel')
    options.dmodel = 's';
end

%%

gamma_c=options.gain;
feedbacktype=options.feedback;
dmodel=options.dmodel;
% Created by Danilo Braghini 05.06.2020
% Spectral element dynamic stiffness matrix for Rods
N=length(wv);
kL_sem_PB=zeros(1,N);kL_sem_SB=kL_sem_PB;
 % This loop runs through the frequency vector
 for n=1:N 
     w=wv(n);
    %% select the damping model
    % viscous damping
    if dmodel == 'v'
        B=rho*c^2;
        k_local=sqrt((w^2*rho - 1i*w*eta)/B);
    % structural damping (histeresis)
    elseif dmodel == 's' 
        k_local= w/c;
        k_local=k_local*(1-1i*eta/2);
    else 
        disp('invalidy damping model')
        return
    end
    %% Segment 1 and 3        
    k11=-1-exp(-2*1i*k_local*L1);
        k12= 2*exp(-1i*k_local*L1);
        k21=k12;
        k22=k11;
    Km1=[ k11   k12
          k21   k22 ].*(A1/(rho*c*(exp(-2*1i*k_local*L1)-1)));

    k11=-1-exp(-2*1i*k_local*L2);
        k12= 2*exp(-1i*k_local*L2);
        k21=k12;
        k22=k11;

    Km2=[ k11   k12
         k21   k22 ].*(A2/(rho*c*(exp(-2*1i*k_local*L2)-1)));
     % auxiliary variables
     C1=Km1(2,2)+Km2(1,1);
     %% select feedback law type
     if feedbacktype=='pid'
         %proportional
         gain(1)=gamma_c(1);
         %integral
         gain(2)=gamma_c(2)/1i/w;
         %derivative
         gain(3)=gamma_c(3)*1i*w;
         C2=Km2(2,1)-norm(gain,1);

     elseif feedbacktype=='r'
        %derivative aproximation
        Nd=1000*wv(end);
        gain=gamma_c*Nd/(1+Nd/(1i*w));
        C2=Km2(2,1)-gain;

%      else
%          gain=gamma_c/1i/w.*function_hwindow(w,114.3*2*pi,415.8*2*pi);
     end
     %% build cell dynamic stiffness matrix
     % combining L1 and L2 and L1
    Kc=[Km1(1,1)   Km1(1,2)       0                 0
         Km1(2,1)    C1          Km2(1,2)           0
            0        C2           C1                Km1(1,2)                    
            0        0          Km1(2,1)          Km1(2,2)];

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
    %% transfer matrix
    Tc=[-Kr(1,1)/Kr(1,2)                                  1/Kr(1,2)
          -Kr(2,1)+Kr(2,2)*Kr(1,1)/Kr(1,2)       -Kr(2,2)/Kr(1,2)];
    %% solving eigenvalue problem
    [~,Lambda]=eig(Tc);
    Fm1=Lambda(1,1);
    mu=-1i*log(Fm1);
    %% resultant normalized wavenumber kL
    kL_sem_PB(n)=real(mu);
    kL_sem_SB(n)=imag(mu);
   
 end
 
end
